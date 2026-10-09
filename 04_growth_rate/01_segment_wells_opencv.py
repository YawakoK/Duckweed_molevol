# Step 01: plate cropping, well splitting and primary frond mask (OpenCV, CIELAB).
#
# For every plate photo (Olympus TG-6, one photo per plate per day, Day0-Day7):
#   1. the plate is located with four ArUco markers (DICT_4X4_50, ids 0-3) and
#      perspective-corrected to 3000 x 2000 px;
#   2. the plate is split into its 2 x 3 wells;
#   3. the primary green-frond mask is thresholded in 8-bit CIELAB
#      (L 0-255, a 90-135, b 155-255) and thin line / crescent-shaped
#      components (reflections) are removed;
#   4. frond area = mask pixels x 0.0016 mm2.
# The cropped well images and binary masks are reused by step 02.
#
# Input : photos/<date>/Day<N>/<PhotoID>.JPG     raw photos (Zenodo dataset; not in this repository)
#         data/image_metadata/photoinfo_<date>.csv   EXIF capture time + plate ID per photo
#         data/image_metadata/sample_list_<date>.csv plate/well -> species, light
# Output: processed/{cropped_plates,marked_aruco,split_wells,binary_wells}/<date>/*.jpg
#         data/frond_area_main_mask_by_well_image.csv   one row per well image (all species)
#         output/image_processing_failures.csv
# Run from this directory:  python 01_segment_wells_opencv.py
# Author: Natsu Katayama
from __future__ import annotations

import re
from pathlib import Path

import cv2
import numpy as np
import pandas as pd


PHOTO_ROOT = Path("photos")
METADATA_DIR = Path("data") / "image_metadata"
DATA_DIR = Path("data")
OUTPUT_DIR = Path("output")
PROCESSED_DIR = Path("processed")
PLATE_DIR = PROCESSED_DIR / "cropped_plates"
SPLIT_DIR = PROCESSED_DIR / "split_wells"
BINARY_DIR = PROCESSED_DIR / "binary_wells"
MARKED_DIR = PROCESSED_DIR / "marked_aruco"

for directory in [OUTPUT_DIR, PLATE_DIR, SPLIT_DIR, BINARY_DIR, MARKED_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# experiment start date (folder name) -> batch label
DATE_CONFIG = {
    "06102025": "batch1L",
    "13102025": "batch2L",
    "20102025": "batch3L",
    "27102025": "batch4L",
    "03112025": "batch5L",
}

PIXEL_TO_MM2 = 0.0016
LIGHT_LEVELS = {"PFD100", "PFD1000"}
THIN_LINE_MIN_AREA = 200
THIN_LINE_MAX_AREA = 50000


def imread_unicode(path: Path):
    data = np.fromfile(str(path), dtype=np.uint8)
    return cv2.imdecode(data, cv2.IMREAD_COLOR)


def imwrite_unicode(path: Path, image) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix if path.suffix else ".jpg"
    ok, encoded = cv2.imencode(suffix, image)
    if not ok:
        raise RuntimeError(f"Could not encode image: {path}")
    encoded.tofile(str(path))


def is_raw_photo(path: Path) -> bool:
    parts = [part.lower() for part in path.parts]
    name = path.name.lower()
    excluded_parts = {
        "kiridashi",
        "split_images",
        "binary_images",
        "marked",
        "marked_aruco",
        "edited binary image",
        "edited binary images",
        "original binary image",
        "original binary images",
        "processed",
    }
    if any(part in excluded_parts for part in parts):
        return False
    if name.startswith(("cropped_", "binary_", "marked_")):
        return False
    return path.suffix.lower() in {".jpg", ".jpeg"}


def build_photo_index(date_dir: Path) -> dict[str, list[Path]]:
    index: dict[str, list[Path]] = {}
    for path in date_dir.rglob("*"):
        if path.is_file() and is_raw_photo(path):
            index.setdefault(path.name.lower(), []).append(path)
    return index


def choose_source_path(paths: list[Path]) -> Path:
    def score(path: Path) -> tuple[int, int, str]:
        parts = [part.lower() for part in path.parts]
        day_score = 0 if any(re.fullmatch(r"day\d+", part) for part in parts) else 1
        original_score = 0 if any("original" in part for part in parts) else 1
        return (day_score, original_score, str(path))

    return sorted(paths, key=score)[0]


def detect_day(path: Path) -> str | None:
    for part in path.parts:
        match = re.fullmatch(r"Day(\d+)", part, flags=re.IGNORECASE)
        if match:
            return f"Day{match.group(1)}"
    return None


def normalized_photo_id(filename: str) -> str:
    return Path(filename).with_suffix(".jpg").name


def crop_plate_with_aruco(image, marker_dict_name=cv2.aruco.DICT_4X4_50):
    aruco_dict = cv2.aruco.getPredefinedDictionary(marker_dict_name)
    aruco_params = cv2.aruco.DetectorParameters()
    detector = cv2.aruco.ArucoDetector(aruco_dict, aruco_params)

    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    corners, ids, _ = detector.detectMarkers(gray)
    if ids is None or len(corners) < 4:
        return None, None, "Insufficient ArUco markers"

    ids = ids.flatten()
    marker_dict = {int(marker_id): corner for marker_id, corner in zip(ids, corners)}
    required = [0, 1, 2, 3]
    if not all(marker_id in marker_dict for marker_id in required):
        return None, None, "Required ArUco markers 0-3 not detected"

    selected_points = np.array(
        [
            marker_dict[0][0][2],
            marker_dict[1][0][3],
            marker_dict[2][0][0],
            marker_dict[3][0][1],
        ],
        dtype=np.float32,
    )
    width, height = 3000, 2000
    target_points = np.array([[0, 0], [width, 0], [width, height], [0, height]], dtype=np.float32)
    transform = cv2.getPerspectiveTransform(selected_points, target_points)
    cropped = cv2.warpPerspective(image, transform, (width, height))

    marked = image.copy()
    cv2.aruco.drawDetectedMarkers(marked, corners, ids)
    return cropped, marked, None


def split_plate(plate_image) -> dict[str, np.ndarray]:
    height, width = plate_image.shape[:2]
    h_step = height // 2
    w_step = width // 3
    wells = {}
    for y in range(2):
        for x in range(3):
            well = f"y{y}-x{x}"
            wells[well] = plate_image[y * h_step : (y + 1) * h_step, x * w_step : (x + 1) * w_step]
    return wells


def make_binary_mask(well_image):
    # primary green-frond threshold: L 0-255, a 90-135, b 155-255 (8-bit CIELAB)
    lab = cv2.cvtColor(well_image, cv2.COLOR_BGR2Lab)
    l_channel, a_channel, b_channel = cv2.split(lab)

    return cv2.bitwise_and(
        cv2.bitwise_and(cv2.inRange(l_channel, 0, 255), cv2.inRange(a_channel, 90, 135)),
        cv2.inRange(b_channel, 155, 255),
    )


def remove_thin_crescent_components(mask: np.ndarray):
    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(mask, 8)
    edited = mask.copy()
    removed_pixels = 0
    removed_components = 0

    for label in range(1, num_labels):
        area = int(stats[label, cv2.CC_STAT_AREA])
        if area < THIN_LINE_MIN_AREA or area > THIN_LINE_MAX_AREA:
            continue

        width = int(stats[label, cv2.CC_STAT_WIDTH])
        height = int(stats[label, cv2.CC_STAT_HEIGHT])
        extent = area / max(1, width * height)
        aspect_ratio = max(width, height) / max(1, min(width, height))

        component_mask = np.zeros_like(mask)
        component_mask[labels == label] = 255
        contours, _ = cv2.findContours(component_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        perimeter = sum(cv2.arcLength(contour, True) for contour in contours)
        circularity = (4 * np.pi * area / (perimeter * perimeter)) if perimeter > 0 else 0

        thin_line_like = aspect_ratio >= 10 and extent <= 0.35
        crescent_like = perimeter >= 200 and circularity <= 0.12 and extent <= 0.22

        if thin_line_like or crescent_like:
            edited[labels == label] = 0
            removed_pixels += area
            removed_components += 1

    return edited, removed_pixels, removed_components


def clean_binary_mask(mask: np.ndarray):
    cleaned, thin_removed_pixels, thin_removed_components = remove_thin_crescent_components(mask)
    return {
        "mask": cleaned,
        "hole_filled_pixels": 0,
        "small_removed_pixels": 0,
        "small_removed_components": 0,
        "thin_removed_pixels": thin_removed_pixels,
        "thin_removed_components": thin_removed_components,
    }


def process_date(date_name: str, batch: str) -> tuple[list[dict], list[dict]]:
    date_dir = PHOTO_ROOT / date_name
    photoinfo = pd.read_csv(METADATA_DIR / f"photoinfo_{date_name}.csv")
    photoinfo["PhotoID"] = photoinfo["FileName"].map(normalized_photo_id)
    photoinfo["photo_time"] = pd.to_datetime(photoinfo["DateTimeOriginal"], format="%Y:%m:%d %H:%M:%S", errors="coerce")
    photoinfo = photoinfo.dropna(subset=["PhotoID", "photo_time", "plate"])
    photoinfo = photoinfo.sort_values("photo_time")
    start_time = photoinfo["photo_time"].min()
    photoinfo["start_time"] = start_time
    photoinfo["elapsed_time"] = (photoinfo["photo_time"] - start_time).dt.total_seconds() / 3600

    sample_list = pd.read_csv(METADATA_DIR / f"sample_list_{date_name}.csv")
    sample_list["sample"] = sample_list["sample"].replace({"M1OE": "M10E"})
    sample_list = sample_list[sample_list["light"].isin(LIGHT_LEVELS)].copy()

    sample_lookup = {
        (row["plate"], row["well"]): row
        for _, row in sample_list.iterrows()
    }

    photo_index = build_photo_index(date_dir)
    measurement_rows: list[dict] = []
    failure_rows: list[dict] = []

    for _, info in photoinfo.iterrows():
        filename_key = str(info["FileName"]).lower()
        if filename_key not in photo_index:
            continue

        source_path = choose_source_path(photo_index[filename_key])
        day = detect_day(source_path)
        image = imread_unicode(source_path)
        if image is None:
            failure_rows.append(
                {"date": date_name, "batch": batch, "PhotoID": info["PhotoID"],
                 "source_path": str(source_path), "reason": "Could not read image"}
            )
            continue

        plate_image, marked_image, error = crop_plate_with_aruco(image)
        if error is not None:
            failure_rows.append(
                {"date": date_name, "batch": batch, "PhotoID": info["PhotoID"],
                 "source_path": str(source_path), "reason": error}
            )
            continue

        photo_stem = Path(info["PhotoID"]).stem
        imwrite_unicode(PLATE_DIR / date_name / f"cropped_{photo_stem}.jpg", plate_image)
        imwrite_unicode(MARKED_DIR / date_name / f"marked_{photo_stem}.jpg", marked_image)

        for well, well_image in split_plate(plate_image).items():
            sample_info = sample_lookup.get((info["plate"], well))
            if sample_info is None:
                continue

            split_name = f"cropped_{photo_stem}_{well}.jpg"
            binary_name = f"binary_{split_name}"
            split_output = SPLIT_DIR / date_name / split_name
            binary_output = BINARY_DIR / date_name / binary_name
            cleanup = clean_binary_mask(make_binary_mask(well_image))
            mask = cleanup["mask"]

            imwrite_unicode(split_output, well_image)
            imwrite_unicode(binary_output, mask)

            pixel_count = int(cv2.countNonZero(mask))
            well_id = f"{info['plate']}_{well}"
            measurement_rows.append(
                {
                    "date": date_name,
                    "batch": batch,
                    "day": day,
                    "PhotoID": info["PhotoID"],
                    "plate": info["plate"],
                    "well": well,
                    "wellID": well_id,
                    "sample": sample_info["sample"],
                    "light": sample_info["light"],
                    "temp": sample_info["temp"],
                    "sampleID": f"{sample_info['sample']}_{sample_info['light']}_{sample_info['temp']}_{well_id}",
                    "start_time": start_time,
                    "photo_time": info["photo_time"],
                    "elapsed_time": info["elapsed_time"],
                    "Processed_File": split_name,
                    "Binary_File": binary_name,
                    "Frond_Pixel_Count": pixel_count,
                    "frond_area": pixel_count * PIXEL_TO_MM2,
                    "Hole_Filled_Pixel_Count": cleanup["hole_filled_pixels"],
                    "Small_Component_Removed_Pixel_Count": cleanup["small_removed_pixels"],
                    "Small_Component_Removed_Count": cleanup["small_removed_components"],
                    "Thin_Crescent_Removed_Pixel_Count": cleanup["thin_removed_pixels"],
                    "Thin_Crescent_Removed_Count": cleanup["thin_removed_components"],
                    "Mask_Cleanup_Profile": "lab_threshold_L0-255_a90-135_b155-255+remove_thin_crescent_components",
                    "source_path": str(source_path),
                    "split_image_path": str(split_output),
                    "binary_image_path": str(binary_output),
                }
            )

    return measurement_rows, failure_rows


def main() -> None:
    all_measurements: list[dict] = []
    all_failures: list[dict] = []
    for date_name, batch in DATE_CONFIG.items():
        measurements, failures = process_date(date_name, batch)
        all_measurements.extend(measurements)
        all_failures.extend(failures)
        print(f"{date_name}: measured {len(measurements)} well-images; failures {len(failures)}")

    measurements = pd.DataFrame(all_measurements).sort_values(["date", "photo_time", "plate", "well"])
    measurements.to_csv(DATA_DIR / "frond_area_main_mask_by_well_image.csv", index=False)
    pd.DataFrame(all_failures).to_csv(OUTPUT_DIR / "image_processing_failures.csv", index=False)
    print(f"Wrote {len(measurements)} well-image measurements; {len(all_failures)} failures")


if __name__ == "__main__":
    main()
