# Step 02: candidate frond masks per well image + browser viewer for manual mask selection.
#
# For each well image of the six focal species (Sp, Lp, Lgp8L, M10E, Wh, Wa; 90 wells x 8 days),
# six candidate masks are computed on the original image and on a CLAHE-brightened copy
# (CLAHE clipLimit 2.5, 8 x 8 tiles on the L channel, then x1.08 + 10), i.e. 12 candidates:
#   Main current   primary mask of step 01 (L 0-255, a 90-135, b 155-255; thin/crescent removal)
#   L36 only       broader mask L 36-255, a 0-138, b 143-255
#   L36 + small200 + crescent                     + components < 200 px removed + standard thin/crescent removal
#   L36 + small500 + stronger reflection          + components < 500 px removed + stronger reflection removal
#   L36 + small500 + stronger reflection + root   + morphological opening (9 px ellipse) to remove root-like lines
#   L36 + redOR + small500 + stronger reflection + root
#                  as above, but L36 mask OR red/orange mask (L 0-255, a 136-255, b 134-244)
# Area = mask pixels x 0.0016 mm2.  The candidate overlays are written to an HTML viewer
# (viewer/index.html); the best-matching mask was chosen manually per well (optionally
# overridden per image) and exported with the viewer's "Export CSV" button as
# data/image_mask_selection.csv.
#
# Input : data/frond_area_main_mask_by_well_image.csv (step 01), processed/split_wells, processed/binary_wells
# Output: data/frond_area_candidate_masks_by_image.csv   720 images x 2 preprocessing x 6 masks
#         viewer/index.html, viewer/thumbs/                 manual selection viewer
# Run from this directory:  python 02_candidate_masks_viewer.py
# Author: Natsu Katayama
from __future__ import annotations

import json
import re
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from PIL import Image


TABLE_PATH = Path("data") / "frond_area_main_mask_by_well_image.csv"
CANDIDATE_TABLE_PATH = Path("data") / "frond_area_candidate_masks_by_image.csv"
OUTPUT_DIR = Path("viewer")
THUMB_DIR = OUTPUT_DIR / "thumbs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
THUMB_DIR.mkdir(parents=True, exist_ok=True)

FOCAL_SAMPLES = ["Sp", "Lp", "Lgp8L", "M10E", "Wh", "Wa"]
PIXEL_TO_MM2 = 0.0016
THUMB_WIDTH = 220
PREPROCESS_CHOICES = ["Original photo", "Brightened photo"]

SMALL_200 = 200
SMALL_500 = 500
THIN_LINE_MIN_AREA_STANDARD = 200
THIN_LINE_MAX_AREA_STANDARD = 50000
THIN_LINE_MIN_AREA_STRONG = 150
THIN_LINE_MAX_AREA_STRONG = 120000
ROOT_OPEN_KERNEL_SIZE = 9

OVERLAY_COLOUR = (0, 220, 255)
OVERLAY_ALPHA = 95

FILTERS = [
    "Main current",
    "L36 only",
    "L36 + small200 + crescent",
    "L36 + small500 + stronger reflection",
    "L36 + small500 + stronger reflection + root",
    "L36 + redOR + small500 + stronger reflection + root",
]


# ---- mask helpers --------------------------------------------------------------------------

def imread_unicode(path: str | Path):
    data = np.fromfile(str(path), dtype=np.uint8)
    return cv2.imdecode(data, cv2.IMREAD_COLOR)


def read_binary_mask(path: str | Path) -> np.ndarray:
    mask = imread_unicode(path)
    if mask is None:
        raise RuntimeError(f"Could not read binary image: {path}")
    if mask.ndim == 3:
        mask = cv2.cvtColor(mask, cv2.COLOR_BGR2GRAY)
    _, mask = cv2.threshold(mask, 127, 255, cv2.THRESH_BINARY)
    return mask


def lab_mask(image_bgr: np.ndarray, l_range, a_range, b_range) -> np.ndarray:
    lab = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2Lab)
    l_channel, a_channel, b_channel = cv2.split(lab)
    return cv2.bitwise_and(
        cv2.bitwise_and(
            cv2.inRange(l_channel, l_range[0], l_range[1]),
            cv2.inRange(a_channel, a_range[0], a_range[1]),
        ),
        cv2.inRange(b_channel, b_range[0], b_range[1]),
    )


def remove_small_components(mask: np.ndarray, min_area: int):
    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(mask, 8)
    if num_labels <= 1:
        return mask.copy(), 0, 0

    areas = stats[:, cv2.CC_STAT_AREA]
    remove_labels = np.flatnonzero((areas < min_area) & (np.arange(num_labels) != 0))
    removed_components = int(len(remove_labels))
    removed_pixels = int(areas[remove_labels].sum()) if removed_components else 0

    edited = mask.copy()
    if removed_components:
        edited[np.isin(labels, remove_labels)] = 0
    return edited, removed_pixels, removed_components


def remove_shape_components(mask: np.ndarray, min_area: int, max_area: int, strong: bool):
    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(mask, 8)
    edited = mask.copy()
    removed_pixels = 0
    removed_components = 0

    for label in range(1, num_labels):
        area = int(stats[label, cv2.CC_STAT_AREA])
        if area < min_area or area > max_area:
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

        if strong:
            remove_component = (
                (aspect_ratio >= 6 and extent <= 0.45)
                or (perimeter >= 150 and circularity <= 0.18 and extent <= 0.30)
                or (aspect_ratio >= 3.5 and extent <= 0.18 and circularity <= 0.25)
            )
        else:
            remove_component = (
                (aspect_ratio >= 10 and extent <= 0.35)
                or (perimeter >= 200 and circularity <= 0.12 and extent <= 0.22)
            )

        if remove_component:
            edited[labels == label] = 0
            removed_pixels += area
            removed_components += 1

    return edited, removed_pixels, removed_components


def remove_root_like_lines(mask: np.ndarray):
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (ROOT_OPEN_KERNEL_SIZE, ROOT_OPEN_KERNEL_SIZE))
    opened = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)
    removed = cv2.bitwise_and(mask, cv2.bitwise_not(opened))
    removed_pixels = int(cv2.countNonZero(removed))
    num_labels, _, _, _ = cv2.connectedComponentsWithStats(removed, 8)
    removed_components = max(0, int(num_labels) - 1)
    return opened, removed_pixels, removed_components


def make_overlay(original: Image.Image, mask: np.ndarray) -> Image.Image:
    binary = Image.fromarray(mask).convert("L")
    if binary.size != original.size:
        binary = binary.resize(original.size, Image.Resampling.NEAREST)
    overlay_colour = Image.new("RGB", original.size, OVERLAY_COLOUR)
    alpha = binary.point(lambda value: OVERLAY_ALPHA if value > 0 else 0)
    overlay = original.copy()
    overlay.paste(overlay_colour, (0, 0), alpha)
    return overlay


# ---- candidate masks -----------------------------------------------------------------------

def safe_id(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_")


def save_jpeg(image: Image.Image, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path, "JPEG", quality=80, optimize=True)


def resize(image: Image.Image, width: int = THUMB_WIDTH) -> Image.Image:
    if image.width <= width:
        return image.copy()
    ratio = width / image.width
    return image.resize((width, int(image.height * ratio)))


def red_orange_lab_mask(image_bgr: np.ndarray) -> np.ndarray:
    return lab_mask(image_bgr, (0, 255), (136, 255), (134, 244))


def brighten_image(image_bgr: np.ndarray) -> np.ndarray:
    lab = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2Lab)
    l_channel, a_channel, b_channel = cv2.split(lab)
    clahe = cv2.createCLAHE(clipLimit=2.5, tileGridSize=(8, 8))
    l_equalized = clahe.apply(l_channel)
    bright_lab = cv2.merge([l_equalized, a_channel, b_channel])
    bright_bgr = cv2.cvtColor(bright_lab, cv2.COLOR_Lab2BGR)
    return cv2.convertScaleAbs(bright_bgr, alpha=1.08, beta=10)


def main_like_mask(image_bgr: np.ndarray) -> np.ndarray:
    mask = lab_mask(image_bgr, (0, 255), (90, 135), (155, 255))
    cleaned, _, _ = remove_shape_components(
        mask, THIN_LINE_MIN_AREA_STANDARD, THIN_LINE_MAX_AREA_STANDARD, strong=False
    )
    return cleaned


def make_masks(image_bgr: np.ndarray, row, use_saved_main: bool) -> dict[str, np.ndarray]:
    masks: dict[str, np.ndarray] = {}
    # original photo: primary mask saved by step 01; brightened photo: recomputed
    masks["Main current"] = read_binary_mask(row["binary_image_path"]) if use_saved_main else main_like_mask(image_bgr)

    l36 = lab_mask(image_bgr, (36, 255), (0, 138), (143, 255))
    masks["L36 only"] = l36

    small200, _, _ = remove_small_components(l36, SMALL_200)
    crescent200, _, _ = remove_shape_components(
        small200, THIN_LINE_MIN_AREA_STANDARD, THIN_LINE_MAX_AREA_STANDARD, strong=False
    )
    masks["L36 + small200 + crescent"] = crescent200

    small500, _, _ = remove_small_components(l36, SMALL_500)
    strong_reflection, _, _ = remove_shape_components(
        small500, THIN_LINE_MIN_AREA_STRONG, THIN_LINE_MAX_AREA_STRONG, strong=True
    )
    masks["L36 + small500 + stronger reflection"] = strong_reflection

    root_removed, _, _ = remove_root_like_lines(strong_reflection)
    masks["L36 + small500 + stronger reflection + root"] = root_removed

    combined = cv2.bitwise_or(l36, red_orange_lab_mask(image_bgr))
    small500, _, _ = remove_small_components(combined, SMALL_500)
    strong_reflection, _, _ = remove_shape_components(
        small500, THIN_LINE_MIN_AREA_STRONG, THIN_LINE_MAX_AREA_STRONG, strong=True
    )
    root_removed, _, _ = remove_root_like_lines(strong_reflection)
    masks["L36 + redOR + small500 + stronger reflection + root"] = root_removed
    return masks


def make_filter_records(image_bgr, display_image, row, row_id, preprocess_id, use_saved_main) -> list[dict]:
    masks = make_masks(image_bgr, row, use_saved_main=use_saved_main)
    records = []
    for filter_name in FILTERS:
        mask = masks[filter_name]
        rel_path = Path("thumbs") / row_id / preprocess_id / (safe_id(filter_name) + ".jpg")
        save_jpeg(resize(make_overlay(display_image, mask)), OUTPUT_DIR / rel_path)
        records.append(
            {
                "name": filter_name,
                "thumb": str(rel_path),
                "area": round(int(cv2.countNonZero(mask)) * PIXEL_TO_MM2, 4),
            }
        )
    return records


def make_record(row, row_index: int) -> dict:
    image_bgr = imread_unicode(row["split_image_path"])
    if image_bgr is None:
        raise RuntimeError(f"Could not read split image: {row['split_image_path']}")

    original_full = Image.fromarray(cv2.cvtColor(image_bgr, cv2.COLOR_BGR2RGB))
    bright_bgr = brighten_image(image_bgr)
    bright_full = Image.fromarray(cv2.cvtColor(bright_bgr, cv2.COLOR_BGR2RGB))

    row_id = safe_id(
        f"{row_index:04d}_{row['date']}_{row['day']}_{row['sample']}_{row['light']}_{row['wellID']}_{Path(row['PhotoID']).stem}"
    )

    original_rel = Path("thumbs") / row_id / "original.jpg"
    save_jpeg(resize(original_full), OUTPUT_DIR / original_rel)
    bright_rel = Path("thumbs") / row_id / "brightened.jpg"
    save_jpeg(resize(bright_full), OUTPUT_DIR / bright_rel)

    return {
        "id": row_id,
        "wellKey": f"{row['batch']}|{row['wellID']}",
        "date": str(row["date"]).zfill(8),
        "day": str(row["day"]),
        "batch": str(row["batch"]),
        "sample": str(row["sample"]),
        "light": str(row["light"]),
        "wellID": str(row["wellID"]),
        "PhotoID": str(row["PhotoID"]),
        "elapsed_time": float(row["elapsed_time"]),
        "original": str(original_rel),
        "brightened": str(bright_rel),
        "filters": make_filter_records(image_bgr, original_full, row, row_id, "original", use_saved_main=True),
        "brightened_filters": make_filter_records(bright_bgr, bright_full, row, row_id, "brightened", use_saved_main=False),
    }


# ---- manual-selection viewer (HTML) --------------------------------------------------------

def build_html(records: list[dict]) -> str:
    data_json = json.dumps(records, ensure_ascii=False, separators=(",", ":"))
    filter_json = json.dumps(FILTERS, ensure_ascii=False)
    sample_json = json.dumps(FOCAL_SAMPLES, ensure_ascii=False)
    preprocess_json = json.dumps(PREPROCESS_CHOICES, ensure_ascii=False)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Fig3 6-species well filter selector</title>
<style>
:root {{
  color-scheme: light;
  font-family: -apple-system, BlinkMacSystemFont, "Helvetica Neue", Arial, sans-serif;
  --bg: #f7f7f4;
  --fg: #20211f;
  --muted: #666960;
  --line: #cfcfc7;
  --card: #ffffff;
  --accent: #0b6b7a;
  --accent-soft: #dff3f6;
}}
body {{
  margin: 0;
  background: var(--bg);
  color: var(--fg);
}}
header {{
  position: sticky;
  top: 0;
  z-index: 10;
  background: rgba(247, 247, 244, 0.96);
  border-bottom: 1px solid var(--line);
  padding: 12px 16px;
}}
h1 {{
  font-size: 18px;
  margin: 0 0 10px;
  font-weight: 600;
}}
.controls, .well-actions {{
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  align-items: end;
}}
label {{
  display: grid;
  gap: 3px;
  font-size: 12px;
  color: var(--muted);
}}
select, input, button, textarea {{
  font: inherit;
  border: 1px solid var(--line);
  border-radius: 6px;
  background: white;
  color: var(--fg);
  padding: 7px 9px;
}}
button {{
  cursor: pointer;
}}
button.primary {{
  background: var(--accent);
  color: white;
  border-color: var(--accent);
}}
main {{
  padding: 14px 16px 28px;
}}
.status {{
  color: var(--muted);
  font-size: 13px;
  margin: 6px 0 12px;
}}
.well-card {{
  background: var(--card);
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 12px;
  margin-bottom: 14px;
  box-shadow: 0 1px 4px rgba(0,0,0,0.06);
}}
.well-head {{
  display: flex;
  flex-wrap: wrap;
  justify-content: space-between;
  gap: 10px;
  align-items: baseline;
  margin-bottom: 10px;
}}
.well-head strong {{
  font-size: 16px;
}}
.well-head span {{
  color: var(--muted);
  font-size: 13px;
}}
.filter-btn[aria-pressed="true"] {{
  background: var(--accent);
  color: white;
  border-color: var(--accent);
}}
.note {{
  width: min(560px, 100%);
  margin: 10px 0;
}}
.image-row {{
  display: grid;
  grid-template-columns: 120px repeat(7, minmax(150px, 1fr));
  gap: 8px;
  align-items: start;
  border-top: 1px solid var(--line);
  padding-top: 8px;
  margin-top: 8px;
}}
.time-meta {{
  color: var(--muted);
  font-size: 12px;
  line-height: 1.35;
}}
.tile {{
  border: 2px solid transparent;
  border-radius: 7px;
  padding: 5px;
  background: #fafafa;
  text-align: left;
}}
.tile.selected {{
  border-color: var(--accent);
  background: var(--accent-soft);
}}
.tile.override {{
  border-style: dashed;
}}
.tile-button {{
  cursor: pointer;
  width: 100%;
}}
.tile img {{
  width: 100%;
  display: block;
  border-radius: 4px;
}}
.tile-title {{
  display: flex;
  justify-content: space-between;
  gap: 6px;
  margin-top: 4px;
  font-size: 11px;
  color: var(--muted);
}}
.tile-title b {{
  color: var(--fg);
  font-weight: 500;
}}
</style>
</head>
<body>
<header>
  <h1>Fig3 6-species well-level filter selector</h1>
  <div class="controls">
    <label>Sample<select id="sampleFilter"></select></label>
    <label>Light<select id="lightFilter"></select></label>
    <label>Search<input id="searchBox" type="search" placeholder="wellID / PhotoID"></label>
    <label>Apply filter to shown wells<select id="bulkFilter"></select></label>
    <label>Apply image to shown wells<select id="bulkPreprocess"></select></label>
    <button id="bulkApplyBtn" type="button">Apply to shown wells</button>
    <button id="exportBtn" class="primary" type="button">Export CSV</button>
  </div>
</header>
<main>
  <div id="status" class="status"></div>
  <section id="viewer"></section>
</main>
<script>
const records = {data_json};
const filterNames = {filter_json};
const samples = {sample_json};
const preprocessChoices = {preprocess_json};
const storageKey = "fig3-6species-well-filter-selection-v3";
const defaultFilter = "Main current";
const defaultPreprocess = "Original photo";
let selections = JSON.parse(localStorage.getItem(storageKey) || "{{}}");

const sampleFilter = document.getElementById("sampleFilter");
const lightFilter = document.getElementById("lightFilter");
const searchBox = document.getElementById("searchBox");
const bulkFilter = document.getElementById("bulkFilter");
const bulkPreprocess = document.getElementById("bulkPreprocess");
const statusEl = document.getElementById("status");
const viewer = document.getElementById("viewer");

function fillSelect(select, values) {{
  select.innerHTML = "";
  values.forEach(value => {{
    const opt = document.createElement("option");
    opt.value = value;
    opt.textContent = value;
    select.appendChild(opt);
  }});
}}
fillSelect(sampleFilter, samples);
fillSelect(lightFilter, ["PFD100", "PFD1000"]);
fillSelect(bulkFilter, filterNames);
fillSelect(bulkPreprocess, preprocessChoices);

function saveSelections() {{
  localStorage.setItem(storageKey, JSON.stringify(selections));
}}

function normalizeSelectionStore() {{
  if (!selections.wells && !selections.images) {{
    selections = {{wells: selections, images: {{}}}};
  }}
  selections.wells = selections.wells || {{}};
  selections.images = selections.images || {{}};
}}
normalizeSelectionStore();

function groupRecords() {{
  const sample = sampleFilter.value;
  const light = lightFilter.value;
  const q = searchBox.value.trim().toLowerCase();
  const subset = records.filter(d => {{
    if (d.sample !== sample || d.light !== light) return false;
    if (q && !(d.wellID.toLowerCase().includes(q) || d.PhotoID.toLowerCase().includes(q))) return false;
    return true;
  }});
  const groups = new Map();
  subset.forEach(d => {{
    if (!groups.has(d.wellKey)) groups.set(d.wellKey, []);
    groups.get(d.wellKey).push(d);
  }});
  return Array.from(groups.entries()).map(([wellKey, rows]) => {{
    rows.sort((a, b) => a.elapsed_time - b.elapsed_time || a.PhotoID.localeCompare(b.PhotoID));
    return {{wellKey, rows}};
  }}).sort((a, b) => a.wellKey.localeCompare(b.wellKey));
}}

function filterListFor(row, preprocess) {{
  return preprocess === "Brightened photo" ? row.brightened_filters : row.filters;
}}

function meanSelectedArea(rows, filterName, preprocess) {{
  filterName = filterName || defaultFilter;
  preprocess = preprocess || defaultPreprocess;
  const values = rows.map(row => filterListFor(row, preprocess).find(f => f.name === filterName)?.area).filter(v => Number.isFinite(v));
  if (!values.length) return "";
  return (values.reduce((a, b) => a + b, 0) / values.length).toFixed(3);
}}

function effectiveSelection(row, wellSelection) {{
  const imageSelection = selections.images[row.id] || {{}};
  return {{
    selected_filter: imageSelection.selected_filter || wellSelection.selected_filter || defaultFilter,
    selected_preprocess: imageSelection.selected_preprocess || wellSelection.selected_preprocess || defaultPreprocess,
    is_image_override: Boolean(imageSelection.selected_filter || imageSelection.selected_preprocess)
  }};
}}

function meanEffectiveArea(rows, wellSelection) {{
  const values = rows.map(row => {{
    const effective = effectiveSelection(row, wellSelection);
    return filterListFor(row, effective.selected_preprocess).find(f => f.name === effective.selected_filter)?.area;
  }}).filter(v => Number.isFinite(v));
  if (!values.length) return "";
  return (values.reduce((a, b) => a + b, 0) / values.length).toFixed(3);
}}

function setSelection(wellKey, filterName, preprocess) {{
  const prev = selections.wells[wellKey] || {{}};
  selections.wells[wellKey] = {{
    selected_filter: filterName,
    selected_preprocess: preprocess || prev.selected_preprocess || "Original photo",
    note: prev.note || "",
    updated_at: new Date().toISOString()
  }};
  saveSelections();
  render();
}}

function setNote(wellKey, note) {{
  const prev = selections.wells[wellKey] || {{}};
  selections.wells[wellKey] = {{
    selected_filter: prev.selected_filter || "",
    selected_preprocess: prev.selected_preprocess || "Original photo",
    note,
    updated_at: new Date().toISOString()
  }};
  saveSelections();
}}

function setPreprocess(wellKey, preprocess) {{
  const prev = selections.wells[wellKey] || {{}};
  selections.wells[wellKey] = {{
    selected_filter: prev.selected_filter || filterNames[0],
    selected_preprocess: preprocess,
    note: prev.note || "",
    updated_at: new Date().toISOString()
  }};
  saveSelections();
  render();
}}

function setImageSelection(imageId, filterName, preprocess) {{
  const prev = selections.images[imageId] || {{}};
  selections.images[imageId] = {{
    selected_filter: filterName,
    selected_preprocess: preprocess || prev.selected_preprocess || defaultPreprocess,
    updated_at: new Date().toISOString()
  }};
  saveSelections();
  render();
}}

function setImagePreprocess(imageId, preprocess) {{
  const prev = selections.images[imageId] || {{}};
  selections.images[imageId] = {{
    selected_filter: prev.selected_filter || defaultFilter,
    selected_preprocess: preprocess,
    updated_at: new Date().toISOString()
  }};
  saveSelections();
  render();
}}

function clearImageSelection(imageId) {{
  delete selections.images[imageId];
  saveSelections();
  render();
}}

function applyToShown() {{
  const groups = groupRecords();
  const filterName = bulkFilter.value;
  const preprocess = bulkPreprocess.value;
  groups.forEach(group => {{
    const prev = selections.wells[group.wellKey] || {{}};
    selections.wells[group.wellKey] = {{
      selected_filter: filterName,
      selected_preprocess: preprocess,
      note: prev.note || "",
      updated_at: new Date().toISOString()
    }};
  }});
  saveSelections();
  render();
}}

function renderWell(group) {{
  const rows = group.rows;
  const first = rows[0];
  const selection = selections.wells[group.wellKey] || {{}};
  const selectedFilter = selection.selected_filter || defaultFilter;
  const selectedPreprocess = selection.selected_preprocess || defaultPreprocess;
  const buttons = filterNames.map(name => `
    <button type="button" class="filter-btn" aria-pressed="${{selectedFilter === name}}" data-well="${{group.wellKey}}" data-filter="${{name}}">
      ${{name}}
    </button>
  `).join("");
  const preprocessButtons = preprocessChoices.map(name => `
    <button type="button" class="filter-btn" aria-pressed="${{selectedPreprocess === name}}" data-well-preprocess="${{group.wellKey}}" data-preprocess="${{name}}">
      ${{name}}
    </button>
  `).join("");
  const imageRows = rows.map(row => {{
    const effective = effectiveSelection(row, selection);
    const selectedName = effective.selected_filter;
    const rowPreprocess = effective.selected_preprocess;
    const filterRecords = filterListFor(row, rowPreprocess);
    const displayOriginal = rowPreprocess === "Brightened photo" ? row.brightened : row.original;
    const rowPreprocessButtons = preprocessChoices.map(name => `
      <button type="button" class="filter-btn" aria-pressed="${{rowPreprocess === name}}" data-image-preprocess="${{row.id}}" data-preprocess="${{name}}">
        ${{name.replace(" photo", "")}}
      </button>
    `).join("");
    const filterTiles = filterRecords.map(f => `
      <button type="button" class="tile tile-button ${{selectedName === f.name ? "selected" : ""}} ${{effective.is_image_override ? "override" : ""}}" data-image="${{row.id}}" data-filter="${{f.name}}" data-preprocess="${{rowPreprocess}}">
        <img src="${{f.thumb}}" alt="${{f.name}} overlay">
        <div class="tile-title"><b>${{f.name.replace("L36 + ", "+ ")}}</b><span>${{f.area.toFixed(1)}}</span></div>
      </button>
    `).join("");
    return `
      <div class="image-row">
        <div class="time-meta">
          <b>${{row.day}}</b><br>${{row.elapsed_time.toFixed(1)}} h<br>${{row.PhotoID}}<br>
          ${{effective.is_image_override ? "<b>day override</b><br>" : "inherits well<br>"}}
          <div class="well-actions">${{rowPreprocessButtons}}</div>
          <button type="button" data-clear-image="${{row.id}}">Clear day</button>
        </div>
        <div class="tile"><img src="${{displayOriginal}}" alt="${{rowPreprocess}}"><div class="tile-title"><b>${{rowPreprocess}}</b><span></span></div></div>
        ${{filterTiles}}
      </div>
    `;
  }}).join("");
  return `
    <article class="well-card">
      <div class="well-head">
        <strong>${{first.sample}} ${{first.light}} | ${{first.batch}} ${{first.wellID}}</strong>
        <span>${{rows.length}} images | well default: ${{selectedPreprocess}} / ${{selectedFilter}} | effective mean area: ${{meanEffectiveArea(rows, selection)}} mm2</span>
      </div>
      <div class="well-actions">${{preprocessButtons}}</div>
      <div class="well-actions">${{buttons}}</div>
      <label class="note">Note
        <textarea data-note="${{group.wellKey}}" placeholder="optional note">${{selection.note || ""}}</textarea>
      </label>
      ${{imageRows}}
    </article>
  `;
}}

function render() {{
  const groups = groupRecords();
  const selectedCount = Object.values(selections.wells).filter(v => v.selected_filter).length;
  const imageOverrideCount = Object.keys(selections.images).length;
  statusEl.textContent = `${{sampleFilter.value}} / ${{lightFilter.value}}: ${{groups.length}} wells shown | ${{selectedCount}} wells selected | ${{imageOverrideCount}} day overrides`;
  viewer.innerHTML = groups.map(renderWell).join("") || "<p>No wells match the current filters.</p>";
  viewer.querySelectorAll("[data-filter]").forEach(btn => {{
    btn.addEventListener("click", () => {{
      const prev = selections.wells[btn.dataset.well] || {{}};
      setSelection(btn.dataset.well, btn.dataset.filter, prev.selected_preprocess || "Original photo");
    }});
  }});
  viewer.querySelectorAll("[data-well-preprocess]").forEach(btn => {{
    btn.addEventListener("click", () => setPreprocess(btn.dataset.wellPreprocess, btn.dataset.preprocess));
  }});
  viewer.querySelectorAll("[data-note]").forEach(textarea => {{
    textarea.addEventListener("input", () => setNote(textarea.dataset.note, textarea.value));
  }});
  viewer.querySelectorAll("[data-image]").forEach(btn => {{
    btn.addEventListener("click", () => setImageSelection(btn.dataset.image, btn.dataset.filter, btn.dataset.preprocess));
  }});
  viewer.querySelectorAll("[data-image-preprocess]").forEach(btn => {{
    btn.addEventListener("click", () => setImagePreprocess(btn.dataset.imagePreprocess, btn.dataset.preprocess));
  }});
  viewer.querySelectorAll("[data-clear-image]").forEach(btn => {{
    btn.addEventListener("click", () => clearImageSelection(btn.dataset.clearImage));
  }});
}}

function exportCsv() {{
  const groupsAll = new Map();
  records.forEach(row => {{
    if (!groupsAll.has(row.wellKey)) groupsAll.set(row.wellKey, []);
    groupsAll.get(row.wellKey).push(row);
  }});
  const rows = Array.from(groupsAll.entries()).map(([wellKey, imageRows]) => {{
    imageRows.sort((a, b) => a.elapsed_time - b.elapsed_time);
    const first = imageRows[0];
    const selection = selections.wells[wellKey] || {{}};
    const selectedFilter = selection.selected_filter || defaultFilter;
    const preprocess = selection.selected_preprocess || defaultPreprocess;
    return [
      wellKey,
      first.batch,
      first.sample,
      first.light,
      first.wellID,
      imageRows.length,
      selectedFilter,
      preprocess,
      meanSelectedArea(imageRows, selectedFilter, preprocess),
      (selection.note || "").replaceAll("\\n", " "),
      selection.updated_at || ""
    ];
  }});
  const header = ["wellKey","batch","sample","light","wellID","n_images","selected_filter","selected_preprocess","selected_mean_area","note","updated_at"];
  const csv = [header, ...rows]
    .map(row => row.map(value => `"${{String(value).replaceAll('"', '""')}}"`).join(","))
    .join("\\n");
  const blob = new Blob([csv], {{type: "text/csv;charset=utf-8"}});
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "fig3_6species_well_filter_selection.csv";
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);

  const imageRows = records.map(row => {{
    const wellSelection = selections.wells[row.wellKey] || {{}};
    const effective = effectiveSelection(row, wellSelection);
    const selectedArea = filterListFor(row, effective.selected_preprocess).find(f => f.name === effective.selected_filter)?.area ?? "";
    const imageSelection = selections.images[row.id] || {{}};
    return [
      row.id,
      row.wellKey,
      row.date,
      row.day,
      row.batch,
      row.sample,
      row.light,
      row.wellID,
      row.PhotoID,
      row.elapsed_time,
      effective.selected_filter,
      effective.selected_preprocess,
      selectedArea,
      effective.is_image_override ? "image" : "well",
      imageSelection.updated_at || wellSelection.updated_at || ""
    ];
  }});
  const imageHeader = ["id","wellKey","date","day","batch","sample","light","wellID","PhotoID","elapsed_time","selected_filter","selected_preprocess","selected_area","selection_scope","updated_at"];
  const imageCsv = [imageHeader, ...imageRows]
    .map(row => row.map(value => `"${{String(value).replaceAll('"', '""')}}"`).join(","))
    .join("\\n");
  const imageBlob = new Blob([imageCsv], {{type: "text/csv;charset=utf-8"}});
  const imageUrl = URL.createObjectURL(imageBlob);
  const imageA = document.createElement("a");
  imageA.href = imageUrl;
  imageA.download = "fig3_6species_image_filter_selection.csv";
  document.body.appendChild(imageA);
  imageA.click();
  imageA.remove();
  URL.revokeObjectURL(imageUrl);
}}

sampleFilter.addEventListener("change", render);
lightFilter.addEventListener("change", render);
searchBox.addEventListener("input", render);
document.getElementById("exportBtn").addEventListener("click", exportCsv);
document.getElementById("bulkApplyBtn").addEventListener("click", applyToShown);
render();
</script>
</body>
</html>
"""


def main() -> None:
    df = pd.read_csv(TABLE_PATH, dtype={"date": str, "PhotoID": str, "day": str})
    df["date"] = df["date"].str.zfill(8)
    df["sample"] = df["sample"].replace({"M1OE": "M10E"})
    df = df[df["sample"].isin(FOCAL_SAMPLES)].copy()
    df = df.sort_values(["sample", "light", "batch", "wellID", "elapsed_time", "PhotoID"]).reset_index(drop=True)

    records = []
    for idx, row in df.iterrows():
        records.append(make_record(row, idx))
        if (idx + 1) % 100 == 0:
            print(f"Prepared {idx + 1} / {len(df)} images")

    with (OUTPUT_DIR / "index.html").open("w", encoding="utf-8") as handle:
        handle.write(build_html(records))

    candidate_rows = []
    for record in records:
        for preprocess, key in [("Original photo", "filters"), ("Brightened photo", "brightened_filters")]:
            for candidate in record[key]:
                candidate_rows.append(
                    {
                        **{k: record[k] for k in ["id", "date", "day", "batch", "sample", "light", "wellID", "PhotoID", "elapsed_time"]},
                        "preprocess": preprocess,
                        "filter": candidate["name"],
                        "area_mm2": candidate["area"],
                    }
                )
    pd.DataFrame(candidate_rows).to_csv(CANDIDATE_TABLE_PATH, index=False)

    print(f"Wrote viewer: {OUTPUT_DIR / 'index.html'}")
    print(f"Wrote {len(candidate_rows)} candidate areas: {CANDIDATE_TABLE_PATH}")


if __name__ == "__main__":
    main()
