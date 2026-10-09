#!/usr/bin/env bash
# Step 01: extract the columns used downstream from the five GBIF occurrence downloads.
#
# Input : gbif_raw/<download key>.csv   GBIF SIMPLE_CSV downloads of 19 September 2024
#         (not shipped; download them from GBIF, see README "Inputs")
#           Spirodela   0022526-240906103802322  https://doi.org/10.15468/dl.tqazwm
#           Landoltia   0022527-240906103802322  https://doi.org/10.15468/dl.cmkp53
#           Lemna       0022523-240906103802322  https://doi.org/10.15468/dl.sp7u55
#           Wolffiella  0022531-240906103802322  https://doi.org/10.15468/dl.wxhh5g
#           Wolffia     0022602-240906103802322  https://doi.org/10.15468/dl.shuzby
# Output: gbif_raw/<Genus>_extracted.txt   tab-separated, 10 columns
# Run from this directory:  bash 01_extract_columns.sh
set -euo pipefail
RAW=gbif_raw

extract() {  # $1 = GBIF download file, $2 = output file
  awk -F "\t" 'BEGIN {OFS="\t"}
    NR==1 {for (i=1; i<=NF; i++) h[$i]=i}
    {print $h["gbifID"], $h["occurrenceID"], $h["genus"], $h["species"], $h["scientificName"],
           $h["locality"], $h["decimalLatitude"], $h["decimalLongitude"], $h["elevation"],
           $h["coordinateUncertaintyInMeters"]}' "$1" > "$2"
}

extract "$RAW/0022526-240906103802322.csv" "$RAW/Spirodela_extracted.txt"
extract "$RAW/0022527-240906103802322.csv" "$RAW/Landoltia_extracted.txt"
extract "$RAW/0022523-240906103802322.csv" "$RAW/Lemna_extracted.txt"
extract "$RAW/0022531-240906103802322.csv" "$RAW/Wolffiella_extracted.txt"
extract "$RAW/0022602-240906103802322.csv" "$RAW/Wolffia_extracted.txt"
