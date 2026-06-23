## 1. Channge to the correct directory
Go to the directory containing the patient-matched rMATS result:
```
cd /your_directory/sample_download/your_dataset/rMATS/results/paired_control
```
 <br>

## 2. Create the summary script
Create a shell script that summarizes the rMATS results:
```
nano create_goi_rmats_summary.sh
```
 <br>
 Paste thr following code:
```
#!/bin/bash

# File contaning the genes of intrest
GOI="goi.txt"

# rMATS comparison folders included in the analysis
COMPARISONS=(
  "P7_Ovary"
  "P13_Ovary"
  "P16_Ovary"
  "P17_Ovary"
  "P13_Ect"
  "P16_Ect"
  "P17_Ect"
)

# rMATS even types
EVENTS=("SE" "A3SS" "A5SS" "MXE" "RI")

# Output file
echo -e "Gene\tComparison\tEvent_type\tRaw_events\tFiltered_events" > goi_rmats_summary.tsv

# Count raw and filtered events for each GOI, comparisionand event type. 
while read gene; do
  for comp in "${COMPARISONS[@]}"; do
    for ev in "${EVENTS[@]}"; do

      raw_file="${comp}/${ev}.MATS.JCEC.txt"
      filt_file="${comp}/${ev}.MATS.JCEC.filtered.txt"

      raw_count=0
      filt_count=0

      if [ -f "$raw_file" ]; then
        raw_count=$(awk -v g="$gene" 'BEGIN{FS="\t"} NR>1 {gsub(/"/, "", $3); if ($3==g) count++} END{print count+0}' "$raw_file")
      fi

      if [ -f "$filt_file" ]; then
        filt_count=$(awk -v g="$gene" 'BEGIN{FS="\t"} NR>1  {gsub(/"/, "", $3); if ($3==g) count++} END{print count+0}' "$filt_file")
      fi

      echo -e "${gene}\t${comp}\t${ev}\t${raw_count}\t${filt_count}" >> goi_rmats_summary.tsv

    done
  done
done < "$GOI"

echo "Done. Output written to goi_rmats_summary.tsv"
```

Save 
