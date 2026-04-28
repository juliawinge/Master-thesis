# Create and run the hotspot analysis using the Ubuntu terminal
## Create the script
Go to a dictionary dedicated to your analysis:
```
cd /your_directory/hot_spot_analysis/
```
Create a file for your .py script:
```
nano run_hotspot.py
```
Paste the following code:
```
import csv
import re
from pathlib import Path
from collections import defaultdict

# ======= SETTINGS =======
# Directory where the script will be run
ROOT = Path("/your_dictionary/hot_spot_analysis")

# Hotspot threshold: same site must appear in  ≥ N unique patients
MIN_UNIQUE_PATIENTS = 3

# Column names expected in the TSV/TXT files
COL_SAMPLE_ID = "sample_id" 
COL_SAMPLE_NAME = "sample_name"  
COL_MUT = "cds_mutation"

# File type to scan under ROOT
ALLOWED_SUFFIXES = {".tsv", ".txt"}

# Output summary file written inside ROOT
OUTPUT_FILENAME = "hotspot_results.txt"
# ==================================================================================================

# Extract the numerical cDNA-position from the HGVS c.-notation to group mutations
# that occur at the same site. Potential variants in the 3'UTR (c.*) are excluded from this analysis. 
def pos_from_hgvs(h: str):

# Ignore empty values and 3'UTR variants (c.*). 
    h = (h or "").strip()
    if not h or h.startswith("c.*"):
        return None

# If there are indels, the positions are represented by the starting position of the region.
    m = re.match(r"^c\.(\d+)([+-])(\d+)_", h)
    if m:
        base, sign, off = int(m.group(1)), m.group(2), int(m.group(3))
        return base + off if sign == "+" else base - off


# If there are point mutations in the intronic region, the position is calculated as the exon position
# +/- intronic offset. 
    m = re.match(r"^c\.(\d+)([+-])(\d+)", h)
    if m:
        base, sign, off = int(m.group(1)), m.group(2), int(m.group(3))
        return base + off if sign == "+" else base - off

# If there are point mutations in the exonic region, the numeric cDNA position directly represents
# the mutation site. 
    m = re.match(r"^c\.(\d+)", h)
    return int(m.group(1)) if m else None
# ==================================================================================================

# Read TSV/TXT file and extract the minimal information needed for hotspot calling:
# patient (sample_name), sample (sample_id), and a numeric cDNA position. 
def read_table(path: Path):
    rows = []

# Open input file (tab separated).   
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

# Check if the file has a header.         
        if not reader.fieldnames:
            raise ValueError("Empty file or missing header.")

# Requires sample_name and cds_mutation. 
        if COL_SAMPLE_NAME not in reader.fieldnames or COL_MUT not in reader.fieldnames:
            raise ValueError(
                f"Missing required column(s). Found: {reader.fieldnames}. "
                f"Required: {COL_SAMPLE_NAME}, {COL_MUT}  (sample_id optional)."
            )
            
#  sample_id is optional but is kept for traceability. 
        has_id = COL_SAMPLE_ID in reader.fieldnames

# Go through rows and extract relevant fields. 
        for r in reader:
            sname = (r.get(COL_SAMPLE_NAME) or "").strip()
            mut = (r.get(COL_MUT) or "").strip()
            sid = (r.get(COL_SAMPLE_ID) or "").strip() if has_id else ""

# Skip rows with missing sample_name or mutation. sample_name defines unique patient -> must exist
            if not sname or not mut:
                continue

# Convert HGVS -> numerical cDNA position. 
            pos = pos_from_hgvs(mut)
            if pos is None:
                continue
                
# Store mutations linked to patient and position. 
            rows.append((pos, sname, sid, mut))

    return rows
# ==================================================================================================

# Identify hotspots by grouping by site and counting unique patients (unique sample_name) per site.
# One patient contributes at most once per site. 
def find_hotspots(rows):

    # pos -> sample_name -> {"sample_ids": set(), "muts": set ()}
    per_pos = defaultdict(lambda: defaultdict(lambda: {"sample_ids": set(), "muts": set()}))

    for pos, sname, sid, mut in rows:
        per_pos[pos][sname]["muts"].add(mut)
        if sid:
            per_pos[pos][sname]["sample_ids"].add(sid)

    hotspots = []
    for pos, patients in per_pos.items():
        n_unique = len(patients)
        if n_unique >= MIN_UNIQUE_PATIENTS:
            hotspots.append({
                "pos": pos,
                "n_unique": n_unique,
                "patients": patients,  # sample_name -> {"sample_ids", "muts"}
            })

# Sort: most patients first, then lowest position
    hotspots.sort(key=lambda h: (-h["n_unique"], h["pos"]))
    return hotspots
# ==================================================================================================

# Find all candidate input files under ROOT
def main():
    files = []
    for suf in ALLOWED_SUFFIXES:
        files.extend(ROOT.rglob(f"*{suf}"))
 
 # Avoid re-reading our own output file if it matches the suffix. 
    files = sorted([p for p in files if p.name != OUTPUT_FILENAME])

    if not files:
        print(f"No input files found {ALLOWED_SUFFIXES} under {ROOT}", flush=True)
        return

    out_path = ROOT / OUTPUT_FILENAME
    summary_lines = []

# Print run settings
    print(f" Hotspot search in {ROOT}", flush=True)
    print(f"   - patient identifier: {COL_SAMPLE_NAME}", flush=True)
    print(f"   - threshold: ≥{MIN_UNIQUE_PATIENTS} unique patients per site\n", flush=True)

    for path in files:
        rel = path.relative_to(ROOT)
        header = f"\n=== {rel} ==="
        print(header, flush=True)
        summary_lines.append(header)

# Read and parse file.
        try:
            rows = read_table(path)
        except Exception as e:
            msg = f"[ERROR] {rel}: {e}"
            print(msg, flush=True)
            summary_lines.append(msg)
            continue

        hotspots = find_hotspots(rows)

        if not hotspots:
            msg = f"No hotspots found (≥{MIN_UNIQUE_PATIENTS} unique patients at same site)."
            print(msg, flush=True)
            summary_lines.append(msg)
            continue

        msg = f"Found {len(hotspots)} hotspot(s):"
        print(msg, flush=True)
        summary_lines.append(msg)

# Print each hotspot and the patients contributing to it
        for i, h in enumerate(hotspots, 1):
            line = f"[{i}] SITE {h['pos']} (unique patients: {h['n_unique']})"
            print(line, flush=True)
            summary_lines.append(line)

# List patients, with their mutations and sample_ids. 
            for sname in sorted(h["patients"].keys()):
                muts = "; ".join(sorted(h["patients"][sname]["muts"]))
                sids = sorted(h["patients"][sname]["sample_ids"])

                if sids:
                    line2 = f"   - {sname} [sample_id: {', '.join(sids)}]: {muts}"
                else:
                    line2 = f"   - {sname}: {muts}"

                print(line2, flush=True)
                summary_lines.append(line2)

# Save the same summary as the terminal output
    out_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print(f"\n Saved summary: {out_path}", flush=True)
                

if __name__ == "__main__":
    main()

```
Save by pressing Ctrl+O, pressing Enter, followed by Ctrl+X.

## Run the Hotspot analysis
Run the program using the following command:
```
python3 -u run_hotspot.py
```
Your output will be saved as a list in the designated folder. 
Note and mark the given information regarding the target site. 
