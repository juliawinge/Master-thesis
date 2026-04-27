# Detta är en rubrik

```
kod
```

text blbla bla
**bold**
*kursiv*

# Hotspot analysis
## Create sample sheet
### Download data in TSV format
#### Based on Tissue type
Download data of interest from COSMIC (https://cancer.sanger.ac.uk/cosmic) by searching on your GOI. Choose the dataset with the highest number of samples. Scroll down to the table **"Tissue distribution"** and click on the bar under point mutations for the tissue of interest. 
<img width="981" height="293" alt="image" src="https://github.com/user-attachments/assets/41ec270f-91aa-44d6-acb1-8222136905e2" />
 
Go to **"Positive data"**, scroll to the right in the table and press on the column header so that the arrow points upwards, to sort the samples based on coordinates. 
<img width="263" height="215" alt="image" src="https://github.com/user-attachments/assets/86dadeb6-3bb4-4154-b965-8a9dc0127d94" /> 

Click on Export: **TSV**. 
<img width="995" height="280" alt="image" src="https://github.com/user-attachments/assets/382b461a-4b58-406e-80fc-5cee5be825bf" /> 

#### Based on Histology
If you instead want to look at a specific "Histology selection" like Endometriosis, you need to start out by searching the term **"Endometriosis"** in the search bar instead of your GOI. Click on the dataset with the most samples. 

You will then get a similar box as below based on your applied settings, here it is **very important** that you press  **"Go"** to actually apply the filter, otherwise it won't apply and the information will be misleading. 
<img width="987" height="362" alt="image" src="https://github.com/user-attachments/assets/af190061-3f4b-49b8-a399-5916522d61a3" />

 You will then scroll down and press on **"Genes with mutations"**
<img width="882" height="608" alt="image" src="https://github.com/user-attachments/assets/1aac4147-d5c1-4798-bc6d-9e6dad9a7e20" />

 Search for your GOI in the search bar, if there are any hits, they will be listed, if nothing is listed data is missing for the specific gene in the specific histology type.
<img width="1362" height="356" alt="image" src="https://github.com/user-attachments/assets/9204b931-dcf6-4d1e-a948-530b9a44cc11" />
If there are hits, download the tsv according to previous instructions. If there are no hits, there is nothing to do about it. 

### Create sample sheet from the downloaded TSV
Open a new Excel document/ new page,  click on *"Data"* in the top menu and then press *"From text/CSV"*.
<img width="936" height="305" alt="image" src="https://github.com/user-attachments/assets/01cd2a02-9de5-48cd-ac4e-addb5a8b9d7a" />

Click on menu and select **"All files"** to see the TSV files, then press on the TSV file of interest, followed by **"Import"**.
<img width="955" height="597" alt="image" src="https://github.com/user-attachments/assets/6c7ce6ec-73e5-4253-847b-0e982a857e85" />

 Use the following settings and press **"Read file"**.
 <img width="917" height="664" alt="image" src="https://github.com/user-attachments/assets/05a07bc6-f587-400e-946c-040e900710c9" />

Create a new Excel file/ page and copy over the columns **"Sample ID"**, **"Sample Name"** and **"CDS Mutations"**, change the corresponding column headers to **"sample_id"**, **"sample_name"** and **"cds_mutations"**, which will be used as input in the script later on. 
<img width="536" height="168" alt="image" src="https://github.com/user-attachments/assets/9b7b8d58-5a7b-423d-9a1f-622faf7b08ad" />

Save the file in **TSV** or **TXT** format in a designated folder **"/your_directory/hot_spot_analysis"**.
<img width="465" height="91" alt="image" src="https://github.com/user-attachments/assets/a086b365-e734-411f-982a-614bfeb7b865" />


## Create and run the hotspot analysis
## Create the script
Go to a dictonary dedicated to your analysis:
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

# Extract the numerical cDNA-position from the HGVS c.-notation in order to group mutations
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


# If there are point mutations in the intronic region, the position is calculated as exon position
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
 
 # Avoid re-reading our own output file if it matches suffix. 
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

# Read and prase file.
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
Save by pressing Ctrl+O, press enter followed by Ctrl+X.

### Run the Hotspot analysis
Run the program using the following command:
```
python3 -u run_hotspot.py
```
Your output will be saved as a list in the designated folder. 
Note and mark the given information regarding the hotspots. 
