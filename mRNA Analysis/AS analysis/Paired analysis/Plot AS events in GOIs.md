# Plot rMATS events in GOIs from matched samples using UBUNTU
*This step summarizes significant alternative splicing (AS) events in selected genes of interest (GOIs) from patient-matched ectopic lesion samples, and visualizes the results as stacked bar plots.*

*Only patients with both ovarian and peritoneal lesion samples, as well as a corresponding eutopic endometrium sample, were included. For each patient, rMATS output folders should already have been generated from paired comparisons between the ovarian lesion and the eutopic endometrium, and between the peritoneal lesion and the eutopic endometrium.*

The workflow contains two main steps:
1. Create a summary table called `goi_rmats_summary.tsv`.
2. Use the summary table to generate a stacked bar plot.

## Requirements


## 1. Go to the analysis  directory
Go to the directory containing the patient-matched rMATS result:
```
cd /your_directory/sample_download/your_datasetrMATS/results/paired_control
```
 <br>

The directory should contain one folder for each paired comparison. For example:
```
P7_Ovary
P7_Ect
P13_Ovary
P13_Ect
P16_Ovary
P16_Ect
```
Here, `P13_Ovary` represents the rMATS comparison between patient 13 ovarian lesion and its 13 eutopic endometrium. Similarly, `P13_Ect` represents the rMATS comparison between patient 13 peritonial lesion and its eutopic endometrium.

Each comparison folder should contain rMATS output files for all event types, for example: 
```
SE.MATS.JCEC.txt 
SE.MATS.JCEC.filtered.txt 
A3SS.MATS.JCEC.txt 
A3SS.MATS.JCEC.filtered.txt 
A5SS.MATS.JCEC.txt 
A5SS.MATS.JCEC.filtered.txt 
MXE.MATS.JCEC.txt 
MXE.MATS.JCEC.filtered.txt 
RI.MATS.JCEC.txt 
RI.MATS.JCEC.filtered.txt
```
## 2. Create a GOI list
Create a text file contaning the genes of intres:

Save the file by pressing `Ctrl+O`, then `Enter`, followed by `Ctrl+X`.

Only patients with both ovarian and peritoneal lesion samples, as well as a corresponding eutopic endometrium sample, were included.* <br>

Allt detta är skapat i en egen mapp för endast paired samples
/your_directory/sample_download/your_datasetrMATS/results/paired_control

Okej för de patienter som har peritoniala och ovariala lesions prover, samt tillhörande interna kontroll dvs eutopic endometium. 

Utgår från goi_rmats_summary.tsv vilkety är en sammanfattning av (paired) filtered as events för  between ovarial och eutopic samt  peritoneal och eutopic i varje GOI.

Skapad via: 
```
#!/bin/bash

GOI="goi.txt"

COMPARISONS=(
  "P7_Ovary"
  "P13_Ovary"
  "P16_Ovary"
  "P17_Ovary"
  "P13_Ect"
  "P16_Ect"
  "P17_Ect"
)

EVENTS=("SE" "A3SS" "A5SS" "MXE" "RI")

echo -e "Gene\tComparison\tEvent_type\tRaw_events\tFiltered_events" > goi_rmats_summary.tsv

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
där "P7_Ovary" är output map från rMATS results för P7_Ovary vs P7_Eutopic, osv för varje patient och provtyp som listas. 

Den använder goi.txt som input, den ska ligga i den desegnedade mappen och är i följande format SLC1A1
STAT3
TCHH
CYP4F11
ENOSF1, dvs endast dina listade GOIs

Det skapas då en goi_rmats_summary.tsv i /your_directory/sample_download/your_datasetrMATS/results/paired_control

Gene	Comparison	Event_type	Raw_events	Filtered_events
SLC1A1	P7_Ovary	SE	2	0
SLC1A1	P7_Ovary	A3SS	0	0
SLC1A1	P7_Ovary	A5SS	0	0
SLC1A1	P7_Ovary	MXE	0	0
SLC1A1	P7_Ovary	RI	0	0
SLC1A1	P13_Ovary	SE	1	0






*The script ....*

## 1. Go to the analysis directory
Go to a directory dedicated to your analysis:
```
cd /your_directory/sample_download/your_dataset/rMATS.....
```
 <br>
 
The directory should contain the simplified input files created from downloaded COSMIC data. 

## 2. Create a Python script
Create a new Python script:
```
nano check_goi_per_person.py
```
 <br>
 
Paste the following code:
```
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# ===== FONT =====
font_path = "/usr/share/fonts/truetype/msttcorefonts/Georgia.ttf"
fm.fontManager.addfont(font_path)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Georgia"],
    "font.size": 18,
    "axes.titlesize": 22,
    "axes.labelsize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 18,
    "legend.title_fontsize": 18,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none"
})

# ===== DATA =====
df = pd.read_csv("goi_rmats_summary.tsv", sep="\t")

events = ["SE", "A3SS", "A5SS", "MXE", "RI"]

genes = ["CYP4F11", "ENOSF1", "SLC1A1", "STAT3", "TCHH"]

# ===== COLORS =====
event_colors = {
    "SE": "#DD8452",
    "A3SS": "#4C72B0",
    "A5SS": "#55A868",
    "MXE": "#B39FDD",
    "RI": "#D62728"
}

# ===== PATIENTS AND COMPARISONS =====
patients = ["13", "16", "17"]

comparisons = {
    "13": {"Ovary": ["P13_Ovary"], "Peritoneal": ["P13_Ect"]},
    "16": {"Ovary": ["P16_Ovary"], "Peritoneal": ["P16_Ect"]},
    "17": {"Ovary": ["P17_Ovary"], "Peritoneal": ["P17_Ect"]}
}

# ===== FUNCTION =====
def plot_panel(ax, data, comparisons, value_col, title, ylim, ylabel=False):
    sub = data[data["Comparison"].isin(comparisons)]

    pivot = (
        sub.groupby(["Gene", "Event_type"])[value_col]
        .sum()
        .unstack(fill_value=0)
    )

    # Ensure all genes are present, even if they have zero events
    pivot = pivot.reindex(genes, fill_value=0)

    # Ensure all event types are present
    for ev in events:
        if ev not in pivot.columns:
            pivot[ev] = 0

    pivot = pivot[events]

    pivot.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        edgecolor="black",
        color=[event_colors[ev] for ev in events],
        width=0.7
    )

    ax.set_title(title, fontsize=18, loc="left")
    ax.set_xlabel("")
    ax.set_ylim(0, ylim)

    if ylabel:
        ax.set_ylabel("Number of events", fontsize=18)
    else:
        ax.set_ylabel("")

    ax.tick_params(axis="x", rotation=45, labelsize=15)
    ax.tick_params(axis="y", labelsize=15)

# ===== FIGURE =====
fig, axes = plt.subplots(
    nrows=3,
    ncols=2,
    figsize=(13, 13),
    sharey=True
)

ylim = 5

for row, patient in enumerate(patients):
    plot_panel(
        axes[row, 0],
        df,
        comparisons[patient]["Ovary"],
        "Filtered_events",
        f"{chr(65 + row*2)}) Patient {patient} - Ovary",
        ylim,
        ylabel=True
    )

    plot_panel(
        axes[row, 1],
        df,
        comparisons[patient]["Peritoneal"],
        "Filtered_events",
        f"{chr(66 + row*2)}) Patient {patient} - Peritoneal",
        ylim,
        ylabel=True
    )

for ax in axes[:,1]:
   ax.tick_params (axis="y",labelleft=True)

# ===== LEGEND =====
handles, labels = axes[0, 0].get_legend_handles_labels()

for ax in axes.flatten():
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()

fig.legend(
    handles,
    labels,
    title="Event type",
    loc="center right",
    bbox_to_anchor=(1.02, 0.85)
)

plt.tight_layout(rect=[0.03, 0.03, 0.88, 0.98])

# ===== SAVE =====
plt.savefig(
    "Figure_GOI_individual_comparisons_ovary_peritoneal_filtered_panel.png",
    dpi=600,
    bbox_inches="tight"
)
plt.show()
```
 <br>
 
Save the file by pressing `Ctrl+O`, then `Enter`, followed by `Ctrl+X`.

## 4. Run the script to create the plot
Run the script:
```
python3 check_goi_per_person.py
```

## 5. Output
The results will be saved in the analysis directory as:
```
Figure_GOI_individual_comparisons_ovary_peritoneal_filtered_panel.png
```
 <br>
 
The plot will be presented when it is done. 
