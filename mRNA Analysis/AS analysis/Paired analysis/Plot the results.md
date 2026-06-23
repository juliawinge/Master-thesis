## 1. Change to the correct directory
Go to the directory containing the patient-matched rMATS result:
```
cd /your_directory/sample_download/your_dataset/rMATS/results/paired_control
```
 <br>

## 2. Create the plotting script
Create a python script to visulize the filtered AS events:
```
nano check_goi_per_person.py
```
 <br>
 
Paste thr following code:
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

# ===== DATA ==========================================================
df = pd.read_csv("goi_rmats_summary.tsv", sep="\t")

events = ["SE", "A3SS", "A5SS", "MXE", "RI"]

genes = ["CYP4F11", "ENOSF1", "SLC1A1", "STAT3", "TCHH"]

# ===== COLORS =========================================================
event_colors = {
    "SE": "#DD8452",
    "A3SS": "#4C72B0",
    "A5SS": "#55A868",
    "MXE": "#B39FDD",
    "RI": "#D62728"
}

# ===== PATIENTS AND COMPARISONS ======================================
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

print("Combined ovary and peritoneal individual comparison panel figure created!")
```
<br>

Save the file by pressing `Ctrl+O`, then `Enter`, followed by `Ctrl+X`.

## 4. Run the plotting script
```
python3 check_goi_per_person.py
```

## 5. Output
The results will be saved in the analysis directory as:
`Figure_GOI_individual_comparisons_ovary_peritoneal_filtered_panel.png`

The figure summarizes filtered AS evenets in the selected GOIs for matched ectopic lesions samples. Each row represents one patient, and the columns show ovarian and preitoneial lesion comparisons. Bars represents GOIs, and stacked colors represent diffrent rMATS event types.
