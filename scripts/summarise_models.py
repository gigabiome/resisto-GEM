"""Summarise every reference GEM: size, exchange reactions and biomass flux.

For each model this reports
  - number of reactions, metabolites and genes
  - number of exchange reactions
  - biomass (objective) flux with the model's own bounds ("as supplied")
  - biomass flux with every exchange unconstrained (bounds -1000..1000)

The objective reaction itself is never opened, even when it is an exchange
(e.g. iLP844 uses the biomass export EX_cpd11416_c0 as objective); otherwise the
model could import "biomass" and re-export it.

Writes summaries/reference_gem_summary.csv and summaries/reference_gem_summary.pdf.
Usage: ~/anaconda3/envs/gem/bin/python scripts/summarise_models.py
"""
import logging
from pathlib import Path

import cobra
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

logging.getLogger("cobra").setLevel(logging.ERROR)

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "reference_gems"
OUT = ROOT / "summaries"
BOUND = 1000.0

# (species, model, role, path relative to reference_gems)
MODELS = [
    ("Escherichia coli", "iML1515", "primary", "Escherichia_coli/models/iML1515.xml"),
    ("Escherichia coli", "Ecoli-GEM v0.2.0 (HSDy fix)", "primary alt.", "Escherichia_coli/models/Ecoli-GEM_v0.2.0/Ecoli-GEM_HSDyfix.xml"),
    ("Escherichia coli", "Ecoli-GEM v0.2.0 (as released)", "unpatched", "Escherichia_coli/models/Ecoli-GEM_v0.2.0/Ecoli-GEM.xml"),
    ("Staphylococcus aureus", "iYS854", "primary", "Staphylococcus_aureus/models/iYS854.xml"),
    ("Acinetobacter baumannii", "iACB23LX", "primary", "Acinetobacter_baumannii/models/Leonidou2024/iACB23LX/iACB23LX.xml"),
    ("Acinetobacter baumannii", "iCN718 (BiGG)", "secondary", "Acinetobacter_baumannii/models/iCN718.xml"),
    ("Acinetobacter baumannii", "iCN718 (Leonidou 2024)", "re-curated", "Acinetobacter_baumannii/models/Leonidou2024/iCN718/model.xml"),
    ("Acinetobacter baumannii", "iAB5075 (Leonidou 2024)", "re-curated", "Acinetobacter_baumannii/models/Leonidou2024/iAB5075/model.xml"),
    ("Acinetobacter baumannii", "iATCC1906 (Leonidou 2024)", "re-curated", "Acinetobacter_baumannii/models/Leonidou2024/iATCC1906/model.xml"),
    ("Acinetobacter baumannii", "iLP844 (Leonidou 2024)", "re-curated", "Acinetobacter_baumannii/models/Leonidou2024/iLP844/model.xml"),
    ("Acinetobacter baumannii", "iJS784 (Leonidou 2024)", "re-curated", "Acinetobacter_baumannii/models/Leonidou2024/iJS784/model.xml"),
    ("Acinetobacter baumannii", "AbyMBEL891 (Leonidou 2024)", "re-curated", "Acinetobacter_baumannii/models/Leonidou2024/AbyMBEL891/model.xml"),
    ("Acinetobacter baumannii", "AGORA AB0057 (Leonidou 2024)", "re-curated", "Acinetobacter_baumannii/models/Leonidou2024/AGORA/model.xml"),
    ("Pseudomonas aeruginosa", "iSD1509", "primary", "Pseudomonas_aeruginosa/models/iSD1509/model.xml"),
    ("Pseudomonas aeruginosa", "iPau21", "secondary", "Pseudomonas_aeruginosa/models/iPau21.xml"),
    ("Klebsiella pneumoniae", "KpSC-pan v2", "primary", "Klebsiella_pneumoniae/models/KpSC_pan_metabolic_model_v2.xml"),
    ("Klebsiella pneumoniae", "iYL1228", "secondary", "Klebsiella_pneumoniae/models/iYL1228.xml"),
    ("Klebsiella pneumoniae", "Jenior clinical", "context", "Klebsiella_pneumoniae/models/Jenior2022_context/clinical_maxfit/model.xml"),
    ("Klebsiella pneumoniae", "Jenior clinical (reps)", "context", "Klebsiella_pneumoniae/models/Jenior2022_context/clinical_maxfit_reps/model.xml"),
    ("Klebsiella pneumoniae", "Jenior laboratory", "context", "Klebsiella_pneumoniae/models/Jenior2022_context/laboratory_maxfit/model.xml"),
    ("Klebsiella pneumoniae", "Jenior laboratory (reps)", "context", "Klebsiella_pneumoniae/models/Jenior2022_context/laboratory_maxfit_reps/model.xml"),
    ("Enterococcus faecium", "iDR478", "primary", "Enterococcus_faecium/models/iDR478/Model_Files/DataSheet_S1_iDR478_model_SBML.xml"),
]

NOTES = {
    "Ecoli-GEM v0.2.0 (as released)": "HSDy blocks Thr/Met/Ile synthesis",
    "iCN718 (Leonidou 2024)": "ships with minimal medium (0.1)",
    "iAB5075 (Leonidou 2024)": "ships with minimal medium (0.1)",
    "iATCC1906 (Leonidou 2024)": "minimal medium (0.1); KEGG IDs",
    "iLP844 (Leonidou 2024)": "minimal medium (0.1); ModelSEED IDs",
    "iJS784 (Leonidou 2024)": "cannot make biomass; ModelSEED IDs",
    "AbyMBEL891 (Leonidou 2024)": "no GPRs; 1 compartment; R### IDs",
    "AGORA AB0057 (Leonidou 2024)": "minimal medium (0.1); VMH IDs",
    "iPau21": "all exchanges closed; ModelSEED IDs",
    "Jenior clinical": "fluxes fixed by transcriptome fit; pruned",
    "Jenior clinical (reps)": "fluxes fixed by transcriptome fit; pruned",
    "Jenior laboratory": "fluxes fixed by transcriptome fit; pruned",
    "Jenior laboratory (reps)": "fluxes fixed by transcriptome fit; pruned",
}


def exchange_reactions(model):
    """cobra's extracellular-compartment heuristic, falling back to EX_ prefixes."""
    ex = list(model.exchanges)
    return ex if ex else [r for r in model.reactions if r.id.startswith("EX_")]


def biomass_flux(model):
    sol = model.optimize()
    return sol.objective_value if sol.status == "optimal" else float("nan"), sol.status


def summarise(species, name, role, rel):
    model = cobra.io.read_sbml_model(str(REF / rel))
    objective = [r for r in model.reactions if r.objective_coefficient]
    exchanges = exchange_reactions(model)

    default_flux, _ = biomass_flux(model)
    with model:
        for r in exchanges:
            if r not in objective:
                r.bounds = (-BOUND, BOUND)
        open_flux, open_status = biomass_flux(model)

    return {
        "Species": species,
        "Model": name,
        "Role": role,
        "Reactions": len(model.reactions),
        "Metabolites": len(model.metabolites),
        "Genes": len(model.genes),
        "Exchanges": len(exchanges),
        "Objective": ", ".join(r.id for r in objective),
        "Biomass (as supplied)": default_flux,
        "Biomass (exchanges open)": open_flux,
        "Status (exchanges open)": open_status,
        "Notes": NOTES.get(name, ""),
        "Path": rel,
    }


def fmt(x):
    if pd.isna(x):
        return "infeasible"
    return f"{x:.3f}"


def write_pdf(df, path):
    cols = ["Species", "Model", "Role", "Reactions", "Metabolites", "Genes", "Exchanges",
            "Biomass (as supplied)", "Biomass (exchanges open)", "Notes"]
    headers = ["Species", "Model", "Role", "Reactions", "Metabolites", "Genes", "Exchange\nreactions",
               "Biomass flux\n(as supplied)", "Biomass flux\n(exchanges open)", "Notes"]
    cells = df[cols].copy()
    cells["Species"] = cells["Species"].str.replace(r"^(\w)\w+ ", r"\1. ", regex=True)
    for c in ["Biomass (as supplied)", "Biomass (exchanges open)"]:
        cells[c] = cells[c].map(fmt)

    fig = plt.figure(figsize=(11.69, 8.27))  # A4 landscape
    fig.text(0.03, 0.965, "Reference genome-scale metabolic models", fontsize=14, fontweight="bold", va="top")
    fig.text(0.03, 0.935, f"{len(df)} models, {df['Species'].nunique()} species", fontsize=9, color="#555555", va="top")
    ax = fig.add_axes([0.03, 0.12, 0.94, 0.79])
    ax.axis("off")
    table = ax.table(cellText=cells.values, colLabels=headers, cellLoc="center", bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(7.5)
    table.auto_set_column_width(list(range(len(cols))))

    shade = {}
    for i, sp in enumerate(df["Species"].unique()):
        shade[sp] = "#f2f2f2" if i % 2 else "white"
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#bbbbbb")
        if row == 0:
            cell.set_facecolor("#3b4a5a")
            cell.set_text_props(color="white", fontweight="bold")
            continue
        r = df.iloc[row - 1]
        cell.set_facecolor(shade[r["Species"]])
        if col in (0, 1, 9):
            cell.set_text_props(ha="left")
            cell.PAD = 0.03
        if r["Role"].startswith("primary") and col != 9:
            cell.set_text_props(fontweight="bold")
        if col in (7, 8) and not (cells.iloc[row - 1, col] != "infeasible" and float(cells.iloc[row - 1, col]) > 1e-6):
            cell.set_text_props(color="#b00020", fontweight="bold")

    note = (
        "Biomass flux = optimal value of each model's own objective (h$^{-1}$). 'As supplied' uses the bounds in the published file; "
        "'exchanges open' sets every exchange reaction to [-1000, 1000]\n"
        "(the objective reaction is never opened). Red = no growth. Bold rows = primary reference per species. "
        "Exchanges = cobrapy's extracellular boundary reactions (fallback: EX_ prefix).\n"
        "Fluxes with all exchanges open are not physiological; they show only whether the network can make biomass at all. "
        f"Solver: {cobra.Configuration().solver.__name__.split('.')[-1]}, cobrapy {cobra.__version__}."
    )
    fig.text(0.03, 0.10, note, fontsize=7, va="top", color="#333333", linespacing=1.5)
    with PdfPages(path) as pdf:
        pdf.savefig(fig)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    rows = []
    for spec in MODELS:
        print(f"loading {spec[1]} ...", flush=True)
        rows.append(summarise(*spec))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "reference_gem_summary.csv", index=False)
    write_pdf(df, OUT / "reference_gem_summary.pdf")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(df.drop(columns=["Path", "Objective"]).to_string(index=False))


if __name__ == "__main__":
    main()
