"""Patch Ecoli-GEM v0.2.0 so it can grow; the original download stays untouched.

v0.2.0 makes HSDy (homoserine dehydrogenase) irreversible in the *oxidative*
direction (hom__L -> aspsa). That blocks Thr/Met/Ile synthesis, so biomass = 0.
The fix mirrors upstream SysBioChalmers/Ecoli-GEM#34 / PR #36 (open as of 2026-10):
keep HSDy irreversible, but in the biosynthetic direction (aspsa -> hom__L).

Usage: python scripts/patch_ecoli_gem.py
"""
import logging
from pathlib import Path

import cobra

logging.getLogger("cobra").setLevel(logging.ERROR)

d = Path(__file__).resolve().parent.parent / "reference_gems/Escherichia_coli/models/Ecoli-GEM_v0.2.0"
model = cobra.io.read_sbml_model(str(d / "Ecoli-GEM.xml"))
print(f"original growth: {model.slim_optimize():.3f}")

model.reactions.HSDy.bounds = (-1000.0, 0.0)
model.id = "Ecoli_GEM_v0_2_0_HSDyfix"

growth = model.slim_optimize()
print(f"patched growth:  {growth:.3f}")
assert growth > 0.1, "patched model still does not grow"

cobra.io.write_sbml_model(model, str(d / "Ecoli-GEM_HSDyfix.xml"))
cobra.io.save_json_model(model, str(d / "Ecoli-GEM_HSDyfix.json"))
