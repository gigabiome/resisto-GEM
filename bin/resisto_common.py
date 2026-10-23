"""Shared helpers for the resisto-GEM model-building scripts (reframed / CarveMe based)."""
import json
import os
import platform
from copy import deepcopy

import carveme
from carveme.reconstruction.utils import load_media_db
from reframed import load_cbmodel, save_cbmodel, set_default_solver
from reframed.core.model import Compartment

# ID of the generic biomass reaction in CarveMe 1.6.6 universe models, as read by reframed
# (verified: it carries no R_ prefix). Code asserts it via model.biomass_reaction.
CARVEME_BIOMASS = 'Growth'
# Metadata key used to tag gap-filled reactions; same key CarveMe uses in merge_models().
GAPFILL_TAG = 'GAP_FILL'
SOURCE_TAG = 'RESISTOGEM_SOURCE'


def rid(reaction_id):
    """Reaction ID in reframed/BiGG-flavour form (R_ prefix)."""
    return reaction_id if reaction_id.startswith('R_') else 'R_' + reaction_id


def universe_path(name):
    path = os.path.join(carveme.project_dir, 'data', 'generated', f'universe_{name}.xml.gz')
    if not os.path.exists(path):
        raise FileNotFoundError(f'CarveMe universe "{name}" not found at {path}')
    return path


def load_universe(name):
    return load_cbmodel(universe_path(name), flavor='bigg')


def load_medium(medium, mediadb=None):
    path = mediadb or os.path.join(carveme.project_dir, 'data', 'input', 'media_db.tsv')
    media = load_media_db(path)
    if medium not in media:
        raise KeyError(f'Medium "{medium}" not in {path}; available: {sorted(media)}')
    return sorted(media[medium]), path


def harmonise_compartments(model, target):
    """Rename compartments 'x' -> 'C_x' (CarveMe convention) in place.

    Compartments already present in `target` are reused so that metabolites shared
    with the CarveMe universe end up in the same compartment object. Returns the mapping.
    """
    mapping = {}
    for c_id, comp in list(model.compartments.items()):
        new_id = c_id if c_id.startswith('C_') else 'C_' + c_id
        mapping[c_id] = new_id
        if new_id in target.compartments:
            new_comp = deepcopy(target.compartments[new_id])
        else:
            new_comp = Compartment(new_id, name=comp.name, external=comp.external, size=comp.size)
        if new_id != c_id:
            del model.compartments[c_id]
        model.compartments[new_id] = new_comp
    for met in model.metabolites.values():
        met.compartment = mapping[met.compartment]
    return mapping


def copy_reaction_into(model, source, r_id, new_id=None):
    """Copy reaction r_id (plus missing metabolites, compartments and genes) from source into model.

    With new_id the copy is added under that ID, leaving any existing reaction r_id untouched.
    """
    rxn = source.reactions[r_id]
    for m_id in rxn.stoichiometry:
        if m_id not in model.metabolites:
            met = deepcopy(source.metabolites[m_id])
            if met.compartment not in model.compartments:
                model.add_compartment(deepcopy(source.compartments[met.compartment]))
            model.add_metabolite(met)
    new = deepcopy(rxn)
    if new_id:
        new.id = new_id
    for g_id in list(_gpr_genes(new)):
        if g_id not in model.genes and g_id in source.genes:
            model.add_gene(deepcopy(source.genes[g_id]))
    model.add_reaction(new)
    return new


def _gpr_genes(rxn):
    gpr = getattr(rxn, 'gpr', None)
    return gpr.get_genes() if gpr is not None else set()


def save_model(model, path):
    """Save as SBML (BiGG flavour, as CarveMe does) without boundaryCondition on any species.

    reframed writes boundaryCondition = compartment.external, and marks C_e as external when it
    loads CarveMe models. COBRApy then auto-creates EX_<met> reactions whose IDs collide with the
    model's own R_EX_<met> exchanges and opens their uptake (verified: growth 50.2 instead of 3.48
    on the CarveMe E. coli benchmark). Exchanges are explicit reactions, so no species needs the flag.
    """
    for comp in model.compartments.values():
        comp.external = False
    save_cbmodel(model, path, flavor='bigg')


def set_solver(name, threads=None):
    """Select the reframed solver and, for CPLEX/Gurobi, cap the thread count.

    reframed 1.6.0 exposes no thread parameter, so solver construction is wrapped.
    Without this, CPLEX uses every core on the node, which oversubscribes SLURM jobs.
    """
    set_default_solver(name)
    if threads is None:
        return
    threads = int(threads)
    if name == 'cplex':
        from reframed.solvers import cplex_solver
        cls = cplex_solver.CplexSolver
        original = cls.__init__

        def __init__(self, *args, **kwargs):
            original(self, *args, **kwargs)
            self.problem.parameters.threads.set(threads)
        cls.__init__ = __init__
    elif name == 'gurobi':
        from reframed.solvers import gurobi_solver
        cls = gurobi_solver.GurobiSolver
        original = cls.__init__

        def __init__(self, *args, **kwargs):
            original(self, *args, **kwargs)
            self.problem.setParam('Threads', threads)
        cls.__init__ = __init__


def solver_version(name):
    try:
        if name == 'cplex':
            import cplex
            return cplex.Cplex().get_version()
        if name == 'gurobi':
            import gurobipy
            return '.'.join(map(str, gurobipy.gurobi.version()))
        if name == 'scip':
            import pyscipopt
            return pyscipopt.__version__
    except Exception as e:  # report, do not hide
        return f'unavailable ({e})'
    return 'unknown'


def software_versions(solver=None):
    import pandas
    import reframed
    versions = {
        'python': platform.python_version(),
        'carveme': carveme.__version__,
        'reframed': getattr(reframed, '__version__', 'unknown'),
        'pandas': pandas.__version__,
    }
    if solver:
        versions[f'solver_{solver}'] = solver_version(solver)
    return versions


def read_fbc_fields(sbml_path):
    """{metabolite id: (formula or None, charge or None)} from the SBML FBC attributes.

    reframed reads formula/charge from the SBML notes when present (which can hold several
    ';'-joined formulas or a charge that disagrees with fbc:charge), and COBRApy turns a missing
    fbc:charge into 0. Reading libSBML directly keeps "not set" distinct from "0".
    """
    import libsbml
    doc = libsbml.readSBMLFromFile(sbml_path)
    model = doc.getModel()
    if model is None:
        raise IOError(f'Cannot read SBML model from {sbml_path}')
    fields = {}
    for species in model.getListOfSpecies():
        fbc = species.getPlugin('fbc')
        formula = fbc.getChemicalFormula() if fbc is not None and fbc.isSetChemicalFormula() else None
        charge = fbc.getCharge() if fbc is not None and fbc.isSetCharge() else None
        fields[species.getId()] = (formula or None, charge)
    return fields


def write_json(obj, path):
    with open(path, 'w') as f:
        json.dump(obj, f, indent=2, sort_keys=True, default=str)
