#!/usr/bin/env python
"""Run CarveMe's unmodified `carve` command with three additions needed on a cluster:

  1. CPLEX/Gurobi thread count capped at --threads (reframed exposes no thread setting, and the
     solvers otherwise use every core on the node).
  2. The carving MILP status/objective recorded in --status-json (CarveMe only prints failures).
  3. Non-zero exit if no model file was written: `carve` returns exit code 0 even when it prints
     "Failed to build model." or cannot run DIAMOND.

Usage: carve_wrapper.py --solver cplex --threads 4 --status-json s.json -- <carve arguments>
"""
import argparse
import os
import sys

from resisto_common import set_solver, software_versions, write_json


def main():
    if '--' not in sys.argv:
        sys.exit('Separate wrapper options from carve arguments with --')
    split = sys.argv.index('--')
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--solver', required=True)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--status-json', required=True)
    args = p.parse_args(sys.argv[1:split])
    carve_args = sys.argv[split + 1:]
    if '--solver' in carve_args:
        sys.exit('Pass --solver to the wrapper, not to carve')
    if '-o' not in carve_args:
        sys.exit('carve output (-o) is required')
    output = carve_args[carve_args.index('-o') + 1]

    set_solver(args.solver, args.threads)

    import carveme.reconstruction.carving as carving
    record = {'solver': args.solver, 'threads': args.threads, 'carve_args': carve_args,
              'software': software_versions(args.solver)}
    original = carving.minmax_reduction

    def minmax_reduction(*a, **kw):
        sol = original(*a, **kw)
        record['milp'] = {'status': sol.status.name, 'message': getattr(sol, 'message', None),
                          'objective': sol.fobj}
        return sol
    carving.minmax_reduction = minmax_reduction

    from carveme.cli import carve
    sys.argv = ['carve'] + carve_args + ['--solver', args.solver]
    carve.main()

    record['model_written'] = os.path.exists(output)
    write_json(record, args.status_json)
    if not record['model_written']:
        sys.exit(f'carve did not write {output}; see log above (MILP: {record.get("milp")})')


if __name__ == '__main__':
    main()
