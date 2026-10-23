#!/usr/bin/env python
"""Check that the requested MILP solver works at the problem sizes CarveMe needs.

CarveMe's carving and gap-filling MILPs have thousands of binary variables. Free/community
editions of commercial solvers fail on these (CPLEX Community Edition: 1,000 variables and
constraints; Gurobi's pip "restricted" licence: 2,000), which would otherwise only show up as
failed reconstructions. This builds a 2,500-variable binary problem and checks the optimum.
"""
import argparse
import sys

from reframed.solvers import solver_instance
from reframed.solvers.solver import VarType

from resisto_common import set_solver, software_versions, write_json

N = 2500


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver', required=True)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--out', default='solver_check.json')
    args = p.parse_args()

    record = {'solver': args.solver, 'threads': args.threads, 'variables': N,
              'software': software_versions(args.solver)}
    try:
        set_solver(args.solver, args.threads)
        solver = solver_instance()
        names = [f'x{i}' for i in range(N)]
        for name in names:
            solver.add_variable(name, 0, 1, vartype=VarType.BINARY)
        solver.update()
        solver.add_constraint('cap', {n: 1 for n in names}, '<', 10)
        for i in range(0, N - 1, 2):
            solver.add_constraint(f'pair{i}', {names[i]: 1, names[i + 1]: 1}, '<', 1)
        solver.update()
        sol = solver.solve({n: 1 for n in names}, minimize=False)
        record['status'] = sol.status.name
        record['objective'] = sol.fobj
        record['ok'] = sol.status.name == 'OPTIMAL' and abs((sol.fobj or 0) - 10) < 1e-6
    except Exception as e:
        record['ok'] = False
        record['error'] = f'{type(e).__name__}: {e}'

    write_json(record, args.out)
    if not record['ok']:
        sys.exit(f'Solver check failed for {args.solver}: {record.get("error") or record}')
    print(f'{args.solver} OK ({record["software"].get("solver_" + args.solver)})')


if __name__ == '__main__':
    main()
