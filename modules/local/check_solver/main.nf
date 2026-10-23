process CHECK_SOLVER {
    label 'process_single'
    label 'resisto'

    input:
    val solver

    output:
    path 'solver_check.json', emit: report

    script:
    """
    check_solver.py --solver ${solver} --threads ${task.cpus} --out solver_check.json
    """

    stub:
    """
    echo '{"ok": true, "stub": true}' > solver_check.json
    """
}
