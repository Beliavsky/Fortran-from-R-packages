## ROI plugin: qpOASES
## based on ecos interface

make_qpoases_signatures <- function()
    ROI_plugin_make_signature( objective = c("Q", "L"),
                               constraints = c("X", "L"),
                               types = c("C"),
                               bounds = c("X", "V"),
                               cones = c("X"),
                               maximum = c(TRUE, FALSE) )

## SOLVER CONTROLS
.add_controls <- function(solver) {
    ## qpOASES
    ROI_plugin_register_solver_control( solver, "printLevel", "verbose" )

    ROI_plugin_register_solver_control( solver, "hessian_type", "X" )
    ROI_plugin_register_solver_control( solver, "max_num_wsr", "X" )
    ROI_plugin_register_solver_control( solver, "cputime", "X" )

    ROI_plugin_register_solver_control( solver, "enableRamping", "X" )
    ROI_plugin_register_solver_control( solver, "enableFarBounds", "X" )
    ROI_plugin_register_solver_control( solver, "enableFlippingBounds", "X" )
    ROI_plugin_register_solver_control( solver, "enableRegularisation", "X" )
    ROI_plugin_register_solver_control( solver, "enableFullLITests", "X" )
    ROI_plugin_register_solver_control( solver, "enableNZCTests", "X" )
    ROI_plugin_register_solver_control( solver, "enableDriftCorrection", "X" )
    ROI_plugin_register_solver_control( solver, "enableCholeskyRefactorisation", "X" )
    ROI_plugin_register_solver_control( solver, "enableEqualities", "X" )
    ROI_plugin_register_solver_control( solver, "terminationTolerance", "X" )
    ROI_plugin_register_solver_control( solver, "boundTolerance", "X" )
    ROI_plugin_register_solver_control( solver, "boundRelaxation", "X" )
    ROI_plugin_register_solver_control( solver, "maxPrimalJump", "X" )
    ROI_plugin_register_solver_control( solver, "epsDen", "X" )
    ROI_plugin_register_solver_control( solver, "maxDualJump", "X" )
    ROI_plugin_register_solver_control( solver, "initialRamping", "X" )
    ROI_plugin_register_solver_control( solver, "finalRamping", "X" )
    ROI_plugin_register_solver_control( solver, "initialFarBounds", "X" )
    ROI_plugin_register_solver_control( solver, "growFarBounds", "X" )
    ROI_plugin_register_solver_control( solver, "rcondSMin", "X" )
    ROI_plugin_register_solver_control( solver, "epsFlipping", "X" )
    ROI_plugin_register_solver_control( solver, "epsRegularisation", "X" )
    ROI_plugin_register_solver_control( solver, "epsIterRef", "X" )
    ROI_plugin_register_solver_control( solver, "epsLITests", "X" )
    ROI_plugin_register_solver_control( solver, "epsNZCTests", "X" )
    ROI_plugin_register_solver_control( solver, "numRegularisationSteps", "X" )
    ROI_plugin_register_solver_control( solver, "numRefinementSteps", "X" )
    ROI_plugin_register_solver_control( solver, "dropBoundPriority", "X" )
    ROI_plugin_register_solver_control( solver, "dropEqConPriority", "X" )
    ROI_plugin_register_solver_control( solver, "dropIneqConPriority", "X" )
    ROI_plugin_register_solver_control( solver, "enableInertiaCorrection", "X" )
    ROI_plugin_register_solver_control( solver, "enableDropInfeasibles", "X" )
    ROI_plugin_register_solver_control( solver, "initialStatusBounds", "X" )

    invisible( TRUE )
}

.onLoad <- function( libname, pkgname ) {
    solver <- "qpoases"
    ## Solver plugin name (based on package name)
    if( ! pkgname %in% ROI_registered_solvers() ){
        ## Register solver methods here.
        ## One can assign several signatures a single solver method
        ROI_plugin_register_solver_method(
            signatures = make_qpoases_signatures(),
            solver = solver,
            method = getFunction( "solve_OP", where = getNamespace(pkgname)) )
        ## Finally, for status code canonicalization add status codes to data base
        .add_status_codes( solver )
        .add_controls( solver )
    }
}

