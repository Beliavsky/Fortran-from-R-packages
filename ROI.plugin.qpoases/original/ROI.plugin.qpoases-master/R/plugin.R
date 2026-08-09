## -----------------------------------------------------------------------------
##
##   qpOases:
##
## -----------------------------------------------------------------------------

to_dense_vector <- function(x, len, default = 0L) {
    y <- rep.int(default, len)
    if ( is.null(x$ind) ) return(y)
    y[x$ind] <- x$val
    return(y)
}

is_slam_zero_matrix <- function(x) length(x$i) == 0L

is_slam_identity_matrix <- function(x) {
    if ( nrow(x) != sum(x$v == 1) )
        return(FALSE)
    if ( all(seq_len(nrow(x)) != x$i) )
        return(FALSE)
    if ( all(seq_len(nrow(x)) != x$j) )
        return(FALSE)
    TRUE
}

#define HST_ZERO             0
#define HST_IDENTITY         1
#define HST_POSDEF           2
#define HST_POSDEF_NULLSPACE 3
#define HST_SEMIDEF          4
#define HST_INDEF            5
#define HST_UNKNOWN          6

## control <- list("printLevel" = 1L, "enableRamping" = 1L)
## 
## control <- list("hessian_type" = 1L)
## control <- list("hessian_type" = 2L)

qpoases_solve <- function(nvariables, nconstraints, hessian_type, 
    H, g, A, lb, ub, lbA, ubA, max_num_working_set_recalculations, cputime, 
    control) {
    m <- qproblem(nvariables, nconstraints, hessian_type)
    set_options(m, control)
    status <- init_qproblem(m, H, g, A, lb, ub, lbA, ubA, 
                            max_num_working_set_recalculations, cputime)
    msg <- list()
    msg$primal_solution <- get_primal_solution(m)
    msg$dual_solution <- get_dual_solution(m)
    msg$status <- status
    msg
}

solve_OP <- function(x, control = list()){
    solver <- "qpoases"

    nvariables <- ncol(constraints(x))
    nconstraints <- nrow(constraints(x))

    Q <- terms(objective(x))[["Q"]]
    L <- terms(objective(x))[["L"]]
    if ( maximum(x) ) {
        Q <- -Q
        L <- -L
    }
    if ( is.null(control$hessian_type) ) {
        if ( is_slam_zero_matrix(Q) ) { 
           hessian_type <- 0L
        } else if ( is_slam_identity_matrix(Q) ) {
            hessian_type <- 1L
        } else {
            hessian_type <- 6L
        }    
    } else {
        hessian_type <- control$hessian_type
        stopifnot(isTRUE(hessian_type %in% c(seq(0, 6))))
    }

    if ( hessian_type == 1L ) {
        control$enableRegularisation <- 0L
        control$numRegularisationSteps <- 0L
    }

    cputime <- if ( is.null(control$cputime) ) 100 else control$cputime
    
    ## m <- qproblem(nvariables, nconstraints, hessian_type)

    cntrl <- default_control()
    nam <- intersect(names(cntrl), names(control))
    cntrl[nam] <- control[nam]
    ## set_options(m, cntrl)

    ## print_options(m)
    H <- as.vector(t(Q))
    g <- as.vector(L)
    A <- as.vector(t(constraints(x)$L))
    lb <- to_dense_vector(bounds(x)$lower, nvariables)
    ub <- to_dense_vector(bounds(x)$upper, nvariables, Inf)
    lbA <- rep.int(-Inf, nconstraints)
    ubA <- rep.int( Inf, nconstraints)

    rhs <- constraints(x)$rhs
    i <- which(constraints(x)$dir == "==")
    lbA[i] <- rhs[i]
    ubA[i] <- rhs[i]

    i <- which(constraints(x)$dir == "<=")
    ubA[i] <- rhs[i]

    i <- which(constraints(x)$dir == ">=")
    lbA[i] <- rhs[i]
 
    max_num_working_set_recalculations <- 2000

    m <- list(qpoases_solve, nvariables = nvariables, nconstraints = nconstraints, 
              hessian_type = hessian_type, H = H, g = g, A = A,
              lb = lb, ub = ub, lbA = lbA, ubA = ubA, 
              max_num_working_set_recalculations = max_num_working_set_recalculations, 
              cputime = cputime, control = cntrl)
    mode(m) <- "call"

    if ( isTRUE(control$dry_run) ) {
        return(m)
    }

    msg <- eval(m)
    ## NOTE: There as been a issue that I needed to solve the
    ##       problem twice till the correct solution was returned.
    ## status <- init_qproblem(m, H, g, A, lb, ub, lbA, ubA, 
    ##                         max_num_working_set_recalculations, cputime)
    msg$objval <- objective(x)(msg$primal_solution)

    ROI_plugin_canonicalize_solution( solution = msg$primal_solution, optimum  = msg$objval,
                                      status   = msg$status, solver   = solver, message = msg )
}

## just
default_control <- function() {
    EPS <- 2.221e-16
    list("printLevel" = 1L,

         "enableRamping" = 1L,
         "enableFarBounds" = 1L,
         "enableFlippingBounds" = 1L,
         "enableRegularisation" = 1L,
         "enableFullLITests" = 0L,
         "enableNZCTests" = 1L,
         "enableDriftCorrection" = 1L,
         "enableCholeskyRefactorisation" = 0L,
         "enableEqualities" = 0L,

         "terminationTolerance" = 5.0e6 * EPS,
         "boundTolerance" = 1.0e6 * EPS,
         "boundRelaxation" = 1.0e4,

         "epsNum" = -2.221000e-13,
         "epsDen" =  2.221000e-13,
         "maxPrimalJump" = 1.0e8,
         "maxDualJump" = 1.0e8,

         "initialRamping" = 0.5,
         "finalRamping" = 1.0,
         "initialFarBounds" = 1.0e6,
         "growFarBounds" = 1.0e3,

         "rcondSMin" = 1.0e-14,

         "epsFlipping" = 1.0e3 * EPS,
         "epsRegularisation" = 1.0e3 * EPS,
         "epsIterRef" = 1.0e2 * EPS,
         "epsLITests" = 1.0e5 * EPS,
         "epsNZCTests" = 3.0e3 * EPS,

         "numRegularisationSteps" = 1L,
         "numRefinementSteps" = 1L,
    
         "dropBoundPriority" = 1L,
         "dropEqConPriority" = 1L,
         "dropIneqConPriority" = 1L,

         "enableInertiaCorrection" = 1L,
         "enableDropInfeasibles" = 0L,
    
         "initialStatusBounds" = -1L)
}

## printLevel is a value in -2:3
## qpoases_control <- function(eps = 2.221e-16, printLevel = 1L, 
##     enableRamping = 1L, enableFarBounds = 1L, enableFlippingBounds = 1L,
##     enableRegularisation = 1L, enableFullLITests = 0L, enableNZCTests = 1L,
##     enableDriftCorrection = 1L, enableCholeskyRefactorisation = 0L,
##     enableEqualities = 0L,
## 
##          terminationTolerance = 5.0e6 * EPS,
##          boundTolerance = 1.0e6 * EPS,
##          boundRelaxation = 1.0e4,
## 
##          epsNum = -2.221000e-13,
##          epsDen =  2.221000e-13,
##          maxPrimalJump = 1.0e8,
##          maxDualJump = 1.0e8,
## 
##          initialRamping = 0.5,
##          finalRamping = 1.0,
##          initialFarBounds = 1.0e6,
##          growFarBounds = 1.0e3,
## 
##          rcondSMin = 1.0e-14,
## 
##          epsFlipping = 1.0e3 * EPS,
##          epsRegularisation = 1.0e3 * EPS,
##          epsIterRef = 1.0e2 * EPS,
##          epsLITests = 1.0e5 * EPS,
##          epsNZCTests = 3.0e3 * EPS,
## 
##          numRegularisationSteps = 1L,
##          numRefinementSteps = 1L,
## 
##          dropBoundPriority = 1L,
##          dropEqConPriority = 1L,
##          dropIneqConPriority = 1L,
## 
##          enableInertiaCorrection = 1L,
##          enableDropInfeasibles = 0L,
## 
##          initialStatusBounds = -1L)
## }

## STATUS CODES
.add_status_codes <- function(solver) {
    for (i in seq_along(qpoases_status_codes)) {
        status_code <- qpoases_status_codes[[i]]
        success <- c("SUCCESSFUL_RETURN", "RET_QP_SOLVED")
        roi_code <- if ( status_code[["symbol"]] %in% success ) 0L else 1L
        ROI_plugin_add_status_code_to_db(solver, as.integer(status_code[["id"]]), 
                                         status_code[["symbol"]], status_code[["msg"]], 
                                         roi_code)
    }
    invisible(TRUE)
}
