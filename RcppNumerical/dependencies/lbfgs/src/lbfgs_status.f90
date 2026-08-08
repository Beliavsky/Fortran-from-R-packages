module lbfgs_status
    implicit none
    private

    integer, parameter, public :: lbfgs_success = 0
    integer, parameter, public :: lbfgs_convergence = 0
    integer, parameter, public :: lbfgs_stop = 1
    integer, parameter, public :: lbfgs_already_minimized = 2

    integer, parameter, public :: lbfgserr_unknownerror = -1024
    integer, parameter, public :: lbfgserr_logicerror = -1023
    integer, parameter, public :: lbfgserr_outofmemory = -1022
    integer, parameter, public :: lbfgserr_canceled = -1021
    integer, parameter, public :: lbfgserr_invalid_n = -1020
    integer, parameter, public :: lbfgserr_invalid_n_sse = -1019
    integer, parameter, public :: lbfgserr_invalid_x_sse = -1018
    integer, parameter, public :: lbfgserr_invalid_epsilon = -1017
    integer, parameter, public :: lbfgserr_invalid_testperiod = -1016
    integer, parameter, public :: lbfgserr_invalid_delta = -1015
    integer, parameter, public :: lbfgserr_invalid_linesearch = -1014
    integer, parameter, public :: lbfgserr_invalid_minstep = -1013
    integer, parameter, public :: lbfgserr_invalid_maxstep = -1012
    integer, parameter, public :: lbfgserr_invalid_ftol = -1011
    integer, parameter, public :: lbfgserr_invalid_wolfe = -1010
    integer, parameter, public :: lbfgserr_invalid_gtol = -1009
    integer, parameter, public :: lbfgserr_invalid_xtol = -1008
    integer, parameter, public :: lbfgserr_invalid_maxlinesearch = -1007
    integer, parameter, public :: lbfgserr_invalid_orthantwise = -1006
    integer, parameter, public :: lbfgserr_invalid_orthantwise_start = -1005
    integer, parameter, public :: lbfgserr_invalid_orthantwise_end = -1004
    integer, parameter, public :: lbfgserr_outofinterval = -1003
    integer, parameter, public :: lbfgserr_incorrect_tminmax = -1002
    integer, parameter, public :: lbfgserr_rounding_error = -1001
    integer, parameter, public :: lbfgserr_minimumstep = -1000
    integer, parameter, public :: lbfgserr_maximumstep = -999
    integer, parameter, public :: lbfgserr_maximumlinesearch = -998
    integer, parameter, public :: lbfgserr_maximumiteration = -997
    integer, parameter, public :: lbfgserr_widthtoosmall = -996
    integer, parameter, public :: lbfgserr_invalidparameters = -995
    integer, parameter, public :: lbfgserr_increasegradient = -994

    integer, parameter, public :: lbfgs_linesearch_default = 0
    integer, parameter, public :: lbfgs_linesearch_morethuente = 0
    integer, parameter, public :: lbfgs_linesearch_backtracking_armijo = 1
    integer, parameter, public :: lbfgs_linesearch_backtracking = 2
    integer, parameter, public :: lbfgs_linesearch_backtracking_wolfe = 2
    integer, parameter, public :: lbfgs_linesearch_backtracking_strong_wolfe = 3

    public :: lbfgs_status_message

contains

    pure function lbfgs_status_message(status) result(message)
        integer, intent(in) :: status
        character(len=:), allocatable :: message

        select case (status)
        case (lbfgs_success)
            message = "L-BFGS reached convergence."
        case (lbfgs_stop)
            message = "The delta-based stopping criterion was satisfied."
        case (lbfgs_already_minimized)
            message = "The initial variables already minimize the objective function."
        case (lbfgserr_unknownerror)
            message = "Unknown error."
        case (lbfgserr_logicerror)
            message = "Logic error."
        case (lbfgserr_outofmemory)
            message = "Insufficient memory."
        case (lbfgserr_canceled)
            message = "The minimization process was canceled by the progress callback."
        case (lbfgserr_invalid_n)
            message = "Invalid number of variables."
        case (lbfgserr_invalid_n_sse)
            message = "Invalid SSE-padded number of variables."
        case (lbfgserr_invalid_x_sse)
            message = "Invalid SSE memory alignment."
        case (lbfgserr_invalid_epsilon)
            message = "Invalid epsilon parameter."
        case (lbfgserr_invalid_testperiod)
            message = "Invalid past parameter."
        case (lbfgserr_invalid_delta)
            message = "Invalid delta parameter."
        case (lbfgserr_invalid_linesearch)
            message = "Invalid line-search parameter."
        case (lbfgserr_invalid_minstep)
            message = "Invalid minimum line-search step."
        case (lbfgserr_invalid_maxstep)
            message = "Invalid maximum line-search step."
        case (lbfgserr_invalid_ftol)
            message = "Invalid ftol parameter."
        case (lbfgserr_invalid_wolfe)
            message = "Invalid Wolfe parameter."
        case (lbfgserr_invalid_gtol)
            message = "Invalid gtol parameter."
        case (lbfgserr_invalid_xtol)
            message = "Invalid xtol parameter."
        case (lbfgserr_invalid_maxlinesearch)
            message = "Invalid maximum line-search evaluation count."
        case (lbfgserr_invalid_orthantwise)
            message = "Invalid OWL-QN coefficient."
        case (lbfgserr_invalid_orthantwise_start)
            message = "Invalid OWL-QN start index."
        case (lbfgserr_invalid_orthantwise_end)
            message = "Invalid OWL-QN end index."
        case (lbfgserr_outofinterval)
            message = "The line-search step left the interval of uncertainty."
        case (lbfgserr_incorrect_tminmax)
            message = "The line-search interval is invalid."
        case (lbfgserr_rounding_error)
            message = "Rounding error prevented further line-search progress."
        case (lbfgserr_minimumstep)
            message = "The line-search step became smaller than min_step."
        case (lbfgserr_maximumstep)
            message = "The line-search step became larger than max_step."
        case (lbfgserr_maximumlinesearch)
            message = "The line search reached its maximum number of evaluations."
        case (lbfgserr_maximumiteration)
            message = "The optimizer reached its maximum number of iterations."
        case (lbfgserr_widthtoosmall)
            message = "The line-search uncertainty interval became too narrow."
        case (lbfgserr_invalidparameters)
            message = "Invalid line-search parameters or nonpositive trial step."
        case (lbfgserr_increasegradient)
            message = "The search direction is not a descent direction."
        case default
            message = "Unrecognized L-BFGS status code."
        end select
    end function lbfgs_status_message

end module lbfgs_status
