#' Solve the elliptic curve \eqn{y(x)}
#'
#' Solve the elliptic curve \eqn{y(x)} by constructing
#' a cubic polynomial from ecd object. Then solve it
#' and take the smallest real root.
#'
#' @method solve ecd
#'
#' @param a An object of ecd class
#' @param b A vector of \eqn{x} values
#' @param ... Not used. Only here to match the generic signature.
#'
#' @return A vector of roots for \eqn{y(x)}
#'
#' @keywords polynomial solve
#'
#' @examples
#' d <- ecd()
#' x <- seq(-100,100,by=0.1)
#' y <- solve(d,x)

### <======================================================================>
#'
#' @export ecd solve
#' @export solve.ecd
#' @export
"solve.ecd" <- function(a, b, ...)
{

    # handle cusp
    cusp <- a@cusp

    # allow cusp error <= 1e-7, validated for alpha up to -100
    # this is okay since cusp should have only single real root everywhere on x.
    eps <- ifelse(cusp > 0, ifelse(a@alpha<100,1e-7,1e-7*log(a@alpha)^2), 0)

    get_real <- function(roots) {
        real_roots <- roots[abs(Im(roots)) <= eps]
        min(Re(real_roots))
    }

    # handle lambda model
    if ( a@lambda != 3 ) {
        if ( a@beta == 0 & a@alpha != 0 & a@gamma == 0 ) {
            # this is +/-y^lambda + alpha = xi^2
            return(solve_sym.ecd(a,b))
        }

        if ( a@alpha != 0 | a@gamma != 0 ) stop("Lambda-model must have zero alpha and gamma!")

        if ( a@beta == 0 ) {
			xi <- (b-a@mu)/a@sigma
			return(-(xi^2)^(1/a@lambda))
		} else {
            # this is handling the cases outside of analytic coverage in ecld.solve()
            # e.g. lambda = 2,4 have already been handled there

            #return(ecld.solve_by_poly(a,b))
            return(ecld.solve_isomorphic(a,b))
        }
    }

	# standard cubic case

    all_roots <- ecd.cubic(a, b)
    simplify2array(lapply(all_roots, get_real))
}
### <---------------------------------------------------------------------->
#' @rdname solve.ecd
setMethod("solve", c("ecd"), solve.ecd)
### <---------------------------------------------------------------------->
