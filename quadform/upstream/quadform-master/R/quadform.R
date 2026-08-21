#' @export
`ht` <- function(x){ t(Conj(x)) }

#' @export
`cprod` <- function(x, y=NULL){
    if(is.null(y)){
        return(crossprod(Conj(x), x))
    } else {
        return(crossprod(Conj(x), y))
    }
}

#' @export
`tcprod` <- function(x, y=NULL){
    if(is.null(y)){
        return(tcrossprod(x, Conj(x)))
    } else {
        return(tcrossprod(x, Conj(y)))
    }
}

#' @export
`quad.form.chol` <- function (chol, x){
        jj <- cprod(chol, x)
        cprod(jj, jj)
}

#' @export
`quad.form` <- function (M, x){ crossprod(crossprod(M, Conj(x)), x) }

#' @export
`quad.form.inv` <- function (M, x){ cprod(x, solve(M, x)) }

#' @export
`quad3.form_ab` <- function(M, left, right){ crossprod(crossprod(M, Conj(left)), right) }

#' @export
`quad3.form_bc` <- function(M, left, right){ cprod(left, (M %*% right)) }

#' @export
`quad3.form` <- function(M,left,right){
    left <- as.matrix(left)
    right <- as.matrix(right)
    if(ncol(left) < ncol(right)){
        quad3.form_ab(M, left, right)
    } else {
        quad3.form_bc(M, left, right)
    }
}

#' @export
`quad3.form.inv` <- function(M, left, right){ cprod(left, solve(M, right)) }

#' @export
`quad3.tform_ab` <- function(M, left, right){ tcprod(left %*% M,right)}
#' @export

`quad3.tform_bc` <- function(M, left, right){ tcrossprod(left, tcrossprod(Conj(right), M)) }

#' @export
`quad3.tform` <- function(M, left, right){
    if(nrow(left) < nrow(right)){
        quad3.tform_ab(M, left, right)
    } else {
        quad3.tform_bc(M, left, right)
    }
}

#' @export
`quad.tform` <- function(M, x){ tcrossprod(x, tcrossprod(Conj(x), M)) }

#' @export
`quad.tform.inv` <- function(M, x){ quad.form.inv(M, ht(x)) }

#' @export
`quad.diag` <- function(M, x){ colSums(crossprod(M, Conj(x)) * x) }

#' @export
`quad.tdiag` <- function(M, x){ rowSums(tcrossprod(Conj(x), M) * x) }

#' @export
`quad3.diag` <- function(M, left, right){ colSums(crossprod(M, Conj(left)) * right) }

#' @export
`quad3.tdiag` <- function(M, left, right){ colSums(t(left) * tcprod(M, right)) }

#' @export
`quad.trace` <- function(M, x){ sum(crossprod(M, Conj(x)) * x) }

#' @export
`quad.ttrace` <- function(M, x){ sum(tcrossprod(Conj(x), M) * x) }

#' @export
cp <- cprod

#' @export
tcp <- tcprod

#' @export
qf <- quad.form

#' @export
qfi <- quad.form.inv

#' @export
q3 <- quad3.form

#' @export
q3i <- quad3.form.inv


#' @export
q3t <- quad3.tform

#' @export
qt <- quad.tform

#' @export
qti <- quad.tform.inv

#' @export
qd <- quad.diag

#' @export
qtd <- quad.tdiag

#' @export
q3d <- quad3.diag

#' @export
q3td <- quad3.tdiag


#' @export
qtr <- quad.trace

#' @export
qttr <- quad.ttrace

