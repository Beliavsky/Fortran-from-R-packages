#' Solve a quadratic program with absolute values in constraints & objective
#'
#' @description solveQPXT  allows for absolute value constraints and absolute values in the
#' objective.  buildQP builds a parameter list that can then be passed to
#' quadprog::solve.QP.compact or quadprog::solve.QP directly if desired by the user.
#' solveQPXT by default implicitly takes advantage of sparsity in the constraint
#' matrix and can improve numerical stability by normalizing the constraint matrix. For the
#' rest of the documentation, assume that Dmat is n x n.\cr
#' 
#' The solver solves the following problem (each * corresponds to matrix multiplication):
#' \preformatted{
#' min:
#' -t(dvec) * b + 1/2 t(b) * Dmat * b +
#' -t(dvecPosNeg) * c(b_positive, b_negative) +
#' -t(dvecPosNegDelta) * c(deltab_positive, deltab_negative)
#' 
#' s.t.
#' t(Amat) * b >= bvec 
#' t(AmatPosNeg) * c(b_positive, b_negative) >= bvecPosNeg
#' t(AmatPosNegDelta) * c(deltab_positive, deltab_negative) >= bvecPosNegDelta
#' b_positive, b_negative >= 0,
#' b = b_positive - b_negative
#' deltab_positive, deltab_negative >= 0,
#' b - b0 = deltab_positive - deltab_negative
#' }
#' 
#' @inheritParams quadprog::solve.QP
#'
#' 
#' @param AmatPosNeg 2n x k matrix of constraints on the positive and negative part of b
#' @param bvecPosNeg k length vector of thresholds to the constraints in AmatPosNeg
#' @param dvecPosNeg k * 2n length vector of loadings on the positive and negative part of b, respectively
#' 
#' @param b0 a starting point that describes the 'current' state of the problem such that
#' constraints and penalty on absolute changes in the decision variable from a starting point can
#' be incorporated.  b0 is an n length vector. Note that b0 is NOT a starting point for the
#' optimization - that is handled implicitly by quadprog.
#' @param AmatPosNegDelta 2n x l matrix of constraints on the positive and negative part of a change in b from a starting point, b0.
#' @param bvecPosNegDelta l length vector of thresholds to the constraints in AmatPosNegDelta
#' @param dvecPosNegDelta l * 2n length vector of loadings in the objective function on the positive and negative part of changes in b from a starting point of b0.
#' @param tol tolerance along the diagonal of the expanded Dmat for slack variables
#' @param compact logical: if TRUE, it is assumed that we want to use solve.QP.compact to solve the problem, which handles sparsity.
#' @param normalize logical: should constraint matrix be normalized
#' @param ... parameters to pass to buildQP when calling solveQPXT
#' 
#' @details In order to handle constraints on b_positive and b_negative, slack variables are introduced.  The total number of parameters in the problem increases by the following amounts: \cr
#' If all the new parameters (those not already used by quadprog) remain NULL, the problem size does not increase and quadprog::solve.QP (.compact) is called after normalizing the constraint matrix and converting to a sparse matrix representation by default.\cr
#' If AmatPosNeg, bvecPosNeg or dvecPosNeg are not null, the problem size increases by n
#' If AmatPosNegDelta or devecPosNegDelta are not null, the problem size increases by n.
#' This results in a potential problem size of up to 3 * n.
#' Despite the potential large increases in problem size, the underlying solver is written in
#' Fortran and converges quickly for problems involving even hundreds of parameters.  Additionally,
#' it has been the author's experience that solutions solved via the convex quadprog are much more
#' stable than those solved by other methods (e.g. a non-linear solver).
#'
#' Note that due to the fact that the constraints are by default normalized, the original constraint values the user passed will may not be returned by buildQP. 
#' @export solveQPXT
#'
#' @examples
#' ##quadprog example"
#' Dmat       <- matrix(0,3,3)
#' diag(Dmat) <- 1
#' dvec       <- c(0,5,0)
#' Amat       <- matrix(c(-4,-3,0,2,1,0,0,-2,1),3,3)
#' bvec       <- c(-8,2,0)
#' qp <- quadprog::solve.QP(Dmat,dvec,Amat,bvec=bvec)
#' qpXT <- solveQPXT(Dmat,dvec,Amat,bvec=bvec)
#' range(qp$solution - qpXT$solution)
#' 
#' N <- 10
#' set.seed(2)
#' cr <- matrix(runif(N * N, 0, .05), N, N)
#' diag(cr) <- 1
#' cr <- (cr + t(cr)) / 2
#' set.seed(3)
#' sigs <- runif(N, min = .02, max = .25)
#' set.seed(5)

#' dvec <- runif(N, -.1, .1)
#' Dmat <- sigs %o% sigs * cr
#' Amat <- cbind(diag(N), diag(N) * -1)
#' bvec <- c(rep(-1, N), rep(-1, N))

#' resBase <- solveQPXT(Dmat, dvec, Amat, bvec)
#' ##absolute value constraint on decision variable:
#' res <- solveQPXT(Dmat, dvec, Amat, bvec,
#' AmatPosNeg = matrix(rep(-1, 2 * N)), bvecPosNeg = -1)
#' sum(abs(res$solution[1:N]))
#'
#' ## penalty of L1 norm
#' resL1Penalty <- solveQPXT(Dmat, dvec, Amat, bvec, dvecPosNeg = -.005 * rep(1, 2 * N))
#' sum(abs(resL1Penalty$solution[1:N]))
#'
#' ## constraint on amount decision variable can vary from a starting point
#' b0 <- rep(.15, N)
#' thresh <- .25
#' res <- solveQPXT(Dmat, dvec, Amat, bvec, b0 = b0,
#' AmatPosNegDelta = matrix(rep(-1, 2 * N)), bvecPosNegDelta = -thresh)
#' sum(abs(res$solution[1:N] - b0))
#'
#' ##use buildQP, then call solve.QP.compact directly
#' qp <- buildQP(Dmat, dvec, Amat, bvec, b0 = b0,
#' AmatPosNegDelta = matrix(rep(-1, 2 * N)), bvecPosNegDelta = -thresh)
#' res2 <- do.call(quadprog::solve.QP.compact, qp)
#' range(res$solution - res2$solution)

solveQPXT <- function(...){                      
    
    qpArgs <- do.call(buildQP, list(...))
    
    if(!is.null(qpArgs$Aind)){
        res <- do.call(quadprog::solve.QP.compact, qpArgs)
    }else{
        res <- do.call(quadprog::solve.QP, qpArgs)
    }
    
    return(res)
}

#' @rdname solveQPXT
#' @export
buildQP <- function(Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE,
                    AmatPosNeg = NULL,
                    bvecPosNeg = NULL,
                    dvecPosNeg = NULL,
                    b0 = NULL,
                    AmatPosNegDelta = NULL,
                    bvecPosNegDelta = NULL,
                    dvecPosNegDelta = NULL,
                    tol = 1e-8,
                    compact = TRUE,
                    normalize = TRUE
                    ){
    
    ##number of decision variables
    N <- nrow(Dmat)
    ##M is equal length(|b|) or 0; L = length(|b - b0|) or 0
    M <- L <- 0
    if(!is.null(AmatPosNeg) || !is.null(dvecPosNeg)){
        M <- N
    }

    if(!is.null(AmatPosNegDelta) || !is.null(dvecPosNegDelta)){
        L <- N
    }
    
    if(is.null(Amat)){
        Amat <- matrix(0, N, 0)
    }
    if(is.null(AmatPosNeg)){
        AmatPosNeg <- matrix(0, 2 * M, 0)
    }
    if(is.null(AmatPosNegDelta)){
        AmatPosNegDelta <- matrix(0, 2 * L, 0)
    }
    
    ##number of constraints
    K <- ncol(Amat)
    K1 <- ncol(AmatPosNeg)
    K2 <- ncol(AmatPosNegDelta)

    ##expand original constraints to problem size
    Amat <- rbind(
        Amat,
        matrix(0, M + L, K)
    )

    ##create slack constraints: decision vector is: [b, |b|, |b - b0|]
    AmatSlack <- cbind(
        rbind(
            diag(x = 1, N, M),
            diag(x = 1, M, M),
            matrix(0, L, M)
        ),
        rbind(
            diag(x = -1, N, M),
            diag(x = 1, M, M),
            matrix(0, L, M)
        ),
        rbind(
            diag(x = 1, N, L),
            matrix(0, M, L),
            diag(x = 1, L, L)            
        ),
        rbind(
            diag(x = -1, N, L),
            matrix(0, M, L),
            diag(x = 1, L, L)
        )
    )
    bvecSlack <- c(rep(0, M * 2), c(b0, b0 * -1))
    
    ##expand abs value constraint matrices
    AmatAbs <- cbind(
        rbind(
            AmatPosNeg,
            matrix(0, 2 * L, K1)
        ),
        rbind(
            matrix(0, 2 * M, K2),
            AmatPosNegDelta
        )
    )

    ##Map the problem as stated in docs to [b, |b|, |b - b0|] to reduce dimensionality
    MAP <- cbind(
        rbind(
            diag(x = .5, N, M),
            diag(x = .5, M, M),
            matrix(0, L, M)
        ),
        rbind(
            diag(x = -.5, N, M),
            diag(x = .5, M, M),
            matrix(0, L, M)
        ),
        rbind(
            diag(x = .5, N, L),
            matrix(0, M, L),
            diag(x = .5, L, L)
        ),
        rbind(
            diag(x = -.5, N, L),
            matrix(0, M, L),
            diag(x = .5, L, L)
        )
    )

    AmatAbs <- MAP %*% AmatAbs

    AMAT <- cbind(
        Amat,
        AmatSlack,
        AmatAbs
    )

    BVEC <- c(
        bvec,
        bvecSlack,
        bvecPosNeg,
        bvecPosNegDelta
    )

    if(normalize){
        normCons <- normalizeConstraints(AMAT, BVEC)
        AMAT <- normCons$Amat
        BVEC <- normCons$bvec
    }
    
    NVAR <- N + M + L
    DMAT <- matrix(0, NVAR, NVAR)
    diag(DMAT) <- tol
    DMAT[1:N, 1:N] <- Dmat

    if(factorized && (M + L) > 0){
        slackIdx <- (N+1):NVAR
        diag(DMAT)[slackIdx] <- 1 / sqrt(diag(DMAT)[slackIdx])
    }

    ##dvec expansions
    DVECPOSNEGDELTA <- DVECPOSNEG <- matrix(0, 2 * M + 2 * L, 1)
    if(!is.null(dvecPosNeg)){
        DVECPOSNEG[1:(2*M), 1] <- dvecPosNeg
    }
    
    if(!is.null(dvecPosNegDelta)){
        DVECPOSNEG[(2 * M + 1):(2 * M + 2 * L), 1] <- dvecPosNegDelta
    }

    dvec <- c(dvec, rep(0, M + L))
    
    DVEC <- dvec + MAP %*% DVECPOSNEG + MAP %*% DVECPOSNEGDELTA

    if(compact){
        comp <- convertToCompact(AMAT)    
        res <- list(
            Dmat = DMAT,
            dvec = DVEC,
            Amat = comp$Amat,
            Aind = comp$Aind,
            bvec = BVEC,
            meq = meq,
            factorized = factorized
        )
    }else{
        res <- list(
            Dmat = DMAT,
            dvec = DVEC,
            Amat = AMAT,
            bvec = BVEC,
            meq = meq,
            factorized = factorized
        )
    }
    
    res
}
