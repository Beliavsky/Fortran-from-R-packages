\name{rrcov-utils}
\title{Different utility functions to be used in rrcov and packages depending on rrcov}

\alias{myscreeplot}      % plots screeplots for two covariance matrices 
\alias{mtxconfusion}     % computes a confusion matrix 

%%
%% MATLAB-like functions
%%
\alias{sqrtm}            % returns the square root of a symetric positive definite matrix
%\alias{ones}             % returns a nxp matrix filled with ones 
%\alias{zeros}            % returns a nxp matrix filled with zeros
\alias{vecnorm}          % returns the norm of a vector
%\alias{repmat}           % replicates a matrix

\usage{
    myscreeplot(rcov, ccov)
    mtxconfusion(actual, predicted, prior = NULL, printit=FALSE)
    sqrtm(A)
%    ones(n=1, p=1)
%    zeros(n=1, p=1)
%    repmat(A, n, p)
}

\arguments{
  \item{A}{A matrix.}
%  \item{tol}{Optional numerical tolerance (compared to singular values).}
%  \item{sv}{A vector of non-increasing singular values of a matrix A, if already computed.}
  \item{rcov, ccov}{Covariance matrices.}
%  \item{n}{Number of rows.}
%  \item{p}{Number of columns.}
  \item{actual}{Actual class membership.}
  \item{predicted}{Predicted class membership.}
  \item{prior}{Prior probabilities.}
  \item{printit}{Wheather to print the results.}
}

\description{
Different utility functions to be used in rrcov and packages depending on rrcov:

\itemize{
    \item \code{myscreeplot} Plots screeplots for two covariance matrices 
    \item \code{mtxconfusion} Computes a confusion matrix between actual and predicted class membership

    \item MATLAB-like functions
    \itemize{ 
        \item \code{sqrtm} Returns the square root of a symetric positive definite matrix
        \item \code{ones} Returns a nxp matrix filled with ones 
        \item \code{zeros} Returns a nxp matrix filled with zeros
        \item \code{vecnorm} Returns the norm of a vector
        \item \code{rankMM} Returns the rank of the matrix
        \item \code{repmat} Replicates a matrix
    }
}
}
%\details{}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\keyword{internal} 
