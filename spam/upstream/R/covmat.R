# HEADER ####################################################
# This is file spam/R/covmat.R.                             #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################


# construct various precision matrices



covmat <- function(h, theta, ... , type="sph") {
  avtype <- c("exponential", "spherical", "nugget",
              "wu1","wu2","wu3","wendland1","wendland2",
              "matern12", "matern32", "matern52", "matern")
  method <- pmatch(tolower(type), avtype)
  if (is.na(method))
     stop("Covariance function not implemented yet. Please ask for.")
  switch(method,
         return(cov.exp(h, theta, ...)),
         return(cov.sph(h, theta, ...)),
         return(cov.nug(h, theta, ...)),
         return(cov.wu1(h, theta, ...)),
         return(cov.wu2(h, theta, ...)),
         return(cov.wu3(h, theta, ...)),
         return(cov.wend1(h, theta, ...)),
         return(cov.wend2(h, theta, ...)),
         return(cov.mat12(h, theta, ...)),
         return(cov.mat32(h, theta, ...)),
         return(cov.mat52(h, theta, ...)),
         return(cov.mat(h, theta, ...)))
}

.par.check.cov <- function(theta,nr=2){
  if (any(theta<0)) {
    warning("Parameters coerced to positive values")
    theta <- abs(theta)
  }
  nt <- length(theta)

  if (nt < nr)
    return( c( theta, rep(1, nr-nt), 0))
  return( c( theta, 0)[1:(nr+1)])
}



cov.sph <- function(h, theta, ..., eps=getOption("spam.eps")) {
# approach: we separate   the cases of h spam versus not.
# we distinguish 'no-nugget'-type and correlation matrix.

  theta <- .par.check.cov(theta)    # checking parameters, range, sill, [nugget]
  if (theta[1] < eps) return( cov.nug(h, theta[2] + theta[3]))
         # range so small, is considered as nugget

  if (is.spam(h)) {
      tmp <- h@entries/theta[1]     # (range==1) is unlikely, we do not test
      if (theta[3] > eps) {         # case of a nugget
          ypos <- tmp <= eps        # nugget cases
          npos <- tmp >= 1          # zero cases
          ypos2 <- !(ypos|npos)     # cov cases

          if (any(ypos)) {
              h@entries[ypos] <- theta[2] + theta[3] }
      }else{
          npos <- tmp >= 1          # zero cases
          ypos2 <- !npos            # cov cases
      }                             # end nugget remaining part remains the same

      if (any(ypos2)) {
          ttmp <- tmp[ypos2]
          h@entries[ypos2] <-  if ( abs(theta[2]-1) > eps) { theta[2]*(1 - 1.5 * ttmp + 0.5 * (ttmp*ttmp)*ttmp)
                               }   else { (1 - 1.5 * ttmp + 0.5 * (ttmp*ttmp)*ttmp)}
      }
      if (any(npos)) {
          h@entries[npos] <- 0 }

  } else {
    tmp <- c(h)/theta[1]
      if (theta[3] > eps) {         # case of a nugget
          ypos <- tmp <= eps        # nugget cases
          npos <- tmp >= 1          # zero cases
          ypos2 <- !(ypos|npos)     # cov cases

          if (any(ypos)) {
              h[ypos] <- theta[2] + theta[3] }
      }else{
          npos <- tmp >= 1          # zero cases
          ypos2 <- !npos            # cov cases
      }                             # end nugget remaining part remains the same

      if (any(ypos2)) {
          ttmp <- tmp[ypos2]
          h[ypos2] <-  if ( abs(theta[2]-1) > eps)  theta[2]*(1 - 1.5 * ttmp + 0.5 * (ttmp*ttmp)*ttmp)   else (1 - 1.5 * ttmp + 0.5 * (ttmp*ttmp)*ttmp)
       }
      if (any(npos)) {
          h[npos] <- 0 }

  }
  return(h)
}

cor.sph <- function(h, range, ..., eps=getOption("spam.eps")) {
# approach: we separate   the cases of h spam versus not.
# we distinguish 'no-nugget'-type and correlation matrix.

  theta <- .par.check.cov(range,1)  # checking parameters, range
  if (theta[1] < eps) return( cov.nug(h, 1))  # range so small, is considered as nugget
  if (is.spam(h)) {
      tmp <- h@entries/theta[1]     # (range==1) is unlikely, we do not test
      npos <- tmp >= 1          # zero cases
      ypos2 <- !npos            # cor cases

      if (any(ypos2)) {
          ttmp <- tmp[ypos2]
          h@entries[ypos2] <- (1 - 1.5 * ttmp + 0.5 * (ttmp*ttmp)*ttmp)
      }
      if (any(npos)) {
          h@entries[npos] <- 0 }

  } else {
      tmp <- c(h)/theta[1]
      npos <- tmp >= 1          # zero cases
      ypos2 <- !npos            # cor cases

      if (any(ypos2)) {
          ttmp <- tmp[ypos2]
          h[ypos2] <-  (1 - 1.5 * ttmp + 0.5 * (ttmp*ttmp)*ttmp)
      }
      if (any(npos)) {
          h[npos] <- 0 }
  }
  return(h)
}

cov.wend1 <- function(h, theta,  ... , eps=getOption("spam.eps")) {
  # is \phi_{3,1} in the 98 paper and \psi_{3,1} in the 95 paper, the latter corresponds also
  # to \phi_{\mu,\kappa} in Bevilacqua.
  # For validity it would only be necessary that \mu >= (d+1)/2 + \kappa. In d=2 we would require here
  # \mu >= 2.5. Most theorem in Bevilacqua require \mu > (d+1)/2 + \kappa +d/2, some  \mu > (d+1)/2 + \kappa + 3
  # here this would mean \mu > 3.5 and > 5.5 !
  theta <- .par.check.cov(theta)
  if (theta[1] < eps) return( cov.nug(h, theta[2] + theta[3]))
  # range so small, is considered as nugget

  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1

    if (any(ypos)) {
      h@entries[ypos] <- theta[2] + theta[3] }
    if (any(ypos2)) {
      h@entries[ypos2] <- theta[2]  * ((1 - tmp[ypos2])^4*(4*tmp[ypos2]+1)) }
    if (any(npos)) {
      h@entries[npos] <- 0 }

  } else {
    tmp <- c(h)/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1

    if (any(ypos)) {
      h[ypos] <- theta[2] + theta[3] }
    if (any(ypos2)) {
      h[ypos2] <- theta[2]  * ((1 - tmp[ypos2])^4*(4*tmp[ypos2]+1)) }
    if (any(npos)) {
      h[npos] <- 0 }
  }

  return(h)
}


cov.wend2 <- function(h, theta,  ..., eps=getOption("spam.eps")) {
  # is \phi_{3,2} in the 98 paper and \psi_{4,2} in the 95 paper
  # See comment for cov.wend1. Here smoothness is increased, k=2 twice mean squared differentiable.
  # Simple add "1" to the values above.
  theta <- .par.check.cov(theta)
  if (theta[1] < eps) return( cov.nug(h, theta[2] + theta[3]))
     # range so small, is considered as nugget


  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1

    if (any(ypos)) {
      h@entries[ypos] <- theta[2] + theta[3]}
    if (any(ypos2)) {
      h@entries[ypos2] <- theta[2] * ((1 - tmp[ypos2])^6*(35*tmp[ypos2]^2+18*tmp[ypos2]+3))/3 }
    if (any(npos)) {
      h@entries[npos] <- 0 }

  } else {
    tmp <- c(h)/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1

    if (any(ypos)) {
      h[ypos] <- theta[2] + theta[3] }
    if (any(ypos2)) {
      h[ypos2] <- theta[2] * ((1 - tmp[ypos2])^6*(35*tmp[ypos2]^2+18*tmp[ypos2]+3))/3 }
    if (any(npos)) {
      h[npos] <- 0 }
  }

  return(h)
}


cov.wu1 <- function(h, theta, ... ,  eps=getOption("spam.eps")) {
  theta <- .par.check.cov(theta)
  if (theta[1] < eps) return( cov.nug(h, theta[2] + theta[3]))
  # range so small, is considered as nugget

  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1

    if (any(ypos)) {
      h@entries[ypos] <- theta[2] + theta[3] }
    if (any(ypos2)) {
      h@entries[ypos2] <- theta[2] * ((1 - tmp[ypos2])^3*(1+3*tmp[ypos2]+tmp[ypos2]^2)) }
    if (any(npos)) {
      h@entries[npos] <- 0 }

  } else {
    tmp <- c(h)/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1

    if (any(ypos)) {
      h[ypos] <- theta[2] + theta[3] }
    if (any(ypos2)) {
      h[ypos2] <- theta[2] * ((1 - tmp[ypos2])^3*(1+3*tmp[ypos2]+tmp[ypos2]^2)) }
    if (any(npos)) {
      h[npos] <- 0 }
  }

  return( h)
}


cov.wu2 <- function(h, theta,  ... , eps=getOption("spam.eps")) {
  theta <- .par.check.cov(theta)
  if (theta[1] < eps) return( cov.nug(h, theta[2] + theta[3]))
  # range so small, is considered as nugget

  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1

    if (any(ypos)) {
      h@entries[ypos] <- theta[2] + theta[3] }
    if (any(ypos2)) {
      h@entries[ypos2] <- theta[2] * ((1 - tmp[ypos2])^4*(4+16*tmp[ypos2]+12*tmp[ypos2]^2+3*tmp[ypos2]^3))/4 }
    if (any(npos)) {
      h@entries[npos] <- 0 }

  } else {
    tmp <- c(h)/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1

    if (any(ypos)) {
      h[ypos] <- theta[2] + theta[3]}
    if (any(ypos2)) {
      h[ypos2] <- theta[2] * ((1 - tmp[ypos2])^4*(4+16*tmp[ypos2]+12*tmp[ypos2]^2+3*tmp[ypos2]^3))/4 }
    if (any(npos)) {
      h[npos] <-  0 }
  }

  return( h)
}


cov.wu3 <- function(h, theta,  ..., eps=getOption("spam.eps")) {

  theta <- .par.check.cov(theta)
  if (theta[1] < eps) return( cov.nug(h, theta[2] + theta[3]))
  # range so small, is considered as nugget
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1

    if (any(ypos)) {
      h@entries[ypos] <- theta[2] + theta[3] }
    if (any(ypos2)) {
      h@entries[ypos2] <- theta[2] * ((1 - tmp[ypos2])^6*(1+6*tmp[ypos2]+41/3*tmp[ypos2]^2+12*tmp[ypos2]^3+
                                                            5*tmp[ypos2]^4+5/6*tmp[ypos2]^5)) }
    if(any(npos)) {
      h@entries[npos] <- 0 }

  } else {
    tmp <- c(h)/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1

    if (any(ypos)) {
      h[ypos] <-  theta[2] + theta[3] }
    if (any(ypos2)) {
      h[ypos2] <-  theta[2] * ((1 - tmp[ypos2])^6*(1+6*tmp[ypos2]+41/3*tmp[ypos2]^2+12*tmp[ypos2]^3+
                                                     5*tmp[ypos2]^4+5/6*tmp[ypos2]^5)) }
    if (any(npos)) {
      h[npos] <- 0 }
  }

  return( h)
}


cov.mat <- function(h, theta,  ..., eps=getOption("spam.eps")) {

  theta <- .par.check.cov(theta, 3)
  if (theta[1] < eps) return( cov.nug(h, theta[2] + theta[3]))
  # range so small, is considered as nugget
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    ypos <- tmp < eps
    npos <- !ypos

    if (any(ypos)) {
      h@entries[ypos] <- theta[2] + theta[4] }
    if (any(npos)) {
      h@entries[npos] <- theta[2] * (((2^(-(theta[3] - 1)))/gamma(theta[3])) *
                                       (tmp[npos]^theta[3]) *
                                       besselK(tmp[npos], nu=theta[3])) }
  } else {
    tmp <- c(h)/theta[1]
    ypos <- tmp < eps
    npos <- !ypos

    if (any(ypos)) {
      h[ypos] <- theta[2] + theta[4] }
    if (any(npos)) {
      h[npos] <- theta[2] * (((2^(-(theta[3] - 1)))/gamma(theta[3])) *
                               (tmp[npos]^theta[3]) *
                               besselK(tmp[npos], nu=theta[3])) }
  }

  return(h)
}


cov.finnmat <- function (h, theta, ..., eps = getOption("spam.eps"))
{
  theta <- .par.check.cov(theta, 3)
  if (theta[1] < eps) return( cov.nug(h, theta[2] + theta[3]))
  if (is.spam(h)) {
    tmp <- sqrt(8 * theta[3])/theta[1] * h@entries
    ypos <- tmp < eps
    npos <- !ypos

    if (any(ypos)) {
      h@entries[ypos] <- theta[2] + theta[4]
    }
    if (any(npos)) {
      h@entries[npos] <- theta[2] * (((2^(-(theta[3] - 1)))/gamma(theta[3])) *
                                       (tmp[npos]^theta[3]) *
                                       besselK(tmp[npos], nu = theta[3]))
    }
  } else {
    tmp <- sqrt(8 * theta[3])/theta[1] * c(h)
    ypos <- tmp < eps
    npos <- !ypos
    if (any(ypos)) {
      h[ypos] <- theta[2] + theta[4]
    }
    if (any(npos)) {
      h[npos] <- theta[2] * (((2^(-(theta[3] - 1)))/gamma(theta[3])) *
                               (tmp[npos]^theta[3]) *
                               besselK(tmp[npos], nu = theta[3])) }
  }

  return(h)
}



cov.exp <- function(h, theta, ..., eps=getOption("spam.eps")) {

  theta <- .par.check.cov(theta,2)
  if (theta[1] < eps) return( cov.nug(h, theta[2] + theta[3]))
  # range so small, is considered as nugget
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    ypos <- tmp < eps
    npos <- !ypos

    if (any(ypos)) {
      h@entries[ypos] <- theta[2] + theta[3] }
    if (any(npos)) {
      h@entries[npos] <- theta[2] * exp( -tmp[npos]) }

  } else {
    tmp <- c(h)/theta[1]
    ypos <- tmp < eps
    npos <- !ypos

    if (any(ypos)) {
      h[ypos] <- theta[2] + theta[3] }
    if (any(npos)) {
      h[npos] <- theta[2] * exp( -tmp[npos])   }
  }

  return(h)
}

cov.mat32 <- function(h, theta, ..., eps=getOption("spam.eps")) {

  theta <- .par.check.cov(theta,2)
  if (theta[1] < eps) return( cov.nug(h, theta[2] + theta[3]))
  # range so small, is considered as nugget
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    ypos <- tmp < eps
    npos <- !ypos

    if (any(ypos)) {
      h@entries[ypos] <- theta[2] + theta[3] }
    if (any(npos)) {
      h@entries[npos] <- theta[2] * exp( -tmp[npos]) * (1+tmp[npos]) }

  } else {
    tmp <- c(h)/theta[1]
    ypos <- tmp < eps
    npos <- !ypos

    if (any(ypos)) {
      h[ypos] <- theta[2] + theta[3] }
    if (any(npos)) {
      h[npos] <- theta[2] * exp( -tmp[npos]) * (1+tmp[npos])  }
  }

  return(h)
}

cov.mat52 <- function(h, theta, ..., eps=getOption("spam.eps")) {

  theta <- .par.check.cov(theta,2)
  if (theta[1] < eps) return( cov.nug(h, theta[2] + theta[3]))
  # range so small, is considered as nugget
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    ypos <- tmp < eps
    npos <- !ypos

    if (any(ypos)) {
      h@entries[ypos] <- theta[2] + theta[3] }
    if (any(npos)) {
      h@entries[npos] <- theta[2] * exp( -tmp[npos]) * (1+tmp[npos]+(tmp[npos]^2)/3) }

  } else {
    tmp <- c(h)/theta[1]
    ypos <- tmp < eps
    npos <- !ypos

    if (any(ypos)) {
      h[ypos] <- theta[2] + theta[3] }
    if (any(npos)) {
      h[npos] <- theta[2] * exp( -tmp[npos]) * (1+tmp[npos]+(tmp[npos]^2)/3) }
  }

  return(h)
}

cov.mat12 <- cov.exp

cov.nug <- function(h, theta, ..., eps=getOption("spam.eps")) {

  theta <- .par.check.cov(theta,0)
  if (is.spam(h)) {
    ypos <- h@entries < eps
    npos <- !ypos

    if (any(ypos)) {
      h@entries[ypos] <- theta[1] }
    if (any(npos)) {
      h@entries[npos] <- 0 }

  } else {
    ypos <- c(h) < eps
    npos <- !ypos

    if (any(ypos)) {
      h[ypos] <- theta[1] }
    if (any(npos)) {
      h[npos] <- 0 }
  }

  return(h)
}

