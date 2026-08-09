##*## SSlogisJN - logistic model for nonlinear regression
## Revision by John Nash 2022-8-30 to avoid plinear call.

SSlogisJN <- # selfStart(~ Asym/(1 + exp((xmid - input)/scal)),
    selfStart(
        function(input, Asym, xmid, scal)
        {
              .expr1 <- xmid - input
              .expr3 <- exp(.e2 <- .expr1/scal)
              .expr4 <- 1 + .expr3
              .value <- Asym/.expr4
              .actualArgs <- as.list(match.call()[c("Asym", "xmid", "scal")])
              if(all(vapply(.actualArgs, is.name, NA)))
              {
		  .expr10 <- .expr4^2
                  .grad <- array(0, c(length(.value), 3L), list(NULL, c("Asym", "xmid", "scal")))
                  .grad[, "Asym"] <- 1/.expr4
		  .grad[, "xmid"] <- - (xm <- Asym * .expr3/scal/.expr10)
		  .grad[, "scal"] <- xm * .e2
                  dimnames(.grad) <- list(NULL, .actualArgs)
                  attr(.value, "gradient") <- .grad
              }
              .value
        },
		    initial = function(mCall, data, LHS, ...) {
		          xy <- sortedXyData(mCall[["input"]], LHS, data)
##		          cat("xy:\n"); print(xy)
		          if(nrow(xy) < 4) {
		            stop("too few distinct input values to fit a logistic model")
		          }
		          z <- xy[["y"]]
		          z[which(z<=0)] <- 1e-9 # Just in case
		          Asym<-2*max(z)
		          xy[["z"]] <- log(Asym/z - 1)
		          aux <- coef(lm(z ~ x, xy))
		          scal <- -1/aux[2]
		          xmid <- aux[1]*scal
		          pars<-c(Asym=Asym, xmid=xmid, scal=scal)
		          names(pars)<-c("Asym", "xmid", "scal")
		          pars
		          },
              parameters = c("Asym", "xmid", "scal"))
