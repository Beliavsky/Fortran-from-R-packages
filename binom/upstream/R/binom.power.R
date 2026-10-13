binom.power <- function(p.alt,
                        n = 100,
                        p = 0.5,
                        alpha = 0.05,
                        phi = 1,
                        alternative = c("two.sided", "greater", "less"),
                        method = c("cloglog", "logit", "probit", "asymp", "lrt", "exact")) {
  args <- cbind(p.alt, n, p, alpha, phi)
  p.alt <- args[, "p.alt"]
  n <- args[, "n"]
  p <- args[, "p"]
  alpha <- args[, "alpha"]
  phi <- args[, "phi"]
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  if(method %in% c("cloglog", "logit", "probit")) {
    linkfun <- eval(substitute(binomial(x)$linkfun, list(x = method)))
    varfun <- get(paste("var", method, sep = "."))
  } else if(method == "asymp") {
    linkfun <- function(mu) mu
    varfun <- get(sprintf("var.%s", method))
  } else {
    linkfun <- varfun <- NULL
  }
  if(alternative == "two.sided") {
    z <- qnorm(1 - alpha/2)
  } else {
    z <- qnorm(1 - alpha)
  }
  if(!is.null(linkfun)) {
    cloglog <- method == "cloglog"
    gamma0 <- linkfun(if(cloglog) 1 - p else p)
    gamma1 <- linkfun(if(cloglog) 1 - p.alt else p.alt)
    sd0 <- sqrt(phi * varfun(p, n))
    sd1 <- sqrt(phi * varfun(p.alt, n))
    pz0 <- pnorm((gamma1 - gamma0 - z * sd0)/sd1) # decrease in p
    pz1 <- pnorm((gamma0 - gamma1 - z * sd0)/sd1) # increase in p
    power <- rep(alpha, length.out = nrow(args))
    p.diff <- p != p.alt
    power[p.diff] <- if (alternative == "less") {
      (if(cloglog) pz0 else pz1)[p.diff]
    } else if (alternative == "greater") {
      (if(cloglog) pz1 else pz0)[p.diff]
    } else {
      (pz0 + pz1)[p.diff]
    }
  } else if(method %in% c("lrt", "exact")) {
    alpha <- if(alternative == "two.sided") alpha else 2 * alpha
    pci <- binom.confint(n * p, n, 1 - alpha, method)[c("lower", "upper")]
    # P(X > n * pci[2] | p = p.alt)
    # pbinom(n * pci[2], n, p.alt, FALSE)
    k <- length(p.alt)
    power <- numeric(k)
    upper <- lower <- 0
    for(i in seq(k)) {
      x <- n * (1 - p.alt[i])
      if(alternative != "less")
        lower <- pbeta(pci$upper, n + 1 - x, x, lower.tail = FALSE)
      if(alternative != "greater")
        upper <- pbeta(pci$lower, n - x, x + 1)
      power[i] <- upper + lower
    }
  }
  power
}

tkbinom.power2 <- function() {
  stopifnot(requireNamespace("tcltk"))
  stopifnot(requireNamespace("lattice"))
  lattice::trellis.par.set(theme = lattice::col.whitebg())
  local({
    tk.p <- tcltk::tclVar(0.5)
    tk.alpha <- tcltk::tclVar(0.05)
    tk.n <- tcltk::tclVar(25)
    tk.alternative <- tcltk::tclVar("two.sided")
    methods <- c("cloglog", "logit", "probit", "asymp", "exact", "lrt")
    tk.method <- lapply(methods, tcltk::tclVar)
    col <- c("#800000", "#000080", "#008000", "#e00e00", "#808000", "#0000ff")
    names(col) <- names(tk.method) <- methods
    p.sav <- 0.5 # in case replot.maybe is called too early
    replot <- function(...) {
      p <- as.numeric(tcltk::tclObj(tk.p))
      n <- as.numeric(tcltk::tclObj(tk.n))
      alpha <- as.numeric(tcltk::tclObj(tk.alpha))
      alternative <- as.character(tcltk::tclObj(tk.alternative))
      p.alt <- seq(0, 1, 0.01)[-c(1, 101)]
      p.sav <<- p <- as.numeric(tcltk::tclObj(tk.p))
      first <- TRUE
      delta <- p.alt - p
      xlim <- range(delta)
      plot(xlim, c(0, 1), type = "n", axes = FALSE,
           xlab = "Prob(Success)", ylab = "Power", lwd = 3)
      axis(side = 1, at = seq(0, 1, 0.1) - p, labels = seq(0, 1, 0.1))
      axis(side = 2)
      axis(side = 3, at = pretty(delta), labels = round(pretty(delta), 2))
      box()
      show <- NULL
      for(i in methods) {
        plot.method <- as.logical(tcltk::tclObj(tk.method[[i]]))
        if(!is.na(plot.method) && plot.method) {
          show <- c(show, i)
          y <- binom.power(p.alt, n, p, alpha, 1, alternative, i)
          lines(delta, y, col = col[i], lwd = 3)
        }
      }
      abline(h = as.numeric(tcltk::tclObj(tk.alpha)), v = 0)
      if(!is.null(show)) {
        legend(x = min(delta), y = 1, legend = show, col = col[show], lty = 1, lwd = 3)
      }
    }
    replot.maybe <- function(...) {
      if(as.numeric(tcltk::tclObj(tk.p)) != p.sav) replot()
    }
    base <- tcltk::tktoplevel()
    tcltk::tkwm.title(base, "Binomial Power Curves")
    spec.frm <- tcltk::tkframe(base, borderwidth = 2)
    left.frm <- tcltk::tkframe(spec.frm)
    right.frm <- tcltk::tkframe(spec.frm)
    bottom.frm <- tcltk::tkframe(spec.frm)
    ## Two left frames:
    frame1 <- tcltk::tkframe(left.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame1, text = "Alternative"))
    for(i in c("two.sided", "greater", "less")) {
      tmp <- tcltk::tkradiobutton(frame1, command = replot, text = i,
                           value = i, variable = tk.alternative)
      tcltk::tkpack(tmp, anchor = "w")
    }
    frame2 <- tcltk::tkframe(left.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame2, text = "Method"))
    for(i in c("cloglog", "logit", "probit", "asymp", "exact", "lrt") ) {
      tmp <- tcltk::tkcheckbutton(frame2, text = i, command = replot,
                           variable = tk.method[[i]])
      tcltk::tkpack(tmp, anchor = "w")
    }
    ## Two right frames:
    frame3 <- tcltk::tkframe(right.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame3, text = "Sample size"))
    for(i in c(25, 50, 75, 100, 125)) {
      tmp <- tcltk::tkradiobutton(frame3, command = replot, text = i,
                           value = i, variable = tk.n)
      tcltk::tkpack(tmp, anchor = "w")
    }
    frame4 <- tcltk::tkframe(right.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame4, text = "Prob(Success)"))
    tcltk::tkpack(tcltk::tkscale(frame4, command = replot.maybe,
                   from = 0.001, to = 0.999,
                   showvalue = TRUE, variable = tk.p,
                   resolution = 0.001, orient = "horiz"))
    ## build layout
    tcltk::tkpack(frame1, frame2, fill = "x")
    tcltk::tkpack(frame3, frame4, fill = "x")
    tcltk::tkpack(left.frm, right.frm, side = "left", anchor = "n")
    ## bottom frame:
    q.but <- tcltk::tkbutton(base, text = "Quit", command = function() tcltk::tkdestroy(base))
    tcltk::tkpack(spec.frm, q.but)
    replot()
  })
  invisible()
}

tkbinom.power <- function() {
  stopifnot(requireNamespace("tcltk"))
  stopifnot(requireNamespace("lattice"))
  lattice::trellis.par.set(theme = lattice::col.whitebg())
  local({
    tk.p <- tcltk::tclVar(0.5)
    tk.alpha <- tcltk::tclVar(0.05)
    tk.n <- tcltk::tclVar(25)
    tk.alternative <- tcltk::tclVar("two.sided")
    methods <- c("cloglog", "logit", "probit", "asymp", "exact", "lrt")
    tk.method <- lapply(methods, tcltk::tclVar)
    col <- lattice::trellis.par.get("superpose.line")$col[1:6]
    names(col) <- names(tk.method) <- methods
    p.sav <- 0.5 # in case replot.maybe is called too early
    replot <- function(...) {
      p <- as.numeric(tcltk::tclObj(tk.p))
      n <- as.numeric(tcltk::tclObj(tk.n))
      alpha <- as.numeric(tcltk::tclObj(tk.alpha))
      alternative <- as.character(tcltk::tclObj(tk.alternative))
      p.alt <- seq(0, 1, 0.01)[-c(1, 101)]
      p.sav <<- p <- as.numeric(tcltk::tclObj(tk.p))
      first <- TRUE
      delta <- p.alt - p
      xlim <- range(delta)
      ylim <- c(0, 1)
      show <- x <- y <- key <- NULL
      for(i in methods) {
        plot.method <- as.logical(tcltk::tclObj(tk.method[[i]]))
        if(!is.na(plot.method) && plot.method) {
          y <- c(y, binom.power(p.alt, n, p, alpha, 1, alternative, i))
          show <- c(show, rep(i, length(p.alt)))
          x <- c(x, delta)
        }
      }
      if(!is.null(show)) {
        show <- factor(show)
        key <- list(space = "right", text = list(levels(show)),
                    lines = list(col = col[levels(show)], lty = 1, lwd = 3))
      }
      xy <- lattice::xyplot(ylim ~ xlim, panel = function(x, y, .x, .y, .g, col, ...) {
                     if(!is.null(.y)) {
                       for(g in levels(.g)) {
                         s <- .g == g
                         lattice::llines(.x[s], .y[s], col = col[g], ...)
                       }
                     }
                     lattice::panel.abline(h = as.numeric(tcltk::tclObj(tk.alpha)), v = 0)
                   },
                   .x = x, .y = y, .g = show,
                   key = key, type = "l", lwd = 3, col = col,
                   xlab = "Delta", ylab = "Power")
      print(xy)
    }
    replot.maybe <- function(...) {
      if(as.numeric(tcltk::tclObj(tk.p)) != p.sav) replot()
    }
    base <- tcltk::tktoplevel()
    tcltk::tkwm.title(base, "Binomial Power Curves")
    spec.frm <- tcltk::tkframe(base, borderwidth = 2)
    left.frm <- tcltk::tkframe(spec.frm)
    right.frm <- tcltk::tkframe(spec.frm)
    bottom.frm <- tcltk::tkframe(spec.frm)
    ## Two left frames:
    frame1 <- tcltk::tkframe(left.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame1, text = "Alternative"))
    for(i in c("two.sided", "greater", "less")) {
      tmp <- tcltk::tkradiobutton(frame1, command = replot, text = i,
                           value = i, variable = tk.alternative)
      tcltk::tkpack(tmp, anchor = "w")
    }
    frame2 <- tcltk::tkframe(left.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame2, text = "Method"))
    for(i in c("cloglog", "logit", "probit", "asymp", "exact", "lrt") ) {
      tmp <- tcltk::tkcheckbutton(frame2, text = i, command = replot,
                           variable = tk.method[[i]])
      tcltk::tkpack(tmp, anchor = "w")
    }
    ## Two right frames:
    frame3 <- tcltk::tkframe(right.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame3, text = "Sample size"))
    for(i in c(25, 50, 75, 100, 125)) {
      tmp <- tcltk::tkradiobutton(frame3, command = replot, text = i,
                           value = i, variable = tk.n)
      tcltk::tkpack(tmp, anchor = "w")
    }
    frame4 <- tcltk::tkframe(right.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame4, text = "Prob(Success)"))
    tcltk::tkpack(tcltk::tkscale(frame4, command = replot.maybe,
                   from = 0.001, to = 0.999,
                   showvalue = TRUE, variable = tk.p,
                   resolution = 0.001, orient = "horiz"))
    ## build layout
    tcltk::tkpack(frame1, frame2, fill = "x")
    tcltk::tkpack(frame3, frame4, fill = "x")
    tcltk::tkpack(left.frm, right.frm, side = "left", anchor = "n")
    ## bottom frame:
    q.but <- tcltk::tkbutton(base, text = "Quit", command = function() tcltk::tkdestroy(base))
    tcltk::tkpack(spec.frm, q.but)
    replot()
  })
  invisible()
}
