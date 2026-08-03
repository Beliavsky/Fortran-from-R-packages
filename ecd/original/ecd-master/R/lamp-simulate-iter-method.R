#' Simulate lambda process from stable distribution iteratively
#' 
#' Simulate lambda process from stable distribution iteratively 
#' until target length of result is reached. It uses multi-core capability to 
#' run lamp.simulate1 in parallel. 
#' If file slot is specified, simulation result will be persisted to it periodically.
#' A plot interface is provided to monitor the progress.
#' A CPU temperature interface is provided to control CPU from overheating.
#'
#' @param object an object of lamp class
#' @param use.mc  numeric, number of cores for parallel simulations. Default is 4.
#' @param sim.length  numeric, number of Z to simulate. Default is 1000.
#' @param reset.cache  logical, to reset simulation cache or not prior the run. Default is FALSE.
#' @param keep.tau  numeric, 0 to clean up, 1 to return unused tau, 2 to return all tau. Default is 1.
#' @param drop  numeric, number of tau to discard at the end per iteration. Default is 10.
#' @param plot.util  function, interface to plot simulation results. Default is lamp.plot_sim4.
#' @param cpu.temperature  numeric, temperature above which is overhead. Default is 68.
#' @param cpu.temperature.util  function, interface to get CPU temperature. Default is NULL.
#'
#' @return an object of lamp class with Z, B, N populated
#'
#' @keywords simulation
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#'
### <======================================================================>
"lamp.simulate_iter" <- function(object, use.mc = 4,
                                 sim.length = 1000,
                                 reset.cache = FALSE,
                                 drop=10, keep.tau=1,
                                 plot.util = lamp.plot_sim6,
                                 cpu.temperature = 68,
                                 cpu.temperature.util = NULL) 
{
    
    lp <- object # keep object original and small, add data to lp
    if (reset.cache) {
        lp@Z_i <- lp@tau_i <- 0
        lp@Z <- lp@B <- lp@N <- lp@tau <- numeric(0)
    }

    lambda <- lp@lambda
    T.inf <- lp@T.inf 

    cpu <- function(s, loop=TRUE) {
        if (!is.null(cpu.temperature.util)) {
            while (cpu.temperature.util() > cpu.temperature) { 
                print(paste("Temperature is high:", cpu.temperature.util()))
                Sys.sleep(s)
                if (!loop) return(1)
            }
            return(2)
        }
        return(0)
    }
    
    tm <- Sys.time()
    while (lp@Z_i <= sim.length) {
        print(paste("simulate Z_i=", lp@Z_i, "T.inf=", lp@T.inf, "time=", tm))
        Sys.sleep(1)
        
        f <- function(i) lamp.simulate1(object, drop=drop, keep.tau=keep.tau)
        yy <- parallel::mclapply(seq(1,use.mc), f)
        for (y in yy) {
            lp@Z <- c(lp@Z, y@Z)
            lp@N <- c(lp@N, y@N)
            lp@B <- c(lp@B, y@B)
            lp@Z_i <- length(lp@Z)
            lp@tm <- Sys.time()
        }
        Sys.sleep(1)
        cpu(2, loop=FALSE)
        
        # ---------------------------------------------------
        if (lp@Z_i > 400 & !is.null(plot.util)) plot.util(lp)
        
        if (length(lp@file) > 0) {
            save(lp, file=lp@file)
            sz <- file.info(lp@file)$size/1024.0
            sz1 <- if (sz>=1024) sprintf("%.1fMB", sz/1024) else sprintf("%.1fkB", sz)
            print(paste("   -> data saved, size", sz1, lp@file))
        }
        
        tm <- Sys.time()
        # pause if temperatture gets too high
        cpu(3, loop=TRUE)
        
    }
    return(lp)
}
### <---------------------------------------------------------------------->
