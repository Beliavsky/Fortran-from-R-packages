#' Domain coloring of Laplace kernel of lambda distribution
#' 
#' Domain coloring on the complex plane of Laplace kernel of lambda distribution, 
#' \code{exp(ixt) P(x)}, where P(x) is the PDF of a lambda distribution.
#' This is a visualization utility to get insight how the Laplace transform works
#' for lambda distribution. The behavior on the complex plane is deeply associated
#' with the MGF, the skew Levy distribution, and the SaS distribution.
#'
#' @param t numeric or complex
#' @param rec numeric, define the rectangle of plot in the order of (x1, x2, y1, y2)
#' @param n numeric, number of points per axis. Default is 200. Use 1000 for better resolution.
#' @param lambda numeric. Default is 4, and is the only value allowed.
#' 
#' @return return value of call to \code{grid.arrange()}
#'
#' @keywords Domain-Coloring
#'
#' @author Stephen H. Lihn
#'
#' @export
#' 
#' @importFrom grDevices hsv
#' @importFrom parallel mclapply
#' @importFrom ggplot2 ggplot 
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 qplot
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 coord_flip
#' @importFrom gridExtra grid.arrange
#' 
#' @examples
#' \dontrun{
#'   levy.domain_coloring(0.1, c(-25, 50, -50, 50))
#'   levy.domain_coloring(0.1i, c(-25, 25, -25, 25))
#' }
#'
### <======================================================================>
"levy.domain_coloring" <- function(t, rec, n=200, lambda=4)
{
    if (lambda!=4) stop(paste("lambda is not supported", lambda))
    if (Arg(t) != 0 & Arg(t)/pi != 0.5) stop("t must be real or imaginary")
        
    t <- t+0 # to ensure branch cut is handled correctly
    f <- function(x) exp(t*x)*levy.dlambda(x, lambda)
    
    # TODO need to formalize this for lambda !?
    cut = Re(1/t^2/4)
    if (Im(t)>0) cut = Re(1/t^2/8)*1i+0
    
    x = seq(rec[1], rec[2], length.out=n) 
    y = seq(rec[3], rec[4], length.out=n)
    
    z = as.vector(outer(x+0i, 0+1i*y, '+'))
    f_z = f(z)
    max_f_z = max(Mod(f_z))
    scale = (exp(1)-1)*0.9/max_f_z
    r0 = log(1+Mod(f_z)*scale) # replace the standard approach, log(1+Mod(f(z)))
    
    levels <- 25 # number of levels in the intensity
    rl <- c(max(r0)*0.9) # array of levels
    for (i in 2:levels) rl[i] <- rl[i-1]*0.8
    
    f_cut <- Mod(f(cut))
    rl_cut_lst <- which(rl<=f_cut)
    if (length(rl_cut_lst) > 0) {
        rl_cut <- min(which(rl<=f_cut))
        delta_rl <- f_cut-rl[rl_cut]
        rl <- rl + delta_rl # normalize one level to f(cut)
    }
    
    # modularize r into levels for better visulization
    get_level <- function(r) {
        i <- which(r >= rl)
        if (length(i)>0) rl[min(i)] else r
    }
    mc <- function(x,y) simplify2array(parallel::mclapply(x,y))
    r = mc(r0, get_level)
    
    max_rad = 2*pi # pi or 2*pi
    z =data.frame(x=Re(z),
                  y=Im(z),
                  h=(Arg(f_z)<0)*1+Arg(f_z)/(2*pi),
                  s=(1+sin(max_rad*r))/2,
                  v=(1+cos(max_rad*r))/2)
    
    z = z[is.finite(apply(z,1,sum)),]
    
    # -------------------------------------------------------
    lay <- rbind(c(3,1,1),
                 c(3,1,1),
                 c(4,2,2))

    # https://aschinchon.wordpress.com/2014/10/01/complex-domain-coloring/
    # http://www.jedsoft.org/fun/complex/fztopng.html
    
    plot1 <- ggplot2::ggplot(data=z, ggplot2::aes(x=x, y=y)) + 
        ggplot2::geom_tile(fill=grDevices::hsv(z$h,z$s,z$v)) 

    f_x <- log(Mod(f(x+0i)))
    f_y <- log(Mod(f(0+1i*y)))
    
    min_f_x <- min(f_x[x>0])
    x_cut <- x[f_x==min_f_x & x>0]
    min_f_y <- min(f_y[y<0])
    y_cut <- y[f_y==min_f_y & y<0]
    
    plot2 <- ggplot2::qplot(x, f_x, geom="line") + 
        ggplot2::ylab("log Mod f(x)")
    plot3 <- ggplot2::qplot(y, f_y, geom="line") + 
        ggplot2::ylab("log Mod f(yi)") + ggplot2::coord_flip()
    
    if (x_cut < max(x)) plot2 <- plot2 + ggplot2::geom_vline(xintercept=x_cut, linetype="dotted")
    if (y_cut > min(y)) plot3 <- plot3 + ggplot2::geom_vline(xintercept=y_cut, linetype="dotted")
    
    plot4 <- grid::rectGrob(gp=grid::gpar(fill="white", col="white"))
    
    gridExtra::grid.arrange(plot1, plot2, plot3, plot4, layout_matrix = lay)
}
### <---------------------------------------------------------------------->
