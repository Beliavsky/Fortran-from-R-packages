### <======================================================================>
#' Plot ES contribution 
#'
#'These functions plot the contribution of each asset to the overall portfolio expected shortfall.
#'
#'
#' @docType methods
#' @importFrom graphics plot
#' @name plot-ghyp.attribution
#' @rdname  plot-ghyp.attribution
#' @aliases  plot,ghyp.attribution,ANY-method
#'
#' @param x A \code{ghyp.attribution} object.
#' @param metrics either the \code{contribution} or \code{sensitivity} will be plotted.
#' @param column.index which column of the object.
#' @param percentage plot contribution or sensitivity in percent.
#' @param colorset vector of colors for the chart.
#' @param horiz plot horizontally.
#' @param unstacked unstacked plot.
#' @param pie.chart should a pie chart be plotted.
#' @param sub subtitle.
#' @param \dots arguments passed to \code{plot} function.
#'
#' @author Marc Weibel
#' @seealso \code{\link{ESghyp.attribution}}.
#' @keywords attribution
#' @examples 
#' \dontrun{
#'  data(smi.stocks)
#'  
#'  ## Fit a NIG model to Novartis, CS and Nestle log-returns
#'  assets.fit <- fit.NIGmv(smi.stocks[, c("Novartis", "CS", "Nestle")], silent = TRUE)
#'  
#'  ## Define Weights of the Portfolio
#'  weights <- c(0.2, 0.5, 0.3)
#'  
#'  ## Confidence level for Expected Shortfall
#'  es.levels <- c(0.01)
#'  portfolio.attrib <- ESghyp.attribution(alpha=es.levels, object=assets.fit, weights=weights)
#'  
#'  ## Plot Risk Contribution for each Asset
#'  plot(portfolio.attrib, metrics='contribution')  
#' }
#' @export
"plot.ghyp.attrib" <- function(x, metrics=c('contribution', 'sensitivity'), 
                               column.index=NULL, percentage=FALSE, colorset=NULL,
                               horiz=FALSE, unstacked=TRUE,
                               pie.chart=FALSE, sub=NULL, ...)
{		
  
  metrics = match.arg(metrics)	
  
  if(metrics!='contribution' && percentage==TRUE)
    stop('Percentage can only be chosen with contribution ! \n')
  if(metrics!='contribution' && pie.chart==TRUE)
    stop('Pie Chart can only be chosen with contribution and percentage set as TRUE ! \n')	
  
  object <- eval(parse(text=paste('x@', metrics, sep="")))	
  colNames <- colnames(object)
  if(!is.null(column.index)) 
  {			
    object <- as.matrix(object[, column.index])
    colnames(object) <- colNames[column.index]		
    if(is.null(sub)) sub <- paste('Probability = ', 
                                  colnames(object), 
                                  sep="")
  }
  n.row <- NROW(object)
  n.col <- NCOL(object)	
  
  ## Stacked.Plot do not make sense for Sensitivity
  ## as it's not additive to overall Portfolio
  if(n.col>1 && metrics=='sensitivity') 
    stop('Only one-dimensional objects for Sensitivity Chart ! \n')
  
  if(n.col>1 && pie.chart==TRUE) 
    stop('Only one-dimensional objects for Pie Chart ! \n')
  
  
  ## If pie chart was chosen, set percentage as TRUE
  if(metrics=='contribution' && pie.chart==TRUE && percentage==FALSE) 
  {
    cat('percentage has been set to TRUE for pie chart ! \n')
    percentage=TRUE
  }
  
  if(metrics=='contribution') 
  {
    metrics <- 'Contribution'
    if(percentage==TRUE) metrics <- 'Contribution (in %)'
  } else {
    metrics <- 'Sensitivity'
  }
  
  ## Contribution in Percent
  if(percentage==TRUE) object <- t(t(object)/colSums(object)) * 100	
  
  ## Produce a bar plot or a pie chart
  if(pie.chart==FALSE)
  {
    plot.StackedBar(t(object), xlab='Probability', 
                    ylab=metrics, 
                    main=paste('Expected Shortfall ', 
                               metrics, sep=""), horiz=horiz, 
                    colorset, sub=sub, unstacked = unstacked, ...)
    
  } else {		
    if(is.null(colorset)) colorset=.my.pal(n.row,'topo')		
    plot.PieChart(object, 
                  labels=paste(rownames(object)," (",round(object,2),"%)",sep=""), 
                  main='Expected Shortfall Contribution (in %)', 
                  colorset=colorset, 
                  sub=paste('Probability = ', 
                            colnames(object), sep=""),
                  ...)
  }		
}
### <---------------------------------------------------------------------->
setMethod("plot", signature(x = "ghyp.attribution"), plot.ghyp.attrib)
### <---------------------------------------------------------------------->

### <======================================================================>
##------ Color palettes ------------
".my.pal" <- function(n, palette=c('blues','rainbow', 'heat', 'terrain', 'topo', 'cm'))
{	
  palette <- match.arg(palette)
  ch.col = c("rainbow(n, start=.7, end=.1)", "heat.colors(n)",
             "terrain.colors(n)", "topo.colors(n)","cm.colors(n)")
  
  nt <- length(ch.col)        
  colors <- matrix(0,nt, n)
  for (k in 1:nt) {        
    colors[k,] = eval(parse(text=ch.col[k]))      
  }
  rownames(colors) <- c('rainbow', 'heat', 'terrain', 'topo', 'cm')
  return(colors[which(rownames(colors)==palette), ])
}
### <---------------------------------------------------------------------->

seqPalette <- function (n, name = c("Blues", "BuGn", "BuPu", "GnBu", "Greens", 
                                    "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", 
                                    "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")) 
{
  Blues = rgb(c(247, 222, 198, 158, 107, 66, 33, 8, 8), c(251, 
                                                          235, 219, 202, 174, 146, 113, 81, 48), c(255, 247, 239, 
                                                                                                   225, 214, 198, 181, 156, 107), maxColorValue = 255)
  BuGn = rgb(c(247, 229, 204, 153, 102, 65, 35, 0, 0), c(252, 
                                                         245, 236, 216, 194, 174, 139, 109, 68), c(253, 249, 230, 
                                                                                                   201, 164, 118, 69, 44, 27), maxColorValue = 255)
  BuPu = rgb(c(247, 224, 191, 158, 140, 140, 136, 129, 77), 
             c(252, 236, 211, 188, 150, 107, 65, 15, 0), c(253, 244, 
                                                           230, 218, 198, 177, 157, 124, 75), maxColorValue = 255)
  GnBu = rgb(c(247, 224, 204, 168, 123, 78, 43, 8, 8), c(252, 
                                                         243, 235, 221, 204, 179, 140, 104, 64), c(240, 219, 197, 
                                                                                                   181, 196, 211, 190, 172, 129), maxColorValue = 255)
  Greens = rgb(c(247, 229, 199, 161, 116, 65, 35, 0, 0), c(252, 
                                                           245, 233, 217, 196, 171, 139, 109, 68), c(245, 224, 192, 
                                                                                                     155, 118, 93, 69, 44, 27), maxColorValue = 255)
  Greys = rgb(c(255, 240, 217, 189, 150, 115, 82, 37, 0), c(255, 
                                                            240, 217, 189, 150, 115, 82, 37, 0), c(255, 240, 217, 
                                                                                                   189, 150, 115, 82, 37, 0), maxColorValue = 255)
  Oranges = rgb(c(255, 254, 253, 253, 253, 241, 217, 166, 127), 
                c(245, 230, 208, 174, 141, 105, 72, 54, 39), c(235, 206, 
                                                               162, 107, 60, 19, 1, 3, 4), maxColorValue = 255)
  OrRd = rgb(c(255, 254, 253, 253, 252, 239, 215, 179, 127), 
             c(247, 232, 212, 187, 141, 101, 48, 0, 0), c(236, 200, 
                                                          158, 132, 89, 72, 31, 0, 0), maxColorValue = 255)
  PuBu = rgb(c(255, 236, 208, 166, 116, 54, 5, 4, 2), c(247, 
                                                        231, 209, 189, 169, 144, 112, 90, 56), c(251, 242, 230, 
                                                                                                 219, 207, 192, 176, 141, 88), maxColorValue = 255)
  PuBuGn = rgb(c(255, 236, 208, 166, 103, 54, 2, 1, 1), c(247, 
                                                          226, 209, 189, 169, 144, 129, 108, 70), c(251, 240, 230, 
                                                                                                    219, 207, 192, 138, 89, 54), maxColorValue = 255)
  PuOr = rgb(c(127, 179, 224, 253, 254, 247, 216, 178, 128, 
               84, 45), c(59, 88, 130, 184, 224, 247, 218, 171, 115, 
                          39, 0), c(8, 6, 20, 99, 182, 247, 235, 210, 172, 136, 
                                    75), maxColorValue = 255)
  PuRd = rgb(c(247, 231, 212, 201, 223, 231, 206, 152, 103), 
             c(244, 225, 185, 148, 101, 41, 18, 0, 0), c(249, 239, 
                                                         218, 199, 176, 138, 86, 67, 31), maxColorValue = 255)
  Purples = rgb(c(252, 239, 218, 188, 158, 128, 106, 84, 63), 
                c(251, 237, 218, 189, 154, 125, 81, 39, 0), c(253, 245, 
                                                              235, 220, 200, 186, 163, 143, 125), maxColorValue = 255)
  RdPu = rgb(c(255, 253, 252, 250, 247, 221, 174, 122, 73), 
             c(247, 224, 197, 159, 104, 52, 1, 1, 0), c(243, 221, 
                                                        192, 181, 161, 151, 126, 119, 106), maxColorValue = 255)
  Reds = rgb(c(255, 254, 252, 252, 251, 239, 203, 165, 103), 
             c(245, 224, 187, 146, 106, 59, 24, 15, 0), c(240, 210, 
                                                          161, 114, 74, 44, 29, 21, 13), maxColorValue = 255)
  YlGn = rgb(c(255, 247, 217, 173, 120, 65, 35, 0, 0), c(255, 
                                                         252, 240, 221, 198, 171, 132, 104, 69), c(229, 185, 163, 
                                                                                                   142, 121, 93, 67, 55, 41), maxColorValue = 255)
  YlGnBu = rgb(c(255, 237, 199, 127, 65, 29, 34, 37, 8), c(255, 
                                                           248, 233, 205, 182, 145, 94, 52, 29), c(217, 177, 180, 
                                                                                                   187, 196, 192, 168, 148, 88), maxColorValue = 255)
  YlOrBr = rgb(c(255, 255, 254, 254, 254, 236, 204, 153, 102), 
               c(255, 247, 227, 196, 153, 112, 76, 52, 37), c(229, 188, 
                                                              145, 79, 41, 20, 2, 4, 6), maxColorValue = 255)
  YlOrRd = rgb(c(255, 255, 254, 254, 253, 252, 227, 189, 128), 
               c(255, 237, 217, 178, 141, 78, 26, 0, 0), c(204, 160, 
                                                           118, 76, 60, 42, 28, 38, 38), maxColorValue = 255)
  name = match.arg(name)
  orig = eval(parse(text = name))
  rgb = t(col2rgb(orig))
  temp = matrix(NA, ncol = 3, nrow = n)
  x = seq(0, 1, , length(orig))
  xg = seq(0, 1, , n)
  for (k in 1:3) {
    hold = spline(x, rgb[, k], n = n)$y
    hold[hold < 0] = 0
    hold[hold > 255] = 255
    temp[, k] = round(hold)
  }
  palette = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
  palette
}

### <======================================================================>
"plot.StackedBar" <- function (object, colorset = NULL, horiz=FALSE,
                               space = 0.2, cex.axis=0.8, 
                               cex.legend = 0.8, cex.lab = 1, 
                               cex.labels = 0.8, cex.main = 1, 
                               xaxis=TRUE, legend.loc="under",  
                               element.color = "darkgray", unstacked = TRUE, 
                               xlab="Date", ylab="Value", ylim=NULL, 
                               date.format = "%b %y", 
                               major.ticks='auto', minor.ticks=TRUE, 
                               las = 0, xaxis.labels = NULL, ... ) 
{
  ## Data should be organized as columns for each category, 
  ## rows for each period or observation    
  object.columns = NCOL(object)
  object.rows = NROW(object)
  
  posn = barplot(t(object), plot=FALSE, space=space)
  
  if(is.null(colnames(object)))
    legend.loc = NULL
  
  if(is.null(colorset))
    colorset=seqPalette(object.columns, 'Blues')
  
  if(is.null(xlab))
    minmargin = 3
  else
    minmargin = 5
  
  if(unstacked & dim(object)[1] == 1){ # only one row is being passed into 'object', unstack the bars
    if(las > 1) {						# set the bottom border to accomodate labels
      bottommargin = max(c(minmargin, (strwidth(colnames(object),units="in")) / 
                             par("cin")[1])) * cex.lab
      par(mar = c(bottommargin, 4, 4, 2) +.1)
    }
    barplot(object, col = '#9ECAE1', las = las, horiz = horiz, 
            space = space, xlab = "", cex.names = cex.lab, axes = FALSE, ylim=ylim, ...)
    if (horiz==TRUE)
      axis(1, col = element.color, las = las, cex.axis = cex.axis)
    else
      axis(2, col = element.color, las = las, cex.axis = cex.axis)
    box(col = element.color)
  }
  else { # multiple columns being passed into 'object', stack the bars and put a legend underneith
    if(!is.null(legend.loc) ){
      if(legend.loc =="under") {	# put the legend under the chart
        op <- par(no.readonly=TRUE)
        layout(rbind(1,2), heights=c(6,1), widths=1)
        par(mar=c(3,4,4,2)+.1) 	# set the margins of the first panel
      }
    }
    # Brute force solution for plotting negative values in the bar charts:
    positives = object
    for(column in 1:ncol(object)){
      for(row in 1:nrow(object)){ 
        positives[row,column]=max(0,object[row,column])
      }
    }
    
    negatives = object
    for(column in 1:ncol(object)){
      for(row in 1:nrow(object)){ 
        negatives[row,column]=min(0,object[row,column])
      }
    }
    # Set ylim accordingly
    if(is.null(ylim)){
      ymax=max(0,apply(positives,FUN=sum,MARGIN=1))
      ymin=min(0,apply(negatives,FUN=sum,MARGIN=1))
      ylim=c(ymin,ymax)
    }
    if (horiz==TRUE)
    {
      barplot(t(positives), col=colorset, horiz=horiz,
              space=space, 
              axisnames = FALSE, axes = FALSE, 
              xlim=ylim, xlab="", ...)
      barplot(t(negatives), add=TRUE , col=colorset,
              horiz=horiz,
              space=space, las = las, xlab = "", 
              cex.names = cex.lab, 
              axes = FALSE, axisnames = FALSE, xlim=ylim, ...)    
      axis(1, col = element.color, las = las, cex.axis = cex.axis)
    } else {
      barplot(t(positives), col=colorset, horiz=horiz,
              space=space, 
              axisnames = FALSE, axes = FALSE, 
              ylim=ylim, xlab="", ...)
      barplot(t(negatives), add=TRUE , col=colorset,
              horiz=horiz,
              space=space, las = las, xlab = "", 
              cex.names = cex.lab, 
              axes = FALSE, axisnames = FALSE, ylim=ylim, ...)   
      axis(2, col = element.color, las = las, cex.axis = cex.axis)
    }
    
    title(ylab = ylab, cex = cex.lab)
    if (xaxis) {
      label.height = .25 + cex.axis * 
        max(strheight(rownames(object), units="in") / 
              par('cin')[2])
      if(is.null(xaxis.labels))
        xaxis.labels = rownames(object)
      if (horiz==TRUE)
        axis(2, at=posn, labels=xaxis.labels, las=las, lwd=1, 
             mgp=c(3,label.height,0), cex.axis = cex.axis)
      else
        axis(1, at=posn, labels=xaxis.labels, las=las, lwd=1, 
             mgp=c(3,label.height,0), cex.axis = cex.axis)
    }
    box(col = element.color)
    
    if(!is.null(legend.loc)){
      if(legend.loc =="under"){ 	# draw the legend under the chart
        par(mar=c(0,2,0,1)+.1) 	# set the margins of the second panel
        plot.new()
        if(object.columns <4)
          ncol= object.columns
        else
          ncol = 4
        legend("center", legend=colnames(object), cex = cex.legend, 
               fill=colorset, ncol=ncol,
               box.col=element.color, border = element.color)
        par(op)
      } # if legend.loc is null, then do nothing
    }
  }
}
### <---------------------------------------------------------------------->


### <======================================================================>
plot.PieChart <- function (x, labels = names(x), labels.loc = "under", 
                           edges = 200, radius = 0.8, 
                           clockwise = FALSE, init.angle = if (clockwise) 90 else 0, 
                           density = NULL, angle = 45, 
                           colorset = NULL, 
                           border = NULL, 
                           lty = NULL, 
                           main = NULL, 
                           ...) 
{
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  object.columns = NCOL(x)
  object.rows = NROW(x)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) 
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(colorset)) 
    colorset <- if (is.null(density)) 
      colorset=seqPalette(object.columns, 'Blues')
  else par("fg")
  colorset <- rep(colorset, length.out = nx)
  border <- rep(border, length.out = nx)
  lty <- rep(lty, length.out = nx)
  angle <- rep(angle, length.out = nx)
  density <- rep(density, length.out = nx)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = colorset[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      label.name <- 
        lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
      text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
           adj = ifelse(P$x < 0, 1, 0), ...)
    }
  }
  title(main = main, ...)
  
  #if(!is.null(labels.loc)){
  #	if(legend.loc =="under"){ 	# draw the legend under the chart
  #        par(mar=c(0,2,0,1)+.1) 	# set the margins of the second panel            
  #		if(object.rows <3)
  #			ncol= object.rows
  #		else
  #			ncol = 3
  #			label.name = paste(rownames(object)," (",round(object,2),"%)",sep="")
  #           legend("bottom", legend=label.name, cex = cex.legend, 
  #				   fill=colorset, ncol=ncol,
  #				   box.col=element.color, border = element.color)                
  #        } # if label.loc is null, then do nothing
  #}
  invisible(NULL)
}
### <---------------------------------------------------------------------->