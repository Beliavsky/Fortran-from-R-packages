pr.curve<-function( scores.class0, scores.class1=scores.class0, weights.class0=NULL, 
					weights.class1 = {if(is.null(weights.class0)){NULL}else{1-weights.class0}}, sorted = FALSE, curve = FALSE, 
					minStepSize=min(1,ifelse(is.null(weights.class0),1,sum(weights.class0)/100)),
					max.compute=F, min.compute=F, rand.compute=F, dg.compute=T){
	if(!sorted){
		o0<-order(scores.class0);
		scores.class0<-scores.class0[o0];
		if(!is.null(weights.class0)){
			weights.class0<-weights.class0[o0];
		}
		o1<-order(scores.class1);
		scores.class1<-scores.class1[o1];
		if(!is.null(weights.class1)){
			weights.class1<-weights.class1[o1];
		}
	}
	compute.pr(scores.class0,scores.class1,weights.class0,weights.class1,curve,minStepSize,max.compute,min.compute,rand.compute,dg.compute);
}




roc.curve<-function( scores.class0, scores.class1=scores.class0, weights.class0=NULL, 
					 weights.class1 = {if(is.null(weights.class0)){NULL}else{1-weights.class0}}, sorted = FALSE, curve = FALSE, 
					 max.compute=F, min.compute=F, rand.compute=F){
	if(!sorted){
		o0<-order(scores.class0);
		scores.class0<-scores.class0[o0];
		if(!is.null(weights.class0)){
			weights.class0<-weights.class0[o0];
		}
		o1<-order(scores.class1);
		scores.class1<-scores.class1[o1];
		if(!is.null(weights.class1)){
			weights.class1<-weights.class1[o1];
		}
	}
	compute.roc(scores.class0,scores.class1,weights.class0,weights.class1,curve,max.compute,min.compute,rand.compute);
}


check <- function( n, weights ) {
	if( !is.null( weights ) ) {
		if( n != length(weights) ) {
			stop( "The weights must have the same length as the scores." );
		}
		if( sum( weights < 0 ) != 0 ) {
			stop( "The weights must be non-negative." );
		}
	}
}



compute.pr <- function( sorted.scores.class0, sorted.scores.class1=sorted.scores.class0, weights.class0 = NULL, 
						weights.class1 = {if(is.null(weights.class0)){NULL}else{1-weights.class0}}, curve = FALSE, 
						minStepSize=min(1,ifelse(is.null(weights.class0),1,sum(weights.class0)/100)),
						max.compute=F, min.compute=F, rand.compute=F, dg.compute=FALSE ){

	check( length(sorted.scores.class0), weights.class0 );
	check( length(sorted.scores.class1), weights.class1 );
	
	if( !is.null(sorted.scores.class1) & ( length(sorted.scores.class0) != length(sorted.scores.class1) | 
										   suppressWarnings( sum(sorted.scores.class0 != sorted.scores.class1) > 0 ) 
	) & is.null(weights.class0) & is.null(weights.class1) ){
		weights.class0<-c(rep(1,length(sorted.scores.class0)),rep(0,length(sorted.scores.class1)));
		sorted.scores.class0<-c(sorted.scores.class0,sorted.scores.class1);
		o0<-order(sorted.scores.class0);
		sorted.scores.class0<-sorted.scores.class0[o0];
		weights.class0<-weights.class0[o0];
		weights.class1<-1-weights.class0;
		sorted.scores.class1<-sorted.scores.class0;
		
		all.scores<-sorted.scores.class0;
		all.weights.pos<-weights.class0;
		all.weights.neg<-weights.class1;
	}else{
		if(is.null(weights.class0)){
			weights.class0<-rep(1,length(sorted.scores.class0))
		}
		if(is.null(weights.class1)){
			weights.class1<-rep(1,length(sorted.scores.class1))
		}
	
		all.scores<-c(sorted.scores.class0,sorted.scores.class1);
		all.weights.pos<-c(weights.class0,rep(0,length(sorted.scores.class1)));
		all.weights.neg<-c(rep(0,length(sorted.scores.class0)),weights.class1);
	}
	
	
	davis.and.goadrich <- dg.compute & ( length(sorted.scores.class0) == length(sorted.scores.class1) & 
								suppressWarnings( sum( sorted.scores.class0 != sorted.scores.class1 ) == 0 ) & 
								length(weights.class0) == length(weights.class1) &
								suppressWarnings( sum( weights.class0 != (1 - weights.class1) ) == 0 ) &
								sum(weights.class0 != 0 & weights.class0 != 1)==0);
	
	o<-order(all.scores,decreasing = T);
	all.scores<-all.scores[o]
	all.weights.pos<-all.weights.pos[o];
	all.weights.neg<-all.weights.neg[o];
	
	cum.weights.pos<-cumsum(all.weights.pos);
	cum.weights.neg<-cumsum(all.weights.neg);
	cum.use<-c(all.scores[-length(all.scores)]!=all.scores[-1],TRUE)
	
	all.scores<-all.scores[cum.use]
	cum.weights.pos<-cum.weights.pos[cum.use];
	cum.weights.neg<-cum.weights.neg[cum.use];
	
	
	r.fg<-sum(all.weights.pos);
	tp<-cum.weights.pos
	fp<-cum.weights.neg;
	tp.prev<-c(0,cum.weights.pos[ -length(cum.weights.pos) ])
	fp.prev<-c(0,cum.weights.neg[ -length(cum.weights.neg) ])

	h<-(fp-fp.prev)/(tp-tp.prev);
	a<-1+h;
	b<-(fp.prev-h*tp.prev)/r.fg;
	h[tp==tp.prev]<-1;
	a[tp==tp.prev]<-1;
	b[tp==tp.prev]<-0;
	
	v<-( tp/r.fg - tp.prev/r.fg - b / a * ( log( a * tp/r.fg + b ) - log( a * tp.prev/r.fg + b ) ) ) / a;
	v2<-( tp/r.fg - tp.prev/r.fg ) / a;
	v[b==0]<-v2[b==0]
	
	vals<-v
	auc.integral<-sum(vals)
	auc.dg<-NA;

	if(davis.and.goadrich){
		
		min.mat<-cbind(tp.prev,tp,fp.prev,fp);

		idxs<-which(tp-tp.prev>1&tp/(tp+fp)!=tp.prev/(tp.prev+fp.prev));
		if(length(idxs)>0){
			m<-matrix(min.mat[-idxs,],ncol=ncol(min.mat));
		}else{
			m<-min.mat
		}
		
		auc.dg<-(m[,2]-m[,1])/r.fg * (m[,1]/(m[,1]+m[,3]) + m[,2]/(m[,2]+m[,4]))/2;
		if(is.nan(auc.dg[1])){
			auc.dg[1]<-(m[1,2]-m[1,1])/r.fg * (m[1,2]/(m[1,2]+m[1,4]));
		}

		diff<-tp-tp.prev;
		h2<-(fp-fp.prev)/diff;
		
		if(length(idxs)>0){
		
			temp.seq<-sapply(1:max(diff[idxs]),function(i){seq(0,i)})
				
			m<-sapply(idxs,function(i){
				x<-temp.seq[[tp[i]-tp.prev[i]]];
				prcs<-(tp.prev[i]+x)/(tp.prev[i]+x+fp.prev[i]+h2[i]*x)
				sum( 1/r.fg*(prcs[-1]+prcs[-length(prcs)])/2 )
	
			})
			auc.dg<-sum(c(auc.dg,m));
		}else{
			auc.dg<-sum(auc.dg);
		}
		
	}
	
	if(curve){
		minStepSize.2<-minStepSize/r.fg
		min.curve<-cbind(tp/r.fg,tp/(tp+fp),all.scores)
	#	print(min.curve)
		idxs<-which((tp-tp.prev)/r.fg>minStepSize.2 & tp/(tp+fp)!=tp.prev/(tp.prev+fp.prev))
	#	idxs<-which((tp-tp.prev)/r.fg>minStepSize)
		idxs<-idxs[ idxs>1 ];
		#print(idxs)
		if(length(idxs)>0){
			m<-sapply(idxs,function(i){
				x<-seq(0,min.curve[i,1]-min.curve[i-1,1],by = minStepSize.2);
				sns<-min.curve[i-1,1]+x;
				prcs<- ( min.curve[i-1,1]+x ) / ( min.curve[i-1,1]+x + fp[i-1]/r.fg + h[i]*x )
				temp<-rbind(sns,prcs,rep(all.scores[i],length(x)))
				temp
			})
			m<-matrix(unlist(m),ncol=3,byrow=T)
			
			m<-rbind(min.curve,m);
			
		}else{
			m<-min.curve
		}

		m<-m[ order(m[,1],-m[,3],decreasing = T), ]
		m<-matrix(m,ncol=3)

		m<-rbind(c(1,m[1,2:3]),m,c(0,m[nrow(m),2:3]))
		
		dimnames(m)<-c(NULL,NULL)
				
		res<-list( type = "PR", auc.integral = auc.integral, auc.davis.goadrich = auc.dg, curve=m );
	}else{
		res<-list( type = "PR", auc.integral = auc.integral, auc.davis.goadrich = auc.dg );
	}

	
	if(max.compute){
		scores0<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores0<-rep(1,length(sorted.scores.class0));
		}else{
			scores0<-weights.class0;
		}
		scores1<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores1<-rep(0,length(sorted.scores.class1));
		}else{
			scores1<-weights.class0;
		}
		
		max.res<-pr.curve( scores.class0=scores0, scores.class1=scores1,weights.class0=weights.class0,
						   weights.class1=weights.class1,curve=curve,minStepSize=minStepSize,dg.compute=dg.compute);
		res<-c(res,list(max=max.res));
	}
	
	if(min.compute){
		scores0<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores0<-rep(0,length(sorted.scores.class0));
		}else{
			scores0<-(-weights.class0);
		}
		scores1<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores1<-rep(1,length(sorted.scores.class1));
		}else{
			scores1<-(-weights.class0);
		}
		
		min.res<-pr.curve( scores.class0=scores0, scores.class1=scores1,weights.class0=weights.class0,
						   weights.class1=weights.class1,curve=curve,minStepSize=minStepSize,dg.compute=dg.compute);
		res<-c(res,list(min=min.res));
	}
	if(rand.compute){
		rand.auc<-NULL;
		if(is.null(weights.class0)){
			rand.auc<-length(sorted.scores.class0)/(length(sorted.scores.class0)+length(sorted.scores.class1));	
		}else{
			rand.auc<-sum(weights.class0)/sum(weights.class0+weights.class1);
		}
		rand.curve<-create.curve( 2 );
		rand.curve<-append.to.curve( rand.curve, c(0,rand.auc,0), 1 );
		rand.curve<-append.to.curve( rand.curve, c(1,rand.auc,0), 2 );
		rand.result<-list( type = "PR", auc.integral = rand.auc, auc.davis.goadrich = rand.auc, curve=rand.curve );
		class(rand.result)<-"PRROC";
		
		res<-c(res,list(rand=rand.result));
	}
	
	class(res)<-"PRROC";
	res
	
	
}




compute.roc<-function( sorted.scores.class0, sorted.scores.class1=sorted.scores.class0, weights.class0 = NULL, 
					   weights.class1 = {if(is.null(weights.class0)){NULL}else{1-weights.class0}}, curve = FALSE,
					   max.compute=F, min.compute=F, rand.compute=F){
	
	if( !is.null(sorted.scores.class1) & ( length(sorted.scores.class0) != length(sorted.scores.class1) | 
										   suppressWarnings( sum(sorted.scores.class0 != sorted.scores.class1) > 0 ) 
	) & is.null(weights.class0) & is.null(weights.class1) ){
		weights.class0<-c(rep(1,length(sorted.scores.class0)),rep(0,length(sorted.scores.class1)));
		sorted.scores.class0<-c(sorted.scores.class0,sorted.scores.class1);
		o0<-order(sorted.scores.class0);
		sorted.scores.class0<-sorted.scores.class0[o0];
		weights.class0<-weights.class0[o0];
		weights.class1<-1-weights.class0;
		sorted.scores.class1<-sorted.scores.class0;
		
		all.scores<-sorted.scores.class0;
		all.weights.pos<-weights.class0;
		all.weights.neg<-weights.class1;
	}else{
		if(is.null(weights.class0)){
			weights.class0<-rep(1,length(sorted.scores.class0))
		}
		if(is.null(weights.class1)){
			weights.class1<-rep(1,length(sorted.scores.class1))
		}
		
		all.scores<-c(sorted.scores.class0,sorted.scores.class1);
		all.weights.pos<-c(weights.class0,rep(0,length(sorted.scores.class1)));
		all.weights.neg<-c(rep(0,length(sorted.scores.class0)),weights.class1);
	}
	

	o<-order(all.scores,decreasing = T);
	all.scores<-all.scores[o]
	all.weights.pos<-all.weights.pos[o];
	all.weights.neg<-all.weights.neg[o];
	
	cum.weights.pos<-cumsum(all.weights.pos);
	cum.weights.neg<-cumsum(all.weights.neg);
	cum.use<-c(all.scores[-length(all.scores)]!=all.scores[-1],TRUE)
	
	all.scores<-all.scores[cum.use]
	cum.weights.pos<-cum.weights.pos[cum.use];
	cum.weights.neg<-cum.weights.neg[cum.use];
	
	r.fg<-sum(all.weights.pos);
	r.bg<-sum(all.weights.neg);
	
	sns<-c(0,cum.weights.pos/r.fg);
	fprs<-c(0,cum.weights.neg/r.bg);
	
	erg<-sum( (fprs[-1]-fprs[-length(fprs)]) * ( sns[-1]+sns[-length(sns)] )/2 );
	
	list.curve<-NULL;
	if(curve){
		list.curve<-cbind(fprs,sns,c(all.scores,all.scores[length(all.scores)]))
		dimnames(list.curve)<-c(NULL,NULL);
	}
	res<-list( type = "ROC", auc = erg, curve=list.curve );
	
	
	if(max.compute){
		scores0<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores0<-rep(1,length(sorted.scores.class0));
		}else{
			scores0<-weights.class0;
		}
		scores1<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores1<-rep(0,length(sorted.scores.class1));
		}else{
			scores1<-weights.class0;
		}
		
		max.res<-roc.curve( scores.class0=scores0, scores.class1=scores1,weights.class0=weights.class0,
							weights.class1=weights.class1,curve=curve);
		res<-c(res,list(max=max.res));
	}
	
	if(min.compute){
		scores0<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores0<-rep(0,length(sorted.scores.class0));
		}else{
			scores0<-(-weights.class0);
		}
		scores1<-NULL;
		if(is.null(weights.class0) | length(sorted.scores.class0)!=length(sorted.scores.class1) | suppressWarnings(sum(sorted.scores.class0!=sorted.scores.class1)) > 0){
			scores1<-rep(1,length(sorted.scores.class1));
		}else{
			scores1<-(-weights.class0);
		}
		
		min.res<-roc.curve( scores.class0=scores0, scores.class1=scores1,weights.class0=weights.class0,
							weights.class1=weights.class1,curve=curve);
		res<-c(res,list(min=min.res));
	}
	if(rand.compute){
		rand.auc<-0.5;
		rand.curve<-create.curve( 2 );
		rand.curve<-append.to.curve( rand.curve, c(0,0,0), 1 );
		rand.curve<-append.to.curve( rand.curve, c(1,1,0), 2 );
		rand.result<-list( type = "ROC", auc=rand.auc, curve=rand.curve );
		class(rand.result)<-"PRROC";
		
		res<-c(res,list(rand=rand.result));
	}
	
	class(res)<-"PRROC";
	res	
}



create.curve <- function( n ){
	m <- matrix( NA, nrow=n, ncol=3 );
	m
}

append.to.curve <- function( curve, p, row ){
	if( row>=nrow( curve ) ){
		curve2 <- matrix( NA, nrow=nrow( curve ) * 2, ncol=3 );
		curve2[ 1:nrow( curve ), ] <- curve;
		curve <- curve2;
	}
	curve[ row, ] <- p;
#	print(c(row,p))
#	if(is.nan(p[2])){
#		traceback(0)
#	}
	curve
}

print.PRROC<-function(x,...){
	if(x$type == "PR"){
		cat("\n  Precision-recall curve\n");
		cat("\n    Area under curve (Integral):\n");
		cat("    ",x$auc.integral,"\n");
		if( !is.null(x$max) & !is.null(x$min) ){
			cat("\n    Relative area under curve (Integral):\n");
			cat("    ",(x$auc.integral - x$min$auc.integral)/(x$max$auc.integral-x$min$auc.integral),"\n");
		}
		cat("\n    Area under curve (Davis & Goadrich):\n");
		if(!is.null(x$auc.davis.goadrich) & !is.na(x$auc.davis.goadrich)){
			cat("    ",x$auc.davis.goadrich,"\n");
			if( !is.null(x$max) & !is.null(x$min) ){
				cat("\n    Relative area under curves (Davis & Goadrich):\n");
				cat("    ",(x$auc.davis.goadrich - x$min$auc.davis.goadrich)/(x$max$auc.davis.goadrich-x$min$auc.davis.goadrich),"\n");
			}
		}else{
			cat("    cannot be computed for weighted data\n");
		}
		
	}else{
		cat("\n  ROC curve\n");
		cat("\n    Area under curve:\n");
		cat("    ",x$auc,"\n");
		if( !is.null(x$max) & !is.null(x$min) ){
			cat("\n    Relative area under curve:\n");
			cat("    ",(x$auc - x$min$auc)/(x$max$auc-x$min$auc),"\n");
		}
	}
	
	if(!is.null(x$curve)){
		cat("\n    Curve for scores from ",min(x$curve[,3])," to ",max(x$curve[,3]),"\n");
		cat("    ( can be plotted with plot(x) )\n\n");
	}else{
		cat("\n    Curve not computed ( can be done by using curve=TRUE )\n");
	}
	
	if(!is.null(x$max)){
		cat("\n\n    Maximum AUC:\n");
		if(x$type == "PR"){
			cat("    ",x$max$auc.integral," ",x$max$auc.davis.goadrich,"\n");
		}else{
			cat("    ",x$max$auc,"\n");
		}
	}
	
	if(!is.null(x$min)){
		cat("\n\n    Minimum AUC:\n");
		if(x$type == "PR"){
			cat("    ",x$min$auc.integral," ",x$min$auc.davis.goadrich,"\n");
		}else{
			cat("    ",x$min$auc,"\n");
		}
	}
	
	if(!is.null(x$rand)){
		cat("\n\n    AUC of a random classifier:\n");
		if(x$type == "PR"){
			cat("    ",x$rand$auc.integral," ",x$rand$auc.davis.goadrich,"\n");
		}else{
			cat("    ",x$rand$auc,"\n");
		}
	}
}


ggprroc<-function(x, auc.main=TRUE, auc.type=c("integral","davis.goadrich"),
                  xlab=NULL,ylab=NULL,
                  main=NULL,
                  max.plot = FALSE, min.plot = FALSE, rand.plot = FALSE,
                  fill.area = (max.plot & min.plot)){
  if( suppressWarnings( requireNamespace("ggplot2",quietly = T) ) ){
    auc.type<-match.arg(auc.type);
    if(is.null(x$curve)){
      stop("Curve is NULL. Use curve=T in pr.curve or roc.curve to obtain one.");
    }
    if(ncol(x$curve) != 3){
      stop("Curve has wrong dimension");
    }
    if(is.null(xlab)){
      my.xlab<-ifelse(x$type=="PR","Recall","FPR");
    }else{
      my.xlab<-xlab;
    }
    if(is.null(ylab)){
      my.ylab<-ifelse(x$type=="PR","Precision","Sensitivity");
    }else{
      my.ylab<-ylab;
    }
    
    if(is.null(main)){
      my.main<-paste(x$type," curve",sep="",collapse="");
    }else{
      my.main<-main;
    }
    if(auc.main){
      my.main<-paste(my.main,"\nAUC = ",format(ifelse(x$type=="PR",ifelse(auc.type=="integral",x$auc.integral,x$auc.davis.goadrich),x$auc)),sep="",collapse="");
    }
    
    curve <- data.frame(x$curve)
    colnames(curve)<-c("x","y","score")
    curve<-curve[nrow(curve):1,]
    
    pl<-ggplot2::ggplot(curve,ggplot2::aes(x=.data$x,y=.data$y))
    
    
    if( fill.area & !is.null(x$max) & !is.null(x$min)){
      xs<-c(x$min$curve[,1],x$max$curve[nrow(x$max$curve):1,1],x$min$curve[1,1]);
      ys<-c(x$min$curve[,2],x$max$curve[nrow(x$max$curve):1,2],x$min$curve[1,2]);
      pol.df<-data.frame(x=xs,y=ys)
      pl <- pl +
        ggplot2::geom_polygon(data=pol.df,ggplot2::aes(x=.data$x,y=.data$y),fill="grey95",linewidth=0)
      #polygon( x = xs, y = ys, density = -1, border = NA, col = fill.color );
    }
    
    if(max.plot & !is.null(x$max)){
      curve.max<-data.frame(x$max$curve)
      curve.max<-curve.max[nrow(curve.max):1,]
      colnames(curve.max)<-c("x","y","score")
      pl <- pl +
        ggplot2::geom_line(data=curve.max,ggplot2::aes(x=.data$x,y=.data$y),color="grey", linewidth=1,linetype="dashed")
    }
    if(min.plot & !is.null(x$min)){
      curve.min<-data.frame(x$min$curve)
      curve.min<-curve.min[nrow(curve.min):1,]
      colnames(curve.min)<-c("x","y","score")
      pl <- pl +
        ggplot2::geom_line(data=curve.min,ggplot2::aes(x=.data$x,y=.data$y),color="grey", linewidth=1,linetype="dotted")
    }
    if(rand.plot & !is.null(x$rand)){
      curve.rand<-data.frame(x$rand$curve)
      curve.rand<-curve.rand[nrow(curve.rand):1,]
      curve.rand<-curve.rand[!is.na(curve.rand[,1]),]
      colnames(curve.rand)<-c("x","y","score")
      pl <- pl +
        ggplot2::geom_line(data=curve.rand,ggplot2::aes(x=.data$x,y=.data$y),color="grey", linewidth=1,linetype="dotdash")
    }
    
    
    pl<- pl + ggplot2::geom_line(ggplot2::aes(color=.data$score),linewidth=1)+
      ggplot2::scale_color_viridis_c(option="H",name="Score")+
      ggplot2::xlim(0,1)+
      ggplot2::ylim(0,1)+
      ggplot2::xlab(my.xlab)+
      ggplot2::ylab(my.ylab)+
      ggplot2::ggtitle(my.main)+
      ggplot2::theme_bw()
    return(pl)
  }else{
    stop("ggplot2 required for calling ggprroc.")
  }
}


plot.PRROC<-function(x, xlim=c(0,1), ylim=c(0,1), auc.main=TRUE, auc.type=c("integral","davis.goadrich"), 
					 legend=ifelse(is.logical(color) & color==TRUE,4,NA), xlab=NULL, ylab=NULL, main=NULL, color=TRUE, lwd=3, 
					 add=FALSE, scale.color=hsv(h=seq(0,1,length=100)*0.8, s=1, v=1), 
					 max.plot = FALSE, min.plot = FALSE, rand.plot = FALSE, fill.area = (max.plot & min.plot),
					 maxminrand.col = grey(0.5), fill.color = grey(0.95),
					 ...){
	auc.type<-match.arg(auc.type);
	if(is.null(x$curve)){
		stop("Curve is NULL. Use curve=T in pr.curve or roc.curve to obtain one.");
	}
	if(ncol(x$curve) != 3){
		stop("Curve has wrong dimension");
	}
	if(is.null(xlab)){
		my.xlab<-ifelse(x$type=="PR","Recall","FPR");
	}else{
		my.xlab<-xlab;
	}
	if(is.null(ylab)){
		my.ylab<-ifelse(x$type=="PR","Precision","Sensitivity");
	}else{
		my.ylab<-ylab;
	}
	
	if(is.null(main)){
		my.main<-paste(x$type," curve",sep="",collapse="");
	}else{
		my.main<-main;
	}
	if(auc.main){
		my.main<-paste(my.main,"\nAUC = ",format(ifelse(x$type=="PR",ifelse(auc.type=="integral",x$auc.integral,x$auc.davis.goadrich),x$auc)),sep="",collapse="");
	}
 
	
	max.curve<-NULL;
	if(!is.null(x$max) & !is.null(x$max$curve)){
		max.curve<-x$max$curve;
	}
	min.curve<-NULL;
	if(!is.null(x$min) & !is.null(x$min$curve)){
		min.curve<-x$min$curve;
	}
	rand.curve<-NULL;
	if(!is.null(x$rand) & !is.null(x$rand$curve)){
		rand.curve<-x$rand$curve;
	}
	
	x<-x$curve;
	
	cols<-1;
	segment=F;
	plotscale.color=F;
	if( is.logical(color) ){
		if(color){
			min<-min(x[,3]);
			max<-max(x[,3]);
			if(min==max){
			  min <- min-min*0.1
			  max <- max+min*0.1
			}
	
			cols<-getColor( scale.color, x[,3], min, max );
			plotscale.color=T;
			segment=T;
		}else{
			cols<-1;
			segment<-F;
		}
	}else {
		cols<-color;
		segment<-F;
	}
	
	if(!add & !is.na(legend) & (is.numeric(legend) | suppressWarnings(legend==TRUE)) & plotscale.color ){
		if(is.logical(legend)){
			legend<-4;
		}
		m<-NULL;widths<-rep(1,2);heights<-rep(1,2)
		if(legend == 1){
			m<-matrix(c(1,2),nrow=2);
			heights<-c(4,lcm(2));
		}else if(legend==2){
			m<-matrix(c(2,1),nrow=1);
			widths=c(lcm(2.5),4);
		}else if(legend==3){
			m<-matrix(c(2,1),nrow=2);
			heights=c(lcm(2),4);
		}else{
			m<-matrix(c(1,2),nrow=1);
			widths=c(4,lcm(2.5));
		}
		layout(mat = m,widths = widths,heights = heights);
		
	}#else if(!add){
	#	layout(1);
	#}
	
	if(!add){
		plot(0,xlim=xlim,ylim=ylim,col=0,xlab=my.xlab,ylab=my.ylab,main=my.main,...);
	}
	
	if( !add ){
		if( fill.area & !is.null(max.curve) & !is.null(min.curve)){
			xs<-c(min.curve[,1],max.curve[nrow(max.curve):1,1],min.curve[1,1]);
			ys<-c(min.curve[,2],max.curve[nrow(max.curve):1,2],min.curve[1,2]);
			polygon( x = xs, y = ys, density = -1, border = NA, col = fill.color );
		}
		
		if(max.plot & !is.null(max.curve)){
			lines(max.curve[,1],max.curve[,2],col=maxminrand.col, lty="dashed", ...);
		}
		
		if(min.plot & !is.null(min.curve)){
			lines(min.curve[,1],min.curve[,2],col=maxminrand.col, lty="dotted", ...);
		}
		
		if(rand.plot & !is.null(rand.curve)){
			lines(rand.curve[,1],rand.curve[,2],col=maxminrand.col, lty="dotdash", ...);
		}
	}
	
	d=nrow(x);
	if( segment ) {
		segments( x[1:(d-1),1], x[1:(d-1),2], x[2:d,1], x[2:d,2], col=cols, lwd=lwd, ...);
	} else {
		lines( x[,1], x[,2], col=cols, lwd=lwd, ...);
	}
	
	if(!add & legend & !is.numeric(color) & color == TRUE){
		scale<-seq( min, max, length = 100 );
		cols<-getColor( scale.color, scale, min, max );
		bak<-par("mar");
		on.exit(par(mar=bak));
		if(legend==2 | legend==4){
			if(legend==4){par(mar=c(5,1,4,2)+0.1);}else{par(mar=c(5,2,4,1)+0.1);}
			image(c(1),scale,matrix(scale,nrow=1),col=cols,xlab="",ylab="",axes=F)
		}else{
			if(legend==1){par(mar=c(2,4,0,2)+0.1);}else{par(mar=c(0,4,2,2)+0.1);}
			image(scale,c(1),matrix(scale,ncol=1),col=cols,xlab="",ylab="",axes=F)
		}
		axis(legend)
		layout(1)
	}
	

}

getColor <- function( scale, x, min=min(x), max=max(x) )  {
	return( scale[round(1 + (length(scale)-1) * (x - min)/(max-min))] );
}
