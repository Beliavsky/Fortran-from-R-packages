
#1)------Value at Risk------

varg<-function(spec, alpha, ...)
{
        FF<-function (theta, ...) {do.call(paste("q",spec,sep=""),list(theta, ...))}
        var=FF(alpha, ...)
        return(var)

}


## EXAMPLE
#varg("norm", 0.9) 


#2)------Expected Shortfall------

esg<-function(spec, alpha, ...)
{
        FF<-function (theta, ...) {do.call(paste("q",spec,sep=""),list(theta, ...))}
        es = alpha
    	for (i in 1:length(alpha)) {
        es[i] = (1/alpha[i]) * integrate(FF, lower = 0, upper = alpha[i], 
            stop.on.error = FALSE)$value}
    	return(es)
}


## EXAMPLE
#esg("norm", 0.9)



#3)---------Tail Conditional Median-----------------

tcm<-function(spec, alpha, ...)
{
        FF<-function (theta, ...) {do.call(paste("p",spec,sep=""),list(theta, ...))}
        tcm=alpha
        for (i in 1:length(alpha))
        {
           FF2<-function (theta, ...) {(FF(theta, ...)-alpha[i])/(1-alpha[i])-0.5}
           tt=varg(spec, alpha[i], ...)
           tcm[i]=uniroot(FF2,lower=tt,upper=1000)$root
        }
    	return(tcm)
}


## EXAMPLE
#tcm("norm", 0.9)



#4)---------Expectiles-----------------

expp<-function(spec, alpha, a, b, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        expp=alpha
        for (i in 1:length(alpha))
        {
           FF4<-function (ell)
           {FF2<-function (theta, ...) {(max(ell-theta,0))**2*FF(theta, ...)}
           FF3<-function (theta, ...) {(max(theta-ell,0))**2*FF(theta, ...)}
           tt=alpha[i]*integrate(FF2,lower=a,upper=b)$value
           tt=tt-(1-alpha[i])*integrate(FF3,lower=a,upper=b)$value
           return(tt)}
           expp[i]=uniroot(FF4,lower=-1000,upper=1000)$root
        }
    	return(expp)
}


## EXAMPLE
#expp("norm", 0.9, a=-Inf, b=Inf)





#5)---------Beyond Value at Risk-----------------

bvar<-function(spec, alpha, a, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        bvar=alpha
        for (i in 1:length(alpha))
        {
           FF4<-function (ell, ...) {ell*FF(ell, ...)}
           tt=varg(spec, alpha[i], ...)
           bvar[i]=(integrate(FF4,lower=a,upper=tt)$value)/alpha[i]
        }
    	return(bvar)
}


## EXAMPLE
#bvar("norm", 0.9, a=-Inf) 




#6)------Expected Proporional Shortfall------

## If X is a continous random variable

epsg<-function(spec, alpha, ...)
{
        eps=alpha
        for (i in 1:length(alpha))
        {
       	   eps[i]=(1-alpha[i])*(esg(spec, alpha[i],...)/varg(spec, alpha[i], ...))
        }
        return(eps)

}

## EXAMPLE
#epsg("norm", 0.9)




#7)------Expectation------

## If X is a continous random variable

expect<-function(spec, a, b, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        FF4<-function (ell, ...) {ell*FF(ell, ...)}
        expect=integrate(FF4,lower=a,upper=b)$value
        return(expect)
}

## EXAMPLE
#expect("norm", -Inf, Inf)


#8)------Elementary risk measure------

## If X is a continous random variable

expvar<-function(spec, alpha, a, b, ...)
{
        tt=expect(spec, a, b, ...)
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        FF4<-function (ell, ...) {ell*ell*FF(ell, ...)}
        expvar=tt+alpha*sqrt((integrate(FF4,lower=a,upper=b)$value)-tt*tt)
        return(expvar)
}

## EXAMPLE
#expvar("norm", 0.9, -Inf, Inf)


#9)------Omega Risk Measure------

omegag<-function(spec, alpha, a, b, ...)
{
        FF<-function (theta, ...) {do.call(paste("p",spec,sep=""),list(theta, ...))}
        FF2<-function (theta, ...) {1-FF(theta,...)}
        omega=alpha
    	for (i in 1:length(alpha)) {
        omega[i] =  integrate(FF2, lower = alpha[i], upper = b, 
            stop.on.error = FALSE)$value /integrate(FF, lower = a, upper = alpha[i], 
            stop.on.error = FALSE)$value }
    	return(omega)
}


## EXAMPLE
#omegag("norm", 2, -Inf, Inf)


#10)------Sortino Ratio risk Measure------



sortinog<-function(spec, alpha, a, b, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
	sortino = alpha
        tt=expect(spec, a, b, ...)
    	for (i in 1:length(alpha)) {
	FF3<-function (theta, ...) {(alpha[i]-theta)**2*FF(theta,...)}
	sortino[i]=(tt-alpha[i])*(integrate(FF3, lower = a, upper = alpha[i], stop.on.error = FALSE)$value)**(-1/2)}
    	return(sortino)
}

## EXAMPLE
#sortinog("norm", 2, -Inf, Inf)


#11)------Kappa risk Measure------



kappag<-function(spec, alpha, n, a, b, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
	kappa = alpha
        tt=expect(spec, a, b, ...)
    	for (i in 1:length(alpha)) {
	FF3<-function (theta, ...) {(alpha[i]-theta)**n*FF(theta,...)}
	kappa[i]=(tt-alpha[i])*(integrate(FF3, lower = a, upper = alpha[i], stop.on.error = FALSE)$value)**(-1/n)}
    	return(kappa)
}

## EXAMPLE
#kappag("norm", 2, 5, -Inf, Inf)



#12)------Wang's risk Measure------

wangg1<-function(spec, alpha, a, b, ...)
{
        FF<-function (theta, ...) {do.call(paste("p",spec,sep=""),list(theta, ...))}
        wang=alpha
        for (i in 1:length(alpha))
        {
          FF2<-function (theta, ...) {(1-FF(theta,...))**(alpha[i])-(1-FF(theta,...))}
          FF3<-function (theta, ...) {(FF(theta,...))**(alpha[i])-FF(theta,...)}
          tt=0.5*integrate(FF2, lower = a, upper = b, stop.on.error = FALSE)$value
          tt=tt+0.5*integrate(FF3, lower = a, upper = b, stop.on.error = FALSE)$value
          wang[i]=tt
        }
    	return(wang)
}


## EXAMPLE
#wangg1("lnorm", 0.9, 0, Inf)



#13)------Wang's risk Measure------

wangg2<-function(spec, alpha, a, b, ...)
{
        FF<-function (theta, ...) {do.call(paste("p",spec,sep=""),list(theta, ...))}
        wang=alpha
        for (i in 1:length(alpha))
        {
          FF2<-function (theta, ...) {(1-FF(theta,...))**(alpha[i])-(1-FF(theta,...))}
          tt=integrate(FF2, lower = a, upper = b, stop.on.error = FALSE)$value
          wang[i]=tt
        }
    	return(wang)
}


## EXAMPLE
#wangg2("lnorm", 0.9, 0, Inf)


#14)------Stones's risk Measures------


stoneg1<-function(spec, x0, k, a, b, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
    	FF2<-function (theta, ...) {(abs(theta-x0))**k*FF(theta,...)}
        stone=integrate(FF2, lower = a, upper = b, stop.on.error = FALSE)$value
    	return(stone)
}

## EXAMPLE
#stoneg1("norm", 8, 3, -Inf, Inf)



#15)------Stones's risk Measures------


stoneg2<-function(spec, x0, k, a, b, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
    	FF2<-function (theta, ...) {(abs(theta-x0))**k*FF(theta,...)}
        stone=(integrate(FF2, lower = a, upper = b, stop.on.error = FALSE)$value)**(1/k)
    	return(stone)
}

## EXAMPLE
#stoneg2("norm", 8, 3, -Inf, Inf)



#16)------Luces risk measure------


luceg1<-function(spec, a, b, aa, bb, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        FF2<-function (theta, ...) {log(FF(theta,...))*FF(theta,...)}
        tt=bb-aa*integrate(FF2,lower = a, upper = b)$value
    	return(tt)
}

## EXAMPLE
#luceg1("unif", 0, 1, 1, 0)


#17)------Luces risk measure------


luceg2<-function(spec, a, b, aa, bb, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        FF2<-function (theta, ...) {(FF(theta,...))**(1-bb)}
        tt=aa*integrate(FF2,lower = a, upper = b)$value
    	return(tt)
}

## EXAMPLE
#luceg2("unif", 0, 1, 1, 0)


#18)------Luces risk measure------


luceg3<-function(spec, a, b, aa, bb, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        FF2<-function (theta, ...) {log(theta)*FF(theta,...)}
        tt=bb+aa*integrate(FF2,lower = a, upper = b)$value
    	return(tt)
}

## EXAMPLE
#luceg3("unif", 0, 1, 1, 0)


#19)------Luces risk measure------


luceg4<-function(spec, a, b, aa, bb, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        FF2<-function (theta, ...) {theta**bb*FF(theta,...)}
        tt=aa*integrate(FF2,lower = a, upper = b)$value
    	return(tt)
}

## EXAMPLE
#luceg4("norm",-Inf, Inf, 1, 0)

#20)------Sarin risk measure------


saring1<-function(spec, a, b, k, c, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        FF2<-function (theta, ...) {exp(c*theta)*FF(theta,...)}
        tt=k*integrate(FF2, lower = a, upper = b)$value
    	return(tt)
}


## EXAMPLE
#saring1("norm", -Inf, Inf, 1, 0)


#21)------Sarin risk measure------


saring2<-function(spec, a, b, aa, bb1, bb2, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        FF2<-function (theta, ...) {log(abs(theta))*FF(theta,...)}
        FF3<-function (theta, ...) {(log(abs(theta)))**2*FF(theta,...)}
        tt=bb1*integrate(FF,lower = 0, upper = b)$value
        tt=tt+bb2*integrate(FF,lower = a, upper = 0)$value
        tt2=integrate(FF2,lower = a, upper = b)$value
        tt3=integrate(FF3,lower = a, upper = b)$value
        tt=tt+aa*tt2-(aa*aa/2)*(tt3-tt2*tt2)
    	return(tt)
}

## EXAMPLE
#saring2("norm",-Inf, Inf, 1, 1, 1)


#22)------Sarin risk measure------


saring3<-function(spec, a, b, aa, bb1, bb2, ...)
{
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        FF2<-function (theta, ...) {(abs(theta))**aa*FF(theta,...)}
        FF3<-function (theta, ...) {(abs(theta))**(2*aa)*FF(theta,...)}
        tt=bb1*integrate(FF2,lower = 0, upper = b)$value
        tt=tt+bb2*integrate(FF2,lower = a, upper = 0)$value
        tt2=integrate(FF2,lower = a, upper = b)$value
        tt3=integrate(FF3,lower = a, upper = b)$value
        tt=tt+(aa/(2*(aa-1)))*tt3-0.5*(tt3-tt2*tt2)
    	return(tt)
}


## EXAMPLE
#saring3("norm",-Inf, Inf, 1, 1, 1)


#23)------Bronshtein and KurelenKova's risk Measures------


BKg1<-function(spec, alpha, a, b, ...)
{
        tt=expect(spec, a, b, ...)
        tt2=varg(spec, alpha, ...)-tt
        return(tt2)
}


## EXAMPLE
#BKg1("norm", 0.9, -Inf, Inf)


#24)------Bronshtein and KurelenKova's risk Measures------


BKg2<-function(spec, alpha, a, b, ...)
{
        tt=expect(spec, a, b, ...)
        tt2=esg(spec, alpha, ...)-tt
        return(tt2)
}


## EXAMPLE
#BKg2("norm", 0.9, -Inf, Inf)



#25)------Bronshtein and KurelenKova's risk Measures------


BKg3<-function(spec, alpha, a, b, beta, ...)
{
        tt=expect(spec, a, b, ...)
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        FF2<-function (theta, ...) {(abs(theta-tt))*FF(theta,...)}
        tt2=esg(spec, alpha, ...)-beta*integrate(FF2,lower = a, upper = b)$value
        return(tt2)
}



## EXAMPLE
#BKg3("norm", 0.9, -Inf, Inf, 1)


#26)------Bronshtein and KurelenKova's risk Measures------


BKg4<-function(spec, alpha, a, b, beta, ...)
{
        tt=expect(spec, a, b, ...)
        FF<-function (theta, ...) {do.call(paste("d",spec,sep=""),list(theta, ...))}
        FF2<-function (theta, ...) {(abs(theta-tt))*FF(theta,...)}
        tt2=varg(spec, alpha, ...)+esg(spec, alpha, ...)-beta*integrate(FF2,lower = a, upper = b)$value
        return(tt2)
}



## EXAMPLE
#BKg4("norm", 0.9, -Inf, Inf, 1)

