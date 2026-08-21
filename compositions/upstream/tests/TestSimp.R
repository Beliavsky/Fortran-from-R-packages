options(warn=1)
# ad-hod function for comparing numbers
qround <- function(x) round(x, dig=3)

library("compositions")

#js <- read.table("juraset.dat",skip=17,header=TRUE)
#js$Land <- factor( c("Wald","Weide","Wiese","Acker")[js$Land] )
#js$Rock <- factor( c("Argovian","Kimmeridgian","Sequanian","Portlandian","Quaternary")[js$Rock] )

#Rock <- js$Rock
#Land <- js$Land

#cdata <- js[,c("Cd","Cu","Pb","Co","Cr","Ni","Zn")]
#cdata <- data.matrix(js[,c("Cd","Cu","Pb","Co","Cr","Ni","Zn")])

#cd1 <- cdata[,1:3]
#cd2 <- cdata[,4:7]
data(SimulatedAmounts)
cdata <- sa.groups5
Land  <- sa.groups5.area
cd1 <- cdata[,1:3]
cd2 <- cdata[,4:5]

# Transformations
# clo 
clo(c(1,2,3))
clo(matrix(1:4,ncol=2))
data(iris)
clo(iris[1:10,1:4])
clo(0.5)
clo(matrix(0.5))
clo(matrix(0.5,nrow=5))
clo(matrix(0.5,ncol=5))

clo(iris[1:10,],c("Sepal.Length","Petal.Length"))
clo(iris[1:10,],c(2,3))

checker <- function(x,y) {
  x<-unclass(x)
  y<-unclass(y)
  if( sum(c(x-y)^2) > 1E-10 )
    stop("Wrong results in ", deparse(substitute(x)))
  x
}

qround( clr(cdata[1:10,]) )

qround( checker( clrInv(clr(cdata)) , clo(cdata) ) )
qround( ilr(cdata) )
qround( checker( ilrInv(ilr(cdata)) , clo(cdata) ) )
qround( alr(cdata) )
qround( checker( alrInv(alr(cdata)) , clo(cdata) ) )
qround( cpt(cdata) )
qround( checker( cptInv(cpt(cdata)) , clo(cdata) ) )
qround( ipt(cdata) )
qround( checker( iptInv(ipt(cdata)) , clo(cdata) ) )
qround( apt(cdata) )
qround( checker( aptInv(apt(cdata)) , clo(cdata) ) )
qround( ilt(cdata) )
qround( checker( iltInv(ilt(cdata)) , cdata ) )
qround( iit(cdata) )
qround( checker( iitInv(iit(cdata)) , cdata ) )

qround( clr(c(a=1,2,3)))
qround( ilr(c(a=1,2,3)))
qround( alr(c(a=1,2,3)))
qround( cpt(c(a=1,2,3)))
qround( ipt(c(a=1,2,3)))
qround( apt(c(a=1,2,3)))
qround( ilt(c(a=1,2,3)))
qround( iit(c(a=1,2,3)))
qround( checker( clrInv(clr(c(a=1,2,3))) , clo(c(1,2,3))) )
qround( checker( ilrInv(ilr(c(a=1,2,3))) , clo(c(1,2,3))) )
qround( checker( alrInv(alr(c(a=1,2,3))) , clo(c(1,2,3))) )
qround( checker( cptInv(cpt(c(a=1,2,3))) , clo(c(1,2,3))) )
qround( checker( iptInv(ipt(c(a=1,2,3))) , clo(c(1,2,3))) )
qround( checker( aptInv(apt(c(a=1,2,3))) , clo(c(1,2,3))) )
qround( checker( iltInv(ilt(c(a=1,2,3))) , c(1,2,3)) )
qround( checker( iitInv(iit(c(a=1,2,3))) , c(1,2,3)) )

# mean

qround( mean(acomp(cdata)) )
qround( mean(rcomp(cdata)) )
qround( mean(aplus(cdata)) )
qround( mean(rplus(cdata)) )

qround( meanCol(cdata) )
qround( meanCol(clo(cdata)) )
qround( clo(meanCol(cdata)) )

# var (Variation Matrix)
qround( var(rcomp(cdata)) )
qround( var(acomp(cdata)) )
qround( var(aplus(cdata)) )
qround( var(rplus(cdata)) )

# clr
qround( clr(mean(acomp(cdata))) )
qround( meanCol(clr(cdata)) )
qround( clrInv(meanCol(clr(cdata))) )
qround( mean(acomp(cdata)) )

# ilr
qround( ilr(mean(acomp(cdata))) )
qround( meanCol(ilr(cdata)) )
qround( ilrInv(meanCol(ilr(cdata))) )
qround( mean(acomp(cdata)) )

# alr
qround( alr(mean(acomp(cdata))) )
qround( meanCol(alr(cdata)) )
qround( alrInv(meanCol(alr(cdata))) )
qround( mean(acomp(cdata)) )

# Operations
mean(acomp(3 * (cdata - mean(acomp(cdata)))))


# barplot
barplot(acomp(cdata[1:10,]))
barplot(rcomp(cdata[1:10,]))
barplot(aplus(cdata[1:10,]))
barplot(rplus(cdata[1:10,]))

barplot(mean(acomp(cdata)))
barplot(mean(rcomp(cdata)))
barplot(mean(aplus(cdata)))
barplot(mean(rplus(cdata)))

# piechart
pie(mean(acomp(cdata)))
pie(mean(rcomp(cdata)))
pie(mean(aplus(cdata)))



# Triangular Diagrams
plot(acomp(cdata[,1:3]))
qround( mean(acomp(cdata[,1:3])) )

# In this file we suppress unnecessary warnings targeting end-users
# Human: if you want to see them, replace suppressWarnings() by I()
suppressWarnings( plot(acomp(cdata),margin="rcomp",pca=TRUE) )
plot(acomp(cdata),margin="acomp",pca=TRUE)
plot(acomp(cdata),margin="Cd",pca=TRUE) # bug corrected
qround( mean(acomp(cdata)) )

boxplot(acomp(cdata))              # boxplotscale
suppressWarnings( boxplot(acomp(cdata),Land,notch=TRUE) )# notch
#boxplot(acomp(cdata),Rock)

boxplot(acomp(cdata),log=FALSE)
boxplot(acomp(cdata),Land,log=FALSE)
#boxplot(acomp(cdata),Rock,log=FALSE)

boxplot(acomp(cdata),log=FALSE,ylim=c(0,5))
boxplot(acomp(cdata),Land,log=FALSE,ylim=c(0,5))
#boxplot(acomp(cdata),Rock,log=FALSE,ylim=c(0,5))

qqnorm(acomp(cdata),alpha=100)
qqnorm(acomp(cdata),alpha=0.05)
qqnorm(acomp(cdata[,-3]),alpha=0.05)
qqnorm(acomp(cdata[,-3]),alpha=0.05)
 
#boxplot(acomp(cdata[,-1]),js$Cd)
plot(Land,data.matrix(cdata)%*% rep(1,ncol(cdata)))


# rcomp.plots

boxplot(rcomp(cdata))
boxplot(rcomp(cdata),Land)
#boxplot(rcomp)(cdata,Rock)

qqnorm(rcomp(cdata))
qqnorm(rcomp(cdata),alpha=0.05)
qqnorm(rcomp(cdata[,1:3]),alpha=0.05)
plot(acomp(cdata[,1:3]))
ellipses(mean(acomp(cdata[,1:3])), var(acomp(cdata,1:3)),col="red",r=2)

ellipses(mean(rcomp(cdata[,1:3])), var(rcomp(cdata[,1:3])),col="blue",r=2)


plot(rplus(cdata[,1:2]))
ellipses(rplus(mean(rplus(cdata[,1:2]))), var(rplus(cdata[,1:2])),col="blue",r=2)
ellipses(aplus(mean(aplus(cdata[,1:2]))), var(aplus(cdata[,1:2])),col="red",r=2)

plot(aplus(cdata[,1:2]))
ellipses(aplus(mean(aplus(cdata[,1:2]))), var(aplus(cdata[,1:2])),col="red",r=2)
ellipses(rplus(mean(rplus(cdata[,1:2]))), var(rplus(cdata[,1:2])),col="blue",r=2)


straight(acomp(c(1,1,1)),c(2,1,3))


boxplot(rcomp(cdata[,-1]),cdata[,"Cd"])

# biplot, princomp

biplot(princomp(cdata))
biplot(princomp(acomp(cdata)))
biplot(princomp(rcomp(cdata)))
biplot(princomp(aplus(cdata)),choice=c(2,3))
biplot(princomp(rplus(cdata)),choice=c(2,3))

summary(princomp(cdata))
summary(princomp(acomp(cdata)))
summary(princomp(rcomp(cdata)))
summary(princomp(aplus(cdata)))
summary(princomp(rplus(cdata)))


# The following lines produce loadings, which are known to be defined
#   only up to a sign; to ensure meaningful comparison, they are
#   converted to absolute values
# Human: if you want to see the actual numbers, replace abs() by I()
abs( loadings(princomp(cdata)) )
abs( loadings(princomp(acomp(cdata))) )
abs( loadings(princomp(rcomp(cdata))) )
abs( loadings(princomp(aplus(cdata))) )
abs( loadings(princomp(rplus(cdata))) )


#names
qround( meanCol(cdata[,1:3]) )
qround( mean(acomp(cdata[,1:3])) )
oneOrDataset(c(a=1,b=2,c=3))

# covariance
qround( cov(acomp(cd1),acomp(cd2)) )
qround( cov(rcomp(cd1),rcomp(cd2)) )
qround( cov(aplus(cd1),aplus(cd2)) )
qround( cov(rplus(cd1),rplus(cd2)) )

# tmp <- princov(acomp(cd1,cd2))
# The following line produce loadings and scores, which are known to be defined
#   only up to a sign; to ensure meaningful comparison, they are
#   commented
# Human: if you want to see the actual numbers, uncomment next lines
#tmp

tmp <- princomp(acomp(cdata))
# tmp
plot(tmp)
biplot(tmp)
biplot(tmp,ch=2:3)

#pplot(acomp(cdata))

