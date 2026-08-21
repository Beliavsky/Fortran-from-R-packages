options(warn=1)
# ad-hod function for comparing numbers
qround <- function(x) round(x, dig=3)

library("compositions")
set.seed(1334345)
data(SimulatedAmounts)
par(pch=20)
mydata <- simulateMissings(sa.groups5,dl=0.01,knownlimit=TRUE,MAR=0.05,MNARprob=0.05,SZprob=0.05)


# standard way of creating a compositional dataset
# (cdata <- acomp(mydata))

# In this file we suppress unnecessary warnings targeting end-users
# Human: if you want to see them, replace suppressWarnings() by I()
cdata <- suppressWarnings( acomp(mydata) )
suppressWarnings( qround( mean(cdata) ) )
qround( mean(acomp(sa.groups5)) )
plot(acomp(sa.groups5))
suppressWarnings( plot(mean(cdata),add=T,col="red") )
plot(mean(acomp(sa.groups5)),add=T,col="green")
suppressWarnings( qround( mean(cdata - mean(cdata)) ) )

mm  <- suppressWarnings( mean(cdata) )
erg <- suppressWarnings(  var(cdata) )

print(erg)
lapply(lapply(svd(erg), abs), qround)

ellipses(mm,erg)
suppressWarnings(
  ellipses(mean(acomp(sa.groups5)),var(acomp(sa.groups5)))
  )

# cdata <- acomp(mydata)
# plot(cdata)
suppressWarnings( plot(cdata) )
suppressWarnings( plot(mm,add=T,col="blue") )
suppressWarnings( plot(mean(acomp(sa.groups5)),add=T,col="green") )
suppressWarnings( mean(cdata - mean(cdata)) )
# mm <-mean(cdata)
# erg <-var(cdata)
# svd(erg)
suppressWarnings( ellipses(mm,erg) )
suppressWarnings( ellipses(mm,erg,r=2) )



#cdata <- rcomp(mydata)
#cdata

#mean(cdata)
#mean(rcomp(sa.groups5))
#plot(rcomp(sa.groups5))
#plot(mean(cdata),add=T,col="red",pch=20)
#plot(mean(rcomp(sa.groups5)),add=T,col="green",pch=20)
#mean(cdata - mean(cdata)) # Nonsense because the difference is noncompositional


cdata <- aplus(mydata)
# cdata
qround( mean(cdata) )
qround( mean(aplus(sa.groups5)) )
plot(aplus(sa.groups5))
plot(mean(cdata),add=T,col="red",pch=20)
plot(mean(aplus(sa.groups5)),add=T,col="green",pch=20)
mean(cdata - mean(cdata))

cdata <- rplus(mydata)
# cdata
suppressWarnings( qround(  mean(cdata) ) )
qround(  mean(rplus(sa.groups5)) )
plot(rplus(sa.groups5))

# Attention: next row should send a warning and produce no result. 
# Human: uncomment if you want to check that manually
# plot(mean(cdata),add=T,col="red",pch=20)

plot(mean(rplus(sa.groups5)),add=T,col="green",pch=20)
# Attention: next row should send a warning and produce no result
# Human: uncomment if you want to check that manually
# mean(cdata - mean(cdata)) # Nonsense because the difference is non compositonal


suppressWarnings( plot(acomp(mydata)) )
plot(aplus(mydata))
suppressWarnings( plot(rcomp(mydata)) )
plot(rplus(mydata))

suppressWarnings( boxplot(acomp(mydata)) )
boxplot(aplus(mydata))
suppressWarnings( boxplot(rcomp(mydata)) )
boxplot(rplus(mydata))

suppressWarnings( barplot(acomp(mydata)) )
barplot(aplus(mydata))
suppressWarnings( barplot(rcomp(mydata)) )
barplot(rplus(mydata))
