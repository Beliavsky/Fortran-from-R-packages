## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----geometries, echo=FALSE---------------------------------------------------
tbl = matrix(c("no", "no", "no", "no", "yes", "no",  "no", "no", "yes", "yes", "yes", "--",   
               "no", "yes","no", "yes","yes", "no",
               "real compositional", "relative compositional", "real amounts", "relative amounts", "count multinomial", "real multivariate",
               "rcomp", "acomp", "rplus", "aplus", "ccomp", "rmult"), ncol=5)
colnames(tbl) = c("counts?", "total?", "relative?", "geometry", "class")
knitr::kable(tbl)

## ----setup--------------------------------------------------------------------
library(compositions)

## ----setupData----------------------------------------------------------------
data("Hydrochem")
dim(Hydrochem)
colnames(Hydrochem)

## ----dataPreliminaries--------------------------------------------------------
xp = aplus(Hydrochem[, c("Ca","K", "Na", "Mg")])
summary(xp)

## ----dataPreliminaryFig, fig.height=10, fig.asp=1-----------------------------
plot(xp, col=Hydrochem$River, cex=0.5)

## -----------------------------------------------------------------------------
is.aplus(xp)    # check with S3 paradigm
is(xp, "aplus") # check with S4 paradigm

## -----------------------------------------------------------------------------
is(xp, "amounts") # only S4 paradigm

## -----------------------------------------------------------------------------
is(rplus(xp), "amounts")
is(ccomp(xp), "amounts")

## ----compPreliminaries--------------------------------------------------------
xc = acomp(Hydrochem[, c("Ca","K", "Na", "Mg")])
summary(xc)

## ----compPreliminaryFig, fig.height=10, fig.asp=1-----------------------------
plot(xc, col=Hydrochem$River, cex=0.5)

## -----------------------------------------------------------------------------
is.acomp(xc)    # check with S3 paradigm
is(xc, "acomp") # check with S4 paradigm

## -----------------------------------------------------------------------------
is(xc, "compositional") 

## -----------------------------------------------------------------------------
is(rcomp(xc), "compositional")
is(ccomp(xc), "compositional")

## ----eval=FALSE---------------------------------------------------------------
#  as(xc, "data.frame")
#  as(xc, "structure")

## -----------------------------------------------------------------------------
summary(xp$Ca)

## -----------------------------------------------------------------------------
summary(xp$clr)

## -----------------------------------------------------------------------------
summary(xp)
summary(cdt(xp))
summary(xp$cdt)
summary(xp$cdt$cdtInv)

## -----------------------------------------------------------------------------
xp0 = xp[1:5,]
xp0
xc0 = xc[1:5,]
xc0
xc[10,]

## -----------------------------------------------------------------------------
xp0[,"Ca"]
xp0[,c("Ca","K")]

## -----------------------------------------------------------------------------
xc0[,c("Ca","K")]

## -----------------------------------------------------------------------------
xc0$Ca
xc0[,"Ca"]
xp0[,"Ca"]

## -----------------------------------------------------------------------------
xc0[,c("Ca","K"), drop=TRUE]
xp0[,"Ca", drop=TRUE]

## -----------------------------------------------------------------------------
Hydrochem$comp <- xc
head(Hydrochem)
dim(Hydrochem)

## ----lm_example---------------------------------------------------------------
head(Hydrochem$comp)
codalm = lm(alr(comp)~River, data=Hydrochem) 
summary(codalm)

## ----lm_predict---------------------------------------------------------------
codalm = lm(comp$idt~River, data=Hydrochem) 
head(predict(codalm)$idtInv)

## ----lm_predict2--------------------------------------------------------------
head(idtInv(predict(codalm), orig=xc))

## ----selfbacktrafo------------------------------------------------------------
head(backtransform(idt(xp)))
head(backtransform(predict(codalm), as=codalm$residuals))

## ----eval=FALSE---------------------------------------------------------------
#  tbHydro = tibble::as_tibble(Hydrochem)
#  tbHydro$comp = xc
#  tbHydro$comp$Ca

## ----eval=FALSE---------------------------------------------------------------
#  data("jura", package="gstat")
#  coords = jura.pred[, 1:2]
#  compo = jura.pred[, 7:13]
#  xc = cdt(acomp(compo))
#  spcompo = SpatialPointsDataFrame(coords=coords, data=xc)
#  summary(spcompo@data)
#  class(spcompo@data)

