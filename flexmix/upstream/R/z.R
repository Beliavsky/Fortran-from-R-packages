#
#  Copyright (C) 2004-2016 Friedrich Leisch and Bettina Gruen
#  $Id: z.R 5079 2016-01-31 12:21:12Z gruen $
#

###**********************************************************
## Backward compatibility

## component model driver
FLXglm <- FLXMRglm
FLXglmFix <- FLXMRglmfix
FLXmclust <- FLXMCmvnorm
FLXbclust <- FLXMCmvbinary

## concomitant model driver
FLXmultinom <- FLXPmultinom
FLXconstant <- FLXPconstant
