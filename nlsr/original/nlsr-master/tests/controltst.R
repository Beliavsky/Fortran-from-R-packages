## Problem in 1 parameter to ensure methods work in trivial case
require(nlsr)
cat("Check nlsr.control\n")
cat("defaults\n")
ctrl <- nlsr.control()
pctrl(ctrl)
ctrl<-NULL
cat("Try control<-list(femax=4, phi=.12)\n")
control<-list(femax=4, phi=.12)
ctrl <- nlsr.control(control)
pctrl(ctrl)
ctrl<-NULL
cat("Try INVALID control<-list(fellmax=4, phi=.12) [bad name]\n")
control<-list(fellmax=4, phi=.12)
ctrl <- try(nlsr.control(control))
pctrl(ctrl)
