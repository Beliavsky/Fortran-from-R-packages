## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 6,
  fig.align = 'center'
)

## -----------------------------------------------------------------------------
# Load libraries
library(mstate)
library(ggplot2)

# Set general ggplot2 theme
theme_set(theme_bw(base_size = 14))

## ----load_dat-----------------------------------------------------------------
# Load data
data("aidssi")
head(aidssi)

# Shorter name
si <- aidssi

# Prepare transition matrix
tmat <- trans.comprisk(2, names = c("event-free", "AIDS", "SI"))
tmat

# Run msprep
si$stat1 <- as.numeric(si$status == 1)
si$stat2 <- as.numeric(si$status == 2)

silong <- msprep(
  time = c(NA, "time", "time"), 
  status = c(NA, "stat1", "stat2"), 
  data = si, 
  keep = "ccr5", 
  trans = tmat
)

# Run cox model
silong <- expand.covs(silong, "ccr5")

c1 <- coxph(Surv(time, status) ~ ccr5WM.1 + ccr5WM.2 + strata(trans),
            data = silong)

## ----msfit_prep---------------------------------------------------------------
# Data to predict
WW <- data.frame(
  ccr5WM.1 = c(0, 0),
  ccr5WM.2 = c(0, 0), 
  trans = c(1, 2), 
  strata = c(1, 2)
)

# Make msfit object
msf.WW <- msfit(
  object = c1, 
  newdata = WW, 
  trans = tmat
)

## ----msfitplot_base_1---------------------------------------------------------
plot(msf.WW)

## ----msfitplot_ggplot2_1------------------------------------------------------
plot(msf.WW, use.ggplot = TRUE)

## ----msfitplot_base_2---------------------------------------------------------
par(mfrow = c(1, 2))
plot(msf.WW, type = "separate")

## ----msfitplot_ggplot_2-------------------------------------------------------
# Fixed scales
plot(msf.WW, type = "separate", use.ggplot = TRUE, scale_type = "fixed")

# Free scales
plot(msf.WW, type = "separate", use.ggplot = TRUE, scale_type = "free", xlim = c(0, 15))

## ----msfitplot_customs--------------------------------------------------------
par(mfrow = c(1, 1))
# A base R customised plot
plot(
  msf.WW, 
  type = "single", 
  cols = c("blue", "black"), # or numeric e.g. c(1, 2)
  xlim = c(0, 15),
  legend.pos = "top",
  lty = c("dashed", "solid"),
  use.ggplot = FALSE
)
title("Cumulative baseline hazards")

# A ggplot2 customised plot
plot(
  msf.WW, 
  type = "single", 
  cols = c("blue", "black"), # or numeric e.g. c(1, 2)
  xlim = c(0, 15),
  lty = c("dashed", "solid"),
  legend.pos = "bottom",
  use.ggplot = TRUE
) +
  
  # Add title and center
  ggtitle("Cumulative baseline hazards") +
  theme(plot.title = element_text(hjust = 0.5))

## ----msfitplot_confint1-------------------------------------------------------
plot(
  msf.WW, 
  type = "single", 
  use.ggplot = TRUE,
  conf.type = "log", 
  conf.int = 0.95
)

## ----msfitplot_confint2-------------------------------------------------------
plot(
  msf.WW, 
  type = "separate", 
  use.ggplot = TRUE,
  conf.type = "log", 
  conf.int = 0.95,
  scale_type = "free_y"
)

## -----------------------------------------------------------------------------
# Run probtrans
pt.WW <- probtrans(msf.WW, predt = 0)

# Example predict at different times
summary(pt.WW, times = c(1, 5, 10))

## ----plot_pt_filled1----------------------------------------------------------
plot(pt.WW, from = 1)

## ----plot_pt_filled2----------------------------------------------------------
# from = 1 implied
plot(pt.WW, use.ggplot = TRUE)

## ----plot_pt_filled3----------------------------------------------------------
plot(pt.WW, use.ggplot = TRUE, label = "annotate", cex = 6)   

## ----plot_pt_filled4----------------------------------------------------------
# Check state order again from transition matrix
tmat

# Plot aids (state 2), then event-free (state 1), and SI on top (state 3)
plot(pt.WW, use.ggplot = TRUE, ord = c(2, 1, 3))   

## ----plot_pt_stacked----------------------------------------------------------
plot(pt.WW, use.ggplot = TRUE, type = "stacked")   

## ----plot_pt_separate1--------------------------------------------------------
par(mfrow = c(1, 3))
plot(pt.WW, type = "separate")   

## ----plot_pt_separate2--------------------------------------------------------
plot(
  pt.WW, 
  use.ggplot = TRUE, 
  type = "separate",
  conf.int = 0.95, # 95% level
  conf.type = "log"
)   

## ----plot_pt_single1----------------------------------------------------------
plot(
  pt.WW, 
  use.ggplot = TRUE, 
  type = "single",
  conf.int = 0.95, # 95% level
  conf.type = "log"
)   

## ----plot_pt_single2----------------------------------------------------------
plot(
  pt.WW, 
  use.ggplot = TRUE, 
  type = "single",
  conf.type = "none",
  lty = c(1, 2, 3), # change the linetype
  lwd = 1.5, 
)   

## ----plot_pt_single3----------------------------------------------------------
# Run plot and extract data using $data
dat_plot <- plot(x = pt.WW, use.ggplot = TRUE, type = "single")$data
  
# Begin new plot - Exclude or select states to be plotted
ggplot(data = dat_plot[state != "event-free", ], aes(
  x = time, 
  y = prob, 
  ymin = CI_low, 
  ymax = CI_upp,
  group = state,
  linetype = state,
  col = state
)) +
  
  # Add CI and lines; change fill = NA to remove CIs
  geom_ribbon(col = NA, fill = "gray", alpha = 0.5) +
  geom_line() +
  
  # Remaining details
  labs(x = "Time", y = "Cumulative Incidence") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 14), expand = 0)

## ----plot.Cuminc1-------------------------------------------------------------
cum_incid <- Cuminc(
  time = "time",
  status = "status",
  data = si
)

plot(
  x = cum_incid,
  use.ggplot = TRUE,
  conf.type = "log",
  lty = 1:2,
  conf.int = 0.95,
)

## ----plot.cuminc_x------------------------------------------------------------
cum_incid_grp <- Cuminc(
  time = "time",
  status = "status",
  group = "ccr5",
  data = si
)

plot(
  x = cum_incid_grp,
  use.ggplot = TRUE,
  conf.type = "none",
  lty = 1:4, 
  facet = FALSE
)

plot(
  x = cum_incid_grp,
  use.ggplot = TRUE,
  conf.type = "none",
  lty = 1:4, 
  facet = TRUE
)

## ----vis_multiple_prep--------------------------------------------------------
# 1. Prepare patient data - both CCR5 genotypes
WW <- data.frame(
  ccr5WM.1 = c(0, 0), 
  ccr5WM.2 = c(0, 0), 
  trans = c(1, 2), 
  strata = c(1, 2)
)

WM <- data.frame(
  ccr5WM.1 = c(1, 0), 
  ccr5WM.2 = c(0, 1),
  trans = c(1, 2), 
  strata = c(1, 2)
)

# 2. Make msfit objects
msf.WW <- msfit(c1, WW, trans = tmat)
msf.WM <- msfit(c1, WM, trans = tmat)

# 3. Make probtrans objects
pt.WW <- probtrans(msf.WW, predt = 0)
pt.WM <- probtrans(msf.WM, predt = 0)

## ----vis_multiple_plot--------------------------------------------------------
vis.multiple.pt(
  x = list(pt.WW, pt.WM), 
  from = 1,
  to = 2, 
  conf.type = "log",
  cols = c(1, 2),
  labels = c("Pat WW", "Pat WM"),
  legend.title = "Ref patients"
)

## ----vis_multiple_plot2-------------------------------------------------------
vis.multiple.pt(
  x = list(pt.WW), 
  from = 1,
  to = 2, 
  conf.type = "log",
  cols = c(1, 2),
  labels = c("Pat WW", "Pat WM"),
  legend.title = "Ref patients"
)

## ----mirror1------------------------------------------------------------------
vis.mirror.pt(
  x = list(pt.WW, pt.WM),
  titles = c("WW", "WM"),
  horizon = 10
)

## ----mirror2------------------------------------------------------------------
vis.mirror.pt(
  x = list(pt.WW, pt.WM),
  titles = c("WW", "WM"),
  size_titles = 8,
  horizon = 5,
  breaks_x_left = c(0, 1, 2, 3, 4, 5),
  breaks_x_right = c(0, 1, 2, 3, 4),
  ord = c(3, 2, 1)
)

## ----saving, eval=FALSE-------------------------------------------------------
#  # Saving a ggplot2 plot
#  p <- plot(pt.WW, use.ggplot = TRUE)
#  ggplot2::ggsave("my_ggplot2_plot.png")
#  
#  # Standard graphics plot
#  png("my_standard_plot.png")
#  plot(pt.WW, use.ggplot = FALSE)
#  dev.off()

## ----session-info, include=TRUE, echo=TRUE, results='markup'------------------
# Date/time
Sys.time()

# Environment
sessionInfo()

