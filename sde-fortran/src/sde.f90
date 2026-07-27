! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde
   use sde_kinds
   use sde_interfaces
   use sde_utils
   use sde_special
   use sde_random
   use sde_distributions
   use sde_linalg
   use sde_optimization
   use sde_models
   use sde_simulation
   use sde_density
   use sde_likelihood
   use sde_nonparametric
   use sde_change_point
   use sde_markov_distance
   use sde_estimating
   use sde_gmm
   use sde_information
   use sde_divergence
   implicit none
   public
end module sde
