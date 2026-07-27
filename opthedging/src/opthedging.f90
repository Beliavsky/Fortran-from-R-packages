! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
module opthedging
   use opthedging_iid, only : call_payoff, hedging_iid, put_payoff
   use opthedging_interpolation, only : interpolation1d, linear_interpolate_uniform
   use opthedging_kinds, only : dp
   use opthedging_rng, only : random_normal, seed_random
   use opthedging_types, only : hedging_result
   implicit none
   public
end module opthedging
