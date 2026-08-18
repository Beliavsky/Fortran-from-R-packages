! BiasedUrn-fortran
! Public umbrella module.
! Upstream BiasedUrn copyright (c) 2002-2024 Agner Fog.
! License: GNU General Public License version 3.
module biasedurn
   use biasedurn_kinds, only : dp
   use biasedurn_math, only : biasedurn_seed
   use biasedurn_fisher, only : dfnchypergeo, pfnchypergeo, qfnchypergeo, &
      rfnchypergeo, meanfnchypergeo, varfnchypergeo, modefnchypergeo, &
      oddsfnchypergeo, numfnchypergeo, minhypergeo, maxhypergeo
   use biasedurn_wallenius, only : dwnchypergeo, pwnchypergeo, qwnchypergeo, &
      rwnchypergeo, meanwnchypergeo, varwnchypergeo, modewnchypergeo, &
      oddswnchypergeo, numwnchypergeo
   use biasedurn_multifisher, only : dmfnchypergeo, rmfnchypergeo, &
      momentsmfnchypergeo, meanmfnchypergeo, varmfnchypergeo, &
      oddsmfnchypergeo, nummfnchypergeo, minmhypergeo, maxmhypergeo
   use biasedurn_multiwallenius, only : dmwnchypergeo, rmwnchypergeo, &
      momentsmwnchypergeo, meanmwnchypergeo, varmwnchypergeo, &
      oddsmwnchypergeo, nummwnchypergeo
   implicit none
   public
end module biasedurn
