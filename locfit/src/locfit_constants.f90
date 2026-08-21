! Constants translated from locfit src/lfcons.h; GPL-2-or-later.
module locfit_constants
  use locfit_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: mxdim = 15, mxdeg = 7
  real(dp), parameter, public :: pi = 3.141592653589793238462643_dp
  real(dp), parameter, public :: sqrt2 = 1.4142135623730950488_dp
  real(dp), parameter, public :: s2pi = 2.506628274631000502415765_dp
  real(dp), parameter, public :: logpi = 1.144729885849400174143427_dp
  real(dp), parameter, public :: golden = 0.61803398874989484820_dp
  real(dp), parameter, public :: hl2pi = 0.91893853320467267_dp
  real(dp), parameter, public :: sqrpi = 1.77245385090552_dp

  ! Kernel constants used by the upstream locfit C implementation.
  real(dp), parameter, public :: gfact = 2.5_dp
  real(dp), parameter, public :: efact = 3.0_dp
  real(dp), parameter, public :: huberc = 2.0_dp

  ! Adaptive fitting criteria.
  integer, parameter, public :: anone=0, acp=1, akat=2, amdi=3, aok=4

  ! Evaluation structures.
  integer, parameter, public :: enull=0, etree=1, ephull=2, edata=3, egrid=4
  integer, parameter, public :: ekdtr=5, ekdce=6, ecros=7, epres=8, exbar=9
  integer, parameter, public :: enone=10, esphr=11, efitp=50, espec=100

  ! Link functions.
  integer, parameter, public :: linit=0, ldefau=1, lcanon=2, lident=3
  integer, parameter, public :: llog=4, llogit=5, linver=6, lsqrt=7, lasin=8

  ! Likelihood work-vector components (Fortran 1-based equivalents).
  integer, parameter, public :: zlik=1, zmean=2, zdll=3, zddll=4, llen=4

  ! Kernels.
  integer, parameter, public :: wrect=1, wepan=2, wbisq=3, wtcub=4, wtrwt=5
  integer, parameter, public :: wgaus=6, wtria=7, wququ=8, w6cub=9, wminm=10
  integer, parameter, public :: wexpl=11, wmacl=12, wparm=13

  ! Multivariate kernel types.
  integer, parameter, public :: ksph=1, kprod=2, kce=3, klm=4, kzeon=5

  ! Predictor styles.
  integer, parameter, public :: stangl=4, stleft=5, strigh=6, stcpar=7

  ! Families.
  integer, parameter, public :: tnul=0, tden=1, trat=2, thaz=3, tgaus=4
  integer, parameter, public :: tlogt=5, tpois=6, tgamm=7, tgeom=8, tcirc=9
  integer, parameter, public :: trobt=10, trbin=11, tweib=12, tcauc=13, tprob=14

  ! Integration types.
  integer, parameter, public :: invld=0, idefa=1, imult=2, iprod=3, imlin=4
  integer, parameter, public :: ihazd=5, isphr=6, imont=7

  ! Prediction targets.
  integer, parameter, public :: pcoef=1, pt0=2, pnlx=3, pband=4, pdegr=5
  integer, parameter, public :: plik=6, prdf=7, pvari=8

  ! Residual types.
  integer, parameter, public :: rdev=1, rpear=2, rraw=3, rldot=4
  integer, parameter, public :: rdev2=5, rlddt=6, rfit=7, rmean=8

  ! Status codes.
  integer, parameter, public :: lf_ok=0, lf_oob=2, lf_pf=3, lf_ncon=4
  integer, parameter, public :: lf_nopt=6, lf_infa=7, lf_demp=10, lf_xoor=11
  integer, parameter, public :: lf_dnop=12, lf_fprob=80, lf_badp=81
  integer, parameter, public :: lf_lnk=82, lf_fam=83, lf_err=99
end module locfit_constants
