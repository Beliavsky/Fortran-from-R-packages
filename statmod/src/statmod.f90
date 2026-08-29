! Umbrella module for the modern Fortran translation of statmod.
module statmod
use statmod_special
use statmod_invgauss
use statmod_gaussquad
use statmod_expected_deviance
use statmod_utils
use statmod_elda
use statmod_models
use statmod_qres
use statmod_glm_misc
implicit none
public
end module statmod
