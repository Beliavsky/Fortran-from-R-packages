
MODULE tweedie_params_mod
  ! Save commonly-used parameters, for access throughout functions and subroutines
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE , C_BOOL
  
  IMPLICIT NONE
  SAVE
  
  ! Global Parameters
  REAL(KIND=C_DOUBLE), ALLOCATABLE  :: Cmu(:), Cphi(:), Cy(:) 
  REAL(KIND=C_DOUBLE)               :: Cp
  LOGICAL(C_BOOL)                   :: CpSmall, Cverbose, Cpdf
  INTEGER(C_INT)                    :: CN
  

  REAL(KIND=C_DOUBLE) :: current_y, current_mu, current_phi
  INTEGER(C_INT)      :: m_shared

  REAL(KIND=C_DOUBLE), PARAMETER :: PI = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)

END MODULE tweedie_params_mod

