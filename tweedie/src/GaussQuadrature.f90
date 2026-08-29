SUBROUTINE GaussQuadrature(i, a, b, integral_result) 
  ! Perform 512-pt Gaussian-Legendre quadrature to evaluate integrals
  
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  USE Integrands_MOD, ONLY: Integrands
  USE gaussian_data_mod, ONLY: absc, wts
  
  IMPLICIT NONE
  
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: a, b      ! Integration limits
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: integral_result
  INTEGER(C_INT), INTENT(IN)        :: i

  INTEGER(C_INT)                                 :: j, npoints
  REAL(KIND=C_DOUBLE)                            :: xu, xl, fl, fu
  
  !!! NOTE: Guassian abscicca and weights in gaussian_data_mod
  
  ! Set up initial parameters
  integral_result = 0.0E0_C_DOUBLE
  npoints = 256 ! For 512-pt quadrature: symmetry
  
  ! Compute
  DO j = 1, npoints
    ! Adjust abscissae
    xl = ( b - a ) / 2.0E0_C_DOUBLE * absc(j) + ( b + a ) / 2.0E0_C_DOUBLE
    xu = ( a - b ) / 2.0E0_C_DOUBLE * absc(j) + ( b + a ) / 2.0E0_C_DOUBLE

    ! Evaluate
    fl = Integrands(i, xl)
    fu = Integrands(i, xu)
    integral_result = integral_result + wts(j) * (fl + fu)
  END DO
  
  integral_result = integral_result * (b - a) / 2.0E0_C_DOUBLE

END SUBROUTINE GaussQuadrature
