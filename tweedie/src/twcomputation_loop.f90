
SUBROUTINE twcomputation_loop(N, p, phi, y, mu, verbose, pdf, funvalue, exitstatus, relerr, Int_Regions)
  ! Calls FORTRAN to compute the integral; set up common parameters
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  USE R_interfaces

  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN)        :: N, verbose, pdf
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: p
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: phi(N), y(N), mu(N)
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: funvalue(N)
  INTEGER(C_INT), INTENT(OUT)       :: exitstatus(N)
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: relerr(N)
  INTEGER(C_INT), INTENT(OUT)       :: Int_Regions(N)
  
  INTEGER(C_INT)        :: i, Int_RegionsTMP, istat
  REAL(KIND=C_DOUBLE)   :: funvalueTMP

  INTERFACE
  
    SUBROUTINE TweedieIntegration(i, funvalueI, exitstatus, relerr, Int_Regions)
      ! Computes the integral in the PDF or CDF expression

      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      USE tweedie_params_mod
      
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN)                :: i
      REAL(KIND=C_DOUBLE), INTENT(OUT)          :: funvalueI
      REAL(KIND=C_DOUBLE), INTENT(OUT)          :: relerr
      INTEGER(C_INT), INTENT(OUT)               :: Int_Regions
      INTEGER(C_INT), INTENT(OUT)               :: exitstatus
    END SUBROUTINE TweedieIntegration

  END INTERFACE


  ! Initialization
  exitstatus = 0_C_INT
  relerr = 0.0_C_DOUBLE
  funvalueTMP = 0.0_C_DOUBLE
  funvalue = 0.0_C_DOUBLE
  Int_Regions = 0_C_INT
  Int_RegionsTMP = 0_C_INT
  istat = 0_C_INT
  
  IF (.NOT. ALLOCATED(Cy)) THEN
      ALLOCATE(Cy(N), Cmu(N), Cphi(N), STAT=istat)
      IF (istat .NE. 0) THEN
        ! CALL INTPR("Allocation failed!", 0)
        RETURN
      END IF
  ELSE IF (SIZE(Cy) .NE. N) THEN
      DEALLOCATE(Cy, Cmu, Cphi)
      ALLOCATE(Cy(N), Cmu(N), Cphi(N), STAT=istat)
      IF (istat .NE. 0) THEN
        ! CALL INTPR("Allocation failed!", 0)
        RETURN
      END IF
  END IF
  
  Cp = p
  Cy = y
  Cmu = mu
  Cphi = phi
  CN = N

  IF (pdf .EQ. 0) THEN 
    Cpdf = .FALSE.        ! Computing the CDF
  ELSE
    Cpdf = .TRUE.         ! Computing the PDF
  END IF
  
  IF (verbose .EQ. 1) THEN
    Cverbose = .TRUE.     ! Verbose feedback
  ELSE
    Cverbose = .FALSE.    ! Minimal feedback
  END IF


  ! Determine case: pSmall = TRUE means 1 < p < 2
  CpSmall = .FALSE.
  IF ( (p .GT. 1.0_C_DOUBLE) .AND. (p .LT. 2.0_C_DOUBLE) ) CpSmall = .TRUE.


  ! Loop over N values
  DO i = 1, N
    CALL TweedieIntegration(i, funvalueTMP, exitstatus(i), relerr(i), Int_RegionsTMP)
    funvalue(i) = funvalueTMP
    Int_Regions(i) = Int_RegionsTMP

  END DO

END SUBROUTINE twcomputation_loop
