
SUBROUTINE TweedieIntegration(i, funvalueI, exitstatus, relerr, count_Integration_Regions) 
  ! Compute the value of the integrals in the Fourier-inversion expressions for the PDF and CDF

  USE Integrands_MOD, ONLY: Integrands
  USE tweedie_params_mod
  USE TweedieIntHelpers
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE rprintf_mod
  USE Calcs_Imag
  USE Calcs_Real
  USE Calcs_K
  USE R_interfaces

  IMPLICIT NONE
  
  INTEGER(C_INT), INTENT(IN)        :: i              ! Observation index
  INTEGER(C_INT), INTENT(OUT)       :: exitstatus     ! Output status
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: funvalueI      ! Computed result
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: relerr         ! Estimate of relative error
  INTEGER(C_INT), INTENT(OUT)       :: count_Integration_Regions ! Num int regions

  ! Local Variables: All local variables defined here
  LOGICAL(C_BOOL)       :: error          ! TRUE if any computational, numerical problems found
  INTEGER(C_INT)        :: mmax, mfirst, mOld, accMax
  INTEGER(C_INT)        :: m, min_Acc_Regions
  LOGICAL(C_BOOL)       :: convergence_Acc
  REAL(KIND=C_DOUBLE)   :: kmax, tmax, aimrerr
  REAL(KIND=C_DOUBLE)   :: epsilon, areaT, West, Wold, Wold2
  REAL(KIND=C_DOUBLE)   :: zeroL, zeroR
  REAL(KIND=C_DOUBLE), ALLOCATABLE   :: Mmatrix(:, :), Nmatrix(:, :), xvec(:), wvec(:)
  REAL(KIND=C_DOUBLE)   :: zeroStartPoint
  LOGICAL(C_BOOL)       :: left_Of_Max, flip_To_Other_Side
  
  INTEGER, PARAMETER :: MAX_ACC = 200
  INTEGER, PARAMETER :: VEC_SIZE = MAX_ACC + 2
  REAL(C_DOUBLE), PARAMETER :: EPS = 1.0E-12_C_DOUBLE

  ! Zone 1: initial region
  REAL(C_DOUBLE)  :: area0

  ! Zone 2: pre-acceleration
  REAL(C_DOUBLE) :: area1, sumA
  INTEGER(C_INT) :: count_PreAcc_Regions
  LOGICAL(C_BOOL):: stop_PreAccelerate, converged_Pre
  REAL(C_DOUBLE) :: leftPreAccZero

  ! Zone 3: acceleration
  REAL(C_DOUBLE) :: areaA, psi
  INTEGER(C_INT) :: count_Acc_Regions
  LOGICAL(C_BOOL):: keep_Accelerating, converged_Accelerating
  REAL(C_DOUBLE) :: leftAccZero


  INTERFACE

    SUBROUTINE GaussQuadrature(i, a, b, area)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      REAL(KIND=C_DOUBLE), INTENT(OUT)      :: area
      REAL(KIND=C_DOUBLE), INTENT(IN)       :: a, b
      INTEGER(C_INT), INTENT(IN)            :: i
    END SUBROUTINE GaussQuadrature

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    SUBROUTINE accelerate(xvec, wvec, nzeros, Mmatrix, NMatrix, West)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
      IMPORT :: VEC_SIZE  ! Bring the parameter into the interface scope
      INTEGER(C_INT), INTENT(IN)          :: nzeros
      REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: xvec(VEC_SIZE), wvec(VEC_SIZE)
      REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: Mmatrix(2, VEC_SIZE), Nmatrix(2, VEC_SIZE)
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: West
    END SUBROUTINE accelerate

  END INTERFACE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  ! Begin assuming no errors, and set to  .TRUE.  if any errors found.
  error = .FALSE.
  
  ! ALLOCATE these arrays onto the HEAP
  ALLOCATE(Mmatrix(2, VEC_SIZE))
  ALLOCATE(Nmatrix(2, VEC_SIZE))
  ALLOCATE(xvec(VEC_SIZE))
  ALLOCATE(wvec(VEC_SIZE))
  


  ! --- Initialization ---
  aimrerr = 1.0E-12_C_DOUBLE
  mOld = 0
  m = 0
  exitstatus = 0
  relerr = 1.0_C_DOUBLE
  convergence_Acc = .FALSE.
  epsilon = 1.0E-12_C_DOUBLE
  mmax = 0_C_INT
  count_Integration_Regions = 0_C_INT   ! Counter for number of integration regions
  zeroStartPoint = 0.0_C_DOUBLE
  flip_To_Other_Side = .FALSE.
  zeroR = 0.0_C_DOUBLE
  zeroL = 0.0_C_DOUBLE
  left_Of_Max = .TRUE.
  area0 = 0.0_C_DOUBLE



  ! Initialise for pre acceleration
  area1 = 0.0_C_DOUBLE
  stop_PreAccelerate = .FALSE.        ! Stop pre-accelerating (and perhaps move to accelerating)
  converged_Pre = .FALSE.             ! If convergence detected during pre-acceleration
  sumA = 0.0_C_DOUBLE                 ! Total pre-acceleration area
  West = 3.0_C_DOUBLE                 ! The current estimate of the tail area
  Wold = 2.0_C_DOUBLE                 ! Wold and Wold2 are the previous two estimates of the tail area,
  Wold2 = 1.0_C_DOUBLE                !   to allow testing for convergence
  count_PreAcc_Regions = 0_C_INT      ! Number of pre-acceleratioin regions


  ! Initialise for acceleration
  count_Acc_Regions = 0_C_INT
  keep_Accelerating = .TRUE.
  Mmatrix = 0.0_C_DOUBLE
  Nmatrix = 0.0_C_DOUBLE
  xvec = 0.0_C_DOUBLE
  wvec = 0.0_C_DOUBLE
  areaA = 0.0_C_DOUBLE
  count_Acc_Regions = 0_C_INT
  accMax = MAX_ACC                      ! Max acceleration regions
  min_Acc_Regions = 3_C_INT             ! Min preacceleration regions


  ! --- Integration initialization ---
  area0 = 0.0_C_DOUBLE ! Initial area



  ! FIND THE VALUES OF  kmax, tmax, mmax
  CALL findKmax(i, kmax, tmax, mmax, mfirst, left_Of_Max, error)
  m = mfirst

  ! INTEGRATION
  ! Three integration zones:
  !
  !   1. The *initial* area, which is not between zeros of Im{k(t)}: area0
  !      area0 is between tL = 0 and tR.
  !
  !   2. The initial area *before* Sidi acceleration is invoked: area1
  !      (For instance, wait until Im{k(t)} is on the downturn.)
  !      So we set tL to be the value of tR is the initial region.
  !      Then, area1 is between tL set above, and the next zero (i.e., value
  !       of t), and repeat until it is time for acceleration.
  !
  !   3. The area thereafter, upon which Sidi acceleration is
  !      applied; the area returned by acceleration is areaA
  !      Acceleration starts with the tL set as the last value of tR
  !      used in pre-acceleration.




  ! ----------------------------------------------------------------------------
  ! --- 1. INTEGRATE FIRST (sometimes non-standard) REGION: area0 ---
  
  ! Find the value of  zeroR  for the initial region (zeroL is always 0.0)
  CALL findInitialZeroR(mfirst, left_Of_Max, tmax, &
                        zeroR, error)
  ! Integrate:
  CALL GaussQuadrature(i, zeroL, zeroR, area0)   ! area0  is the area of the initial region

  ! Update
  CALL updateTM( i, tmax, mmax, left_Of_Max, &
                 m, zeroL, zeroR, error)


  ! ----------------------------------------------------------------------------
  ! --- 2. INTEGRATE: the PRE-ACCELERATION regions: area1 ---

  ! Update to get the next values of zeroL and zeroR, and the value of m
  ! corresponding to zeroR
  !


  ! The next integration region is between zeroL and zeroR 
  ! (and  m  corresponds to this value of zeroR)

  ! Retain the value of  t  where pre-acceleration starts
  leftPreAccZero = zeroL

  DO WHILE ( .NOT.(stop_PreAccelerate ) )
    count_PreAcc_Regions = count_PreAcc_Regions + 1_C_INT    ! Update count

    ! Integrate
    CALL GaussQuadrature(i, zeroL, zeroR, sumA)
    area1 = area1 + sumA

    ! Update (zeroL, zeroR and m)
    CALL updateTM( i, tmax, mmax, left_Of_Max, &
                   m, zeroL, zeroR, error)

    ! Check for convergence
    CALL checkStopPreAcc(tmax, zeroR, stop_PreAccelerate, converged_Pre, error)
    IF (count_Acc_Regions .GT. accMax) THEN
      stop_PreAccelerate = .TRUE.
      converged_Pre      = .FALSE.
    END IF

  END DO




  ! ----------------------------------------------------------------------------
  ! --- 3. INTEGRATE: the ACCELERATION regions: areaA ---

  ! Retain the value of  t  where acceleration starts
  leftAccZero = zeroL
  
  ! Initialise the acceleration
  count_Acc_Regions = 0_C_INT
  xvec(1) = zeroL
  
  IF ( .NOT.(converged_Pre) ) THEN
    converged_Accelerating = .FALSE.

    DO WHILE ( keep_Accelerating )
      count_Acc_Regions = count_Acc_Regions + 1_C_INT
      
      xvec(count_Acc_Regions + 1_C_INT) = zeroR  ! Update xvec

      ! Integrate
      CALL GaussQuadrature(i, zeroL, zeroR, psi)

      ! Prepare for acceleration
      wvec(count_Acc_Regions)     = psi   ! Update past estimates
                                          ! wvec contains the sequence of integration areas, 
                                          ! starting with the first and up to the limit
      Wold2 = Wold
      Wold  = West  

      ! Accelerate (i.e., update West, the best area estimate)
      CALL accelerate(xvec, wvec, count_Acc_Regions, Mmatrix, Nmatrix, West)

      ! Update (zeroL, zeroR and m)
      CALL updateTM( i, tmax, mmax, left_Of_Max, &
                     m, zeroL, zeroR, error)

      ! Check for convergence
      relerr = ( DABS(West - Wold) + DABS(West - Wold2)) / (DABS(West) + epsilon)
      IF ( (count_Acc_Regions .GE. min_Acc_Regions) .AND. &
           (relerr .LT. aimrerr) ) THEN
        keep_Accelerating = .FALSE.
        converged_Accelerating = .TRUE.
        convergence_Acc = .TRUE.
      END IF
      
      IF (count_Acc_Regions .GE. accMax) THEN
        keep_Accelerating = .FALSE.
        converged_Accelerating = .FALSE.
      END IF
! IF (m .LT. -3) STOP

    END DO
    
    areaA = West
  END IF



  ! --- WIND THINGS UP ---
  count_Integration_Regions = 1_C_INT  +              &   ! Initial zone has one integration region
                              count_PreAcc_Regions +  &   ! Pre-acc regions
                              count_Acc_Regions           ! Acc regions
  areaT = area0 + area1 + areaA

  ! We have the value of the integral in the PDF/CDF calculation, so now work out the actual PDF/CDF
  IF ( Cpdf ) THEN
    funvalueI = areaT/PI 
  ELSE
    funvalueI =  0.5_C_DOUBLE - areaT/PI
  END IF
  
  IF (error) exitstatus = 1_C_INT
  
  
  ! Print things if requested
  IF ( Cverbose ) THEN
    ! -------- Preparatory
    CALL DBLEPR("  -            kmax:", -1, kmax, 1 )
    CALL DBLEPR("  -            tmax:", -1, tmax, 1 )
    CALL INTPR( "  -            mmax:", -1, mmax, 1 )
    CALL INTPR( "  - first zero at m:", -1, mfirst, 1 )

    ! -------- Initial zone
    CALL DBLEPR("Initial region area:", -1, area0, 1)
    CALL DBLEPR("      between 0 and:", -1, zeroR, 1)
    CALL INTPR( "      using right m:", -1, m, 1)
    CALL INTPR( " # pre-acc regions: ", -1, count_PreAcc_Regions, 1)

    ! -------- Pre-acceleration zone
    CALL DBLEPR("       Pre-acc AREA:", -1, area1, 1)
    CALL DBLEPR("            between:", -1, leftPreAccZero, 1)
    CALL DBLEPR("                and:", -1, zeroR, 1)
    CALL INTPR( "      using right m:", -1, m,     1)
    CALL INTPR( "     # acc regions: ", -1, count_Acc_Regions, 1)

    ! -------- Acceleration zone
    IF (converged_Pre) THEN
      CALL DBLEPR(" Accelerating not needed; convergence by t =", -1, zeroR, 1)
    ELSE
      CALL DBLEPR(" Accelerating starting after t =", -1, zeroR, 1)
      CALL DBLEPR("         Acc area:", -1, areaA, 1)
      CALL DBLEPR("          between:", -1, leftAccZero, 1)
      CALL DBLEPR("              and:", -1, zeroR, 1)
      CALL INTPR( "         up to m:", -1, m,     1)
    END IF
    
    ! -------- Summary
    CALL DBLEPR("*** Initial area0: ", -1, area0, 1)
    CALL DBLEPR("*** Pre-acc area1: ", -1, area1, 1)
    CALL DBLEPR("***     Acc area!: ", -1, areaA, 1)
    CALL DBLEPR("***         TOTAL: ", -1, areaT, 1)
    CALL INTPR( "   over regions: ", -1, count_Integration_Regions, 1)
    
    ! -------- Results
    CALL DBLEPR("***    Fun. value:", -1, funvalueI, 1)

  END IF


  ! Tidy up: deallocate the arrays
  DEALLOCATE(Mmatrix, Nmatrix, xvec, wvec)
  
END SUBROUTINE TweedieIntegration

