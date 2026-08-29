
MODULE Calcs_K

  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE R_interfaces

  IMPLICIT NONE
  
  ! --- Shared Module Variables ---
  ! These are "private" to the module so they don't leak elsewhere, 
  ! but "global" to every subroutine inside.

  PUBLIC :: findKmax, findKmaxSP, improveKZeroBounds, findExactZeros, advanceM
  
CONTAINS

  
  SUBROUTINE findKmax(i, kmax, tmax, mmax, mfirst, left_Of_Max, error)
    ! Finds the value of Kmax, Tmax, and Mmax.
    ! Also return the first value of m (mfirst) amd whether this is to the left of the max (left_Of_Max).
  
    USE tweedie_params_mod
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
    USE Calcs_Imag
    USE Calcs_Solvers
    
    IMPLICIT NONE
    
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: kmax, tmax
    INTEGER(C_INT), INTENT(OUT)         :: mmax, mfirst
    LOGICAL(C_BOOL), INTENT(OUT)        :: left_Of_Max
    INTEGER(C_INT), INTENT(IN)          :: i
    LOGICAL(C_BOOL)                     :: error
  
    REAL(KIND=C_DOUBLE)     :: t_Start_Point, slope_At_Zero
    REAL(KIND=C_DOUBLE)     :: aimrerr, tmaxL, tmaxR, ratio, threshold, t_small
    LOGICAL(C_BOOL)         :: errorHere
    
    ! Grab the relevant scalar values for this iteration:
    current_y    = Cy(i)    ! Access y value for index i
    current_mu   = Cmu(i)   ! Access mu value for index i
    current_phi  = Cphi(i)  ! Access phi value for index i

  
    ! A. If Im k(t) heads down initially: easy: kmax = tmax = mmax = 0
    ! B. Otherwise, we need to find kmax etc by solving Im k'(t) = 0
    !    with a suitable starting point:
    !    1. If tmax appears that it will be very large, find approx start pt
    !    2. Otherwise find start point numerically:
    !       a. If p near 1, approximate starting point by p=1 case analytically
    !       b. Otherwise, use standard starting point
  
    
    ! --- Initialization ---
    aimrerr = 1.0E-09_C_DOUBLE
    threshold = 1.0E5_C_DOUBLE ! Threshold for a "large" tmax
  
  
    ! Find starting points
    errorHere = .FALSE._C_BOOL
    CALL evaluateImkd(0.0_C_DOUBLE, slope_At_Zero, errorHere)
    IF (errorHere) THEN
      error = .TRUE.
      IF (Cverbose) CALL DBLEPR("ERROR: Evaluating Imkd: t", -1, slope_At_Zero, 1)
    END IF

    IF (slope_At_Zero .LE. 0.0_C_DOUBLE) THEN
      ! A. Im k(t) initially heads DOWNWARDS
      !
      !    This includes the case 'IF y >= mu'
      !    Nothing to do; easy-peasy:
      mmax = 0
      mfirst = -1
      kmax = 0.0_C_DOUBLE
      tmax = 0.0_C_DOUBLE
      left_Of_Max = .FALSE.
      
      RETURN
    ELSE
      ! B. CASE: IF slope is initially UPWARDS: trickier, esp. with 1 < p < 2
      !    Need to solve  Im k'(t) = 0.

      IF ( ( current_y ** (1.0_C_DOUBLE - Cp) / current_phi) .GT.  &
           (10.0_C_DOUBLE * threshold ) ) THEN 
        ! B.1. Sometimes kmax and co are MASSIVE, making things slow and difficult, and 
        !      sometimes convergence is difficult also.
        ! In some cases, tmax, kmax and mmax are HUGE; e.g., when we have y=0.001, 
        !     power=6, mu=.01, phi=.01, verbose=TRUE), we can see that 
        ! tmax > 1,000,000, kmax > 80,000 and mmax > 25,000...
        !
        ! When p near 1, we can solve Im k'(t) = 0 when p = 1:
        !    t = acos(y/mu)/phi. 
        ! Much more accurate!
        !
        ! When p near 2, do similar with p = 2; need to check, but I get:
        !    t = 1/(phi*mu) * sqrt( (mu - y)/y )
        ! provided mu - y (which checks out).
        !
        ! This hopefully will flag potentially problematically large kmax/tmax/mmax:
  
        tmax = threshold
        t_Start_Point = tInitialGuess()   ! This computes large-t approx
        CALL evaluateImk(tmax, kmax, errorHere)    ! Now find the corresponding value of kmax
        IF (errorHere) THEN
          error = .TRUE.
          IF (Cverbose) CALL DBLEPR("ERROR: Evaluating Imk: t", -1, tmax, 1)
        END IF

      ELSE
        ! Compute ratio
        ratio = current_y / current_mu
        
        ! 1. Try small-t approximation first (for ratio < 1)
        IF (ratio < 1.0_C_DOUBLE) THEN
            t_small = (current_mu**(1.0_C_DOUBLE - Cp) / current_phi) *    &
                      DSQRT( (2.0_C_DOUBLE/Cp) * (1.0_C_DOUBLE - ratio) )
            ! Clamp tiny values
            IF (t_small < 1.0E-12_C_DOUBLE) t_small = 1.0E-12_C_DOUBLE
            t_Start_Point = t_small
        ELSE
            ! 2. Small-t not applicable → standard start
            IF (Cp < 1.1) THEN
                t_Start_Point = DACOS(current_y/current_mu)/current_phi
            ELSE
                ! Use large-t blended guess
                t_Start_Point = tInitialGuess()   ! This computes large-t approx
                IF (t_Start_Point <= 0.0_C_DOUBLE) t_Start_Point = 1.0E-12_C_DOUBLE
            END IF
        END IF
        ! Now find kmax and tmax using this t_Start_Point
  !      IF (Cpsmall) THEN
          tmaxL = 0.0_C_DOUBLE       ! Since Left bound can be zero
          tmaxR = t_Start_Point * 2.0_C_DOUBLE

          CALL improveKmaxSPBounds(t_Start_Point, tmaxL, tmaxR, error)

          ! Crudely improve the bounds that bracket the starting point for finding Kmax.
          CALL rtsafe(  evaluateImkdZero, tmaxL, tmaxR, aimrerr, tmax, errorHere)
          IF (errorHere) THEN
            error = .TRUE.
            IF (Cverbose) CALL DBLEPR("ERROR: rtsafe: tmax", -1, tmax, 1)
          END IF

      END IF  
  
      ! Find mmax, which depends on whether we are working with the PDF or the CDF.
      ! The PDF uses cos Im k(t) in the integrand; the CDF has sin Im k(t) in the integrand.
      ! Thus, the PDF has integrand zeros at Im k(t) = pi/2 + m pi/y;
      !       the CDF has integrand zeros at Im k(t) =        m pi/y.
  
      ! Find kmax from tmax
      CALL evaluateImk(tmax, kmax, error)
      IF (errorHere) THEN
        error = .TRUE.
        IF (Cverbose) CALL DBLEPR("ERROR: Evaluating Imk: t", -1, tmax, 1)
      END IF

      ! Find mmax from kmax
      IF (Cpdf) THEN
        mmax = FLOOR(2.0_C_DOUBLE * kmax / PI) - 1
      ELSE
        mmax = FLOOR(kmax / PI)
      END IF

      ! Establish the first value of m to use, and whether the first zero is to the left of kmax
      IF (mmax .GT. 0) THEN
        mfirst = 1
        left_Of_Max = .TRUE.
      ELSE
        IF (mmax .EQ. 0 ) THEN
          mfirst = 0
          left_Of_Max = .FALSE.
        ELSE
          ! That is, mmax is LESS THAN 0
          mfirst = -1
          left_Of_Max = .FALSE.
        ENDIF 
      END IF
    END IF

  
    CONTAINS
      
      
      FUNCTION tInitialGuess() RESULT(t0)
      
        IMPLICIT NONE
    
        REAL(KIND=C_DOUBLE)     :: t0
        REAL(KIND=C_DOUBLE)     :: ratio, r, t_small, t_large, c, angle
      

        ! Weighting exponent for interpolation
        r = 4.0E0_C_DOUBLE ! Artibtrary
      
        ! Ratio y/mu
        ratio = current_y / current_mu
      
        ! --- Small-t approximation ---
        IF (ratio < 1.0E0_C_DOUBLE) THEN
           t_small = (current_mu**(1.0E0_C_DOUBLE - Cp) / current_phi) *      &
                      DSQRT( (2.0E0_C_DOUBLE/Cp) * (1.0E0_C_DOUBLE - ratio) )
        ELSE
           ! If y >= mu, no small-t root: set small approx to zero
           t_small = 0.0E0_C_DOUBLE
        ENDIF
      
        ! --- Large-t approximation ---
        angle   = 0.5D0 * PI / (Cp - 1.0E0_C_DOUBLE)
        c       = DCOS(angle)**(Cp - 1.0E0_C_DOUBLE)
        t_large = ( current_mu**(2.0_C_DOUBLE*(1.0_C_DOUBLE - Cp)) /  &
                  ((1.0_C_DOUBLE - Cp) * current_phi) ) *             &
                  ( current_y**(Cp - 1.0_C_DOUBLE) )
      
        ! --- Blended initial guess ---
        t0 = ratio**r * t_small + (1.0E0_C_DOUBLE - ratio**r) * t_large
      
      END FUNCTION tInitialGuess

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      
      SUBROUTINE improveKmaxSPBounds(startTKmax, tmaxL, tmaxR, error)
        ! Crudely improve the bounds that bracket the starting point for finding Kmax.
        ! Sometime a decent starting point is crucial for timely convergence.
        USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
      
        IMPLICIT NONE
        REAL(KIND=C_DOUBLE), INTENT(IN)     :: startTKmax
        REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: tmaxL, tmaxR
      
        ! --- Local Variables ---
        REAL(KIND=C_DOUBLE)     :: boundL, boundR, slopeL, slopeR
        REAL(KIND=C_DOUBLE)     :: oldBoundL, oldBoundR
        REAL(KIND=C_DOUBLE)     :: Imk_value
        INTEGER(C_INT)          :: max_Search, search_Its
        LOGICAL(C_BOOL)         :: keep_Searching, error, errorHere
        

        max_Search = 10
        errorHere = .FALSE._C_BOOL
        
        !!!!! LOWER BOUND
        !   - If slope at SP is *positive* (which it should be), only need to creep to the right
        boundL = tmaxL
        boundR = tmaxR
  
        CALL evaluateImkd(boundL, slopeL, errorHere)
        IF (errorHere) THEN
          error = .TRUE.
          IF (Cverbose) CALL DBLEPR("ERROR: integrand zero =", -1, boundL, 1)
          RETURN
        END IF

        CALL evaluateImk(boundL, Imk_value, errorHere)
        IF (errorHere) THEN
          error = .TRUE.
          IF (Cverbose) CALL DBLEPR("ERROR: derivative zero =", -1, boundL, 1)
          RETURN
        END IF
        
        ! Slope at starting point should always be positive if we are here:
        ! if the slope is negative, Im k(t) heads down and m = -1, -2, ...
  
        ! - Now creep to the right from boundL (refine the lower bound)
        keep_Searching = .TRUE.
        DO WHILE (keep_Searching)
          oldBoundL = boundL
          boundL = (boundL + 1.0E-2_C_DOUBLE) * 1.250E0_C_DOUBLE
          CALL evaluateImkd(boundL, slopeL, error)
          CALL evaluateImk(boundL, Imk_value, error)
          IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, boundL, 1)
  
          IF ( (slopeL .LT. 0.0E0_C_DOUBLE ) .OR.   &
               (Imk_value .LT. 0.0_C_DOUBLE) ) THEN
            ! - Gone too far, so keep previous bound
            keep_Searching = .FALSE.
            boundL = oldBoundL
          END IF
        END DO
        tmaxL = boundL
  
  
        !!!!! UPPER BOUND
        
        ! If the improved lower bouond is LARGER than the original upper bound... FIX! 
        IF (tmaxL .GT. boundR) boundR = tmaxL * 2.0_C_DOUBLE 
        
        ! Find the slope at boundR the starting point (SP)
        CALL evaluateImkd(boundR, slopeR, errorHere)
        IF (errorHere) THEN
          error = .TRUE.
          IF (Cverbose) CALL DBLEPR("ERROR: integrand zero =", -1, boundR, 1)
          RETURN          
        END IF
        
        ! Find the value if kmax at the starting point (SP)
        CALL evaluateImk(boundR, Imk_value, errorHere)
        IF (errorHere) THEN
          error = .TRUE.
          IF (Cverbose) CALL DBLEPR("ERROR: integrand zero =", -1, boundR, 1)
          RETURN          
        END IF
        
  
        ! A valid upper bound:
        ! - must have a negative SLOPE (must be negative)i.e., heading down), *AND*
        ! - must have a positive vale of kmax (or we may have a local maximum only) 
      
        search_Its = 0
        ! - If slope at SP is *positive*, and kmax positive, need bold steps to the right
        IF ( (slopeR .GT. 0.0E0_C_DOUBLE) .AND.  &
             (Imk_value .GT. 0.0_C_DOUBLE) ) THEN
          boundR = startTKmax
          keep_Searching = .TRUE.
          
          DO WHILE (keep_Searching)
            ! - If slope at SP is positive, take bold steps right to find upper bound
            search_Its = search_Its + 1
            boundR = (boundR + 0.1_C_DOUBLE) * 2.0E0_C_DOUBLE
      
            CALL evaluateImkd(boundR, slopeR, errorHere)
            IF (errorHere) THEN
              error = .TRUE.
              IF (Cverbose) CALL DBLEPR("ERROR: integrand zero =", -1, boundR, 1)
              RETURN          
            END IF

            CALL evaluateImk(boundR, Imk_value, errorHere)
            IF (errorHere) THEN
              error = .TRUE.
              IF (Cverbose) CALL DBLEPR("ERROR: integrand zero =", -1, boundR, 1)
              RETURN          
            END IF

            IF ( (slopeR .LT. 0.0E0_C_DOUBLE ) .AND.  &
                 (Imk_value .GT. 0.0_C_DOUBLE) ) THEN
              ! - Found an upper bound where the slope is negative, and kmax is positive
              keep_Searching = .FALSE.
            END IF 
            IF ( search_Its .GT. max_Search) keep_Searching = .FALSE.
          END DO
        END IF
  
        ! - Now creep to the left from boundR (refine the upper bound)
        keep_Searching = .TRUE.
        DO WHILE (keep_Searching)
          oldBoundR = boundR
          boundR = boundR * 0.90E0_C_DOUBLE
          CALL evaluateImkd(boundR, slopeR, errorHere)
          IF (errorHere) THEN
            error = .TRUE.
            IF (Cverbose) CALL DBLEPR("ERROR: integrand zero =", -1, boundR, 1)
            RETURN          
          END IF

          CALL evaluateImk(boundR, Imk_value, errorHere)
          IF (errorHere) THEN
            error = .TRUE.
            IF (Cverbose) CALL DBLEPR("ERROR: integrand zero =", -1, boundR, 1)
            RETURN          
          END IF


          IF ( (slopeR .GT. 0.0E0_C_DOUBLE) .AND.   &
               (Imk_value .GT. 0.0_C_DOUBLE) ) THEN
            ! - Gone too far, so keep previous bound
            keep_Searching = .FALSE.
            boundR = oldBoundR
          END IF
        END DO
        tmaxR = boundR
  
      END SUBROUTINE improveKmaxSPBounds
  
  
  END SUBROUTINE findKmax
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    
  REAL(KIND=C_DOUBLE) FUNCTION findKmaxSP() 
    ! Find the starting point for finding Kmax
    
    USE tweedie_params_mod
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
    USE Calcs_Imag
  
    IMPLICIT NONE
    
    REAL(KIND=C_DOUBLE)           :: tsmall, tlarge, abs1mp
    REAL(KIND=C_DOUBLE)           :: omegaInf, slope
    LOGICAL(C_BOOL)               :: error

    ! Initialize
    abs1mp = ABS(1.0_C_DOUBLE - Cp)
    findKmaxSP = 0.0_C_DOUBLE
    
    ! We find a small-t approx, a large-t approx, and a combined approx.
    ! The SP should be the FIRST of these whose slope is NEGATIVE
  
    IF (CpSmall .AND. (current_y < current_mu) ) THEN
      ! CASE: 1 < p < 2, AND y < mu
  
      ! Try small-t approximation
      tsmall = DSQRT(2.0_C_DOUBLE * (current_mu - current_y)/current_mu) * &
               current_mu**(1.0_C_DOUBLE - Cp) / current_phi
      findKmaxSP = tsmall
    ELSE
      ! This should be the CASE: p > 2, AND y < mu
      omegaInf = (PI / 2.0_C_DOUBLE) * &
                 (1.0_C_DOUBLE - Cp)/(2.0_C_DOUBLE*Cp - 1.0_C_DOUBLE)
      tsmall = current_mu**(1.0_C_DOUBLE - Cp) / ( (1.0_C_DOUBLE - Cp)) * &
               DTAN(omegaInf)
      
      CALL evaluateImkd(tsmall, slope, error)
      IF (error) RETURN
    
      IF (slope .LE. 0.0_C_DOUBLE) THEN
        findKmaxSP = tsmall
        RETURN
      END IF
    END IF
    
    ! Large-t approximation
    tlarge = (current_mu / current_y)**(Cp - 1.0_C_DOUBLE) * &
              current_mu**(1.0_C_DOUBLE - Cp) / (current_phi * abs1mp)
  
    CALL evaluateImkd(tlarge, slope, error)
    IF (error) RETURN
  
    IF (slope .LE. 0.0_C_DOUBLE) THEN
      findKmaxSP = tlarge
      RETURN
    END IF
  
    ! Sum of small + large contributions
    findKmaxSP = tsmall + tlarge
    
  END FUNCTION findKmaxSP
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  
  SUBROUTINE improveKZeroBounds(m, left_Of_Max, zeroMid, zeroL, zeroR, error)
    ! Improve the bounds that bracket the zero of Im k(t).
    ! A decent starting point is sometimes crucial to timely convergence.
    
    USE tweedie_params_mod
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
    USE Calcs_Imag, ONLY: evaluateImkM
    
    IMPLICIT NONE
    
    REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: zeroMid
    INTEGER(C_INT), INTENT(IN)          :: m
    REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: zeroL, zeroR
    LOGICAL(C_BOOL), INTENT(IN)         :: left_Of_Max
    LOGICAL(C_BOOL), INTENT(INOUT)      :: error
  
    REAL(KIND=C_DOUBLE)     :: valueL, valueR, multiplier
    REAL(KIND=C_DOUBLE)     :: df, valueMid
    INTEGER(C_INT)          :: maxSearch,its
    LOGICAL(C_BOOL)         :: stopIterating
  
  

    ! Initialisation
    maxSearch = 10  ! Don't spend too long, so set limit
    valueL = 0.0_C_DOUBLE
    valueR = 0.0_C_DOUBLE
    df = 0.0_C_DOUBLE
    valueMid = 0.95_C_DOUBLE
    
    ! Set multipier: this adjust the sign depending on whether we are
    ! left of the max (so left bound is negative) or to the right of
    ! the max (so left bound is positive)
    IF (left_Of_Max) THEN
      multiplier = -1.0E0_C_DOUBLE
    ELSE
      multiplier = 1.0E0_C_DOUBLE
    END IF
  
  
    ! TWO ROLES:
    ! a) if the bounds actually do bound the zero (i.e, bound give opposite signs for the function
    !    for which zeros are sought): improve if we can, Smetimes, a good start pt is necessary.
    ! b) if the bounds actually do NOT bound the zero: find bounds that do!
  
    ! FIND the function value of the starting point (SP)
    zeroMid = (zeroL + zeroR) / 2.0_C_DOUBLE
    CALL evaluateImkM(zeroMid, valueMid, df, m, error)
    CALL evaluateImkM(zeroL,   valueL,   df, m, error)
    CALL evaluateImkM(zeroR,   valueR,   df, m, error)

    ! CHECK IF BOUNDS REALLY DO BOUND THE ZERO:
    IF ( (valueL * valueR) .GE. 0.0_C_DOUBLE) THEN
      ! Bounds DO NOT trap the zero
  
      ! The solution depends on what side of the max we are.
      ! If to the LEFT of the max of Im k(t), zeroL should give a -ive value; zeroR a +ive value.
      IF ( left_Of_Max) THEN
        its = 0_C_INT
        stopIterating = .FALSE.
        DO WHILE ( .NOT.(stopIterating) )  
          its = its + 1_C_INT
          ! We are on the LEFT of the maximum of Im k(t), but the L bound gives a +ive value.
          ! So we need to go LEFT a little.
          zeroL = (zeroL - 0.1_C_DOUBLE) * 0.95_C_DOUBLE
          CALL evaluateImkM(zeroL, valueL, df, m, error)
          
          IF (its .GE. 100_C_INT) stopIterating = .TRUE.
          IF( valueL .LE. 0.0_C_DOUBLE) stopIterating = .TRUE.
        END DO
  
  
  
        ! The solution depends on what side of the max we are.
        ! If to the LEFT of the max of Im k(t), zeroL should give a -ive value; zeroR a +ive value.
        its = 0_C_INT
        stopIterating = .FALSE.
        DO WHILE ( valueR .LT. 0.0_C_DOUBLE) 
          ! We are on the LEFT of the maximum of Im k(t), but the R bound gives a -ive value.
          ! So we need to go RIGHT a little.
          zeroR = (zeroR + 0.1_C_DOUBLE) * 1.05_C_DOUBLE
          CALL evaluateImkM(zeroR, valueR, df, m, error)

          IF (its .GE. 100_C_INT) stopIterating = .TRUE.
          IF( valueR .GE. 0.0_C_DOUBLE) stopIterating = .TRUE.
        END DO
        
      END IF

      ! The solution depends on what side of the max we are.
      ! If to the RIGHT of the max of Im k(t), zeroL should give a +ive value; zeroR a -ive value.
      IF ( .NOT.(left_Of_Max) ) THEN
        DO WHILE ( valueL .LT. 0.0_C_DOUBLE) 
          ! We are on the RIGHT of the maximum of Im k(t), but the L bound gives a -ive value.
          ! So we need to go LEFT a little.
          zeroL = (zeroL - 0.1_C_DOUBLE) * 0.95_C_DOUBLE
          CALL evaluateImkM(zeroL, valueL, df, m, error)
        END DO

        ! The solution depends on what side of the max we are.
        ! If to the RIGHT of the max of Im k(t), zeroL should give a -ive value; zeroR a +ive value.
        DO WHILE ( valueR .GT. 0.0_C_DOUBLE) 
          ! We are on the RIGHT of the maximum of Im k(t), but the R bound gives a +ive value.
          ! So we need to go RIGHT a little.
          zeroR = (zeroR + 0.1_C_DOUBLE) * 1.05_C_DOUBLE
          CALL evaluateImkM(zeroR, valueR, df, m, error)
        END DO
      END IF

    ELSE
      ! Bounds DO trap the zero, so improve a little
    
      ! Find a point halfway between bounds.
      ! If the new point has same sign as L/R bound, make that the new L/R bounds.
      IF ( (valueMid * valueL) .GE. 0.0_C_DOUBLE) THEN
        zeroL = zeroMid
      ELSE IF ( (valueMid * valueR) .GE. 0.0_C_DOUBLE) THEN
        zeroR = zeroMid
      END IF
        
      ! And once more ONLY
      zeroMid = (zeroL + zeroR) / 2.0_C_DOUBLE
      CALL evaluateImkM(zeroMid, valueMid, df, m, error)
      CALL evaluateImkM(zeroL, valueL, df, m, error)
      CALL evaluateImkM(zeroR, valueR, df, m, error)

      ! Find a point halfway between bounds.
      ! If the new point has same sign as L/R bound, make that the new L/R bounds.
      IF ( (valueMid * valueL) .GE. 0.0_C_DOUBLE) THEN
        zeroL = zeroMid
      ELSE IF ( (valueMid * valueR) .GE. 0.0_C_DOUBLE) THEN
        zeroR = zeroMid
      END IF
      
      zeroMid = (zeroL + zeroR) / 2.0_C_DOUBLE
  
    END IF
    
    RETURN
  
  END SUBROUTINE improveKZeroBounds
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  SUBROUTINE advanceM(m, mmax, mOld, left_Of_Max)
      ! Determine the next value of m, for solving the zeros of the integrand
      
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
    
      IMPLICIT NONE
      
      INTEGER(C_INT), INTENT(IN)        :: mmax       ! Maximum value m can take, and current index
      INTEGER(C_INT), INTENT(INOUT)     :: m          ! M index (used for calculation and C-binding)
      INTEGER(C_INT), INTENT(OUT)       :: mOld       ! Previous index value
      LOGICAL(C_BOOL), INTENT(INOUT)    :: left_Of_Max  ! True if on the left side of kmax

      ! Local vars
      LOGICAL(C_BOOL)                   :: flip_To_Other_Side
        ! flip_To_Other_Side. is  .TRUE.  if the value of  m  crosses from the 
        ! left of  mmax   to the right of  mmax
    

      mOld = m
      flip_To_Other_Side = .FALSE.
        ! flip_To_Other_Side signals that m has just crossed mmax;
        ! update integration side exactly once

      IF (current_y .GE. current_mu) THEN
        ! Always heading downwards (away from kmax), so easy
        m = m - 1 
      ELSE
        ! We have a maximum (kmax) to consider
        IF (left_Of_Max) THEN
          IF (m == mmax) THEN 
            ! Move to the other side of the maximum
            left_Of_Max = .FALSE.
            flip_To_Other_Side = .TRUE.
            ! NOTE: The value of m does not change for this iteration; 
            ! same value, just on other side of max
          ELSE
            ! Continue towards the maximum
            m = m + 1
            ! mOld is already saved before the IF block
            left_Of_Max = .TRUE. ! Still on the left side
          END IF
        ELSE
          ! When on the RIGHT of the maximum, can always just reduce m by one
          m = m - 1 
          left_Of_Max = .FALSE.
        END IF
      END IF

    END SUBROUTINE advanceM
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    SUBROUTINE findExactZeros(m, tL, tR, tStart, tZero, left_Of_Max, error) 
  ! Find the exact zeros of the integrand
  
    USE tweedie_params_mod
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
    USE Calcs_Solvers
    USE Calcs_Imag
    
    IMPLICIT NONE
  
    INTEGER(C_INT), INTENT(IN)          :: m
    REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: tL, tR, tStart
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: tZero
    LOGICAL(C_BOOL), INTENT(IN)         :: left_Of_Max
    LOGICAL(C_BOOL), INTENT(INOUT)      :: error
  
    REAL(KIND=C_DOUBLE)   :: xacc, fL, fR, dfL, dfR, tstart_update, tMid
    LOGICAL(C_BOOL)       :: errorHere
  
    ! INITIALIZE ALL LOCAL VARIABLES
    errorHere = .FALSE.
    xacc = 1.0E-11_C_DOUBLE
    fL = 0.0_C_DOUBLE
    fR = 0.0_C_DOUBLE
    dfL = 0.0_C_DOUBLE
    dfR = 0.0_C_DOUBLE
    tstart_update = 0.0_C_DOUBLE
    tMid = 0.0_C_DOUBLE
    tZero = 0.0_C_DOUBLE  ! Initialize OUT parameter too
    
    ! Sync the local value of  m  to the shared value.
    m_shared = m

    ! Set the accuracy
    xacc = 1.0E-11_C_DOUBLE
  
    ! Find mmax, which depends on whether we are working with the PDF or the CDF.
    ! The PDF uses cos Im k(t) in the integrand; the CDF has sin Im k(t) in the integrand.
    ! Thus, the PDF has integrand zeros at Im k(t) = pi/2 + m * pi/y;
    !       the CDF has integrand zeros at Im k(t) =        m * pi/y.
    ! Ensure the bounds actually bound the zero
  
    CALL evaluateImkM(tL, fL, dfL, m, errorHere)
    IF (errorHere) THEN
      error = .TRUE.
      RETURN
    END IF

    CALL evaluateImkM(tR, fR, dfR, m, errorHere)
    IF (errorHere) THEN
      error = .TRUE.
      RETURN
    END IF

    IF ( (fL * fR) .GT. 0.0_C_DOUBLE ) THEN
      ! Then bounds do not bound the zero
       tMid = (tL + tR) / 2.0_C_DOUBLE
   
      CALL improveKZeroBounds(m, left_Of_Max, tMid, tL, tR, errorHere)
      IF (errorHere) THEN
        IF (Cverbose) CALL DBLEPR("Bounds do not bracket the zero (findExactZeros)", -1, fR, 1)
        error = .TRUE.
        RETURN
      END IF
    END IF

    ! For robustness, use rtsafe when the  distance between zeros 
    ! is expected to be small (i.e., in the tail).
    IF ( m .LE. -3 ) THEN 
    
      CALL rtsafe(evaluateImkM_wrapper, tL, tR, xacc, tZero, errorHere)
      IF (errorHere) THEN
        error = .TRUE.
        IF (Cverbose) CALL DBLEPR("ERROR: cannot solve", -1, tZero, 1)
        RETURN
      END IF  
    ELSE IF ( (Cpsmall) .AND. (current_y .LT. current_mu) ) THEN
    
      ! When small p and small y, fight harder for good starting bounds
      CALL improveKZeroBounds(m, left_Of_Max, tStart, tL, tR, error)
!      IF (error) THEN
!        IF (Cverbose) CALL DBLEPR("ERROR: cannot solve", -1, tZero, 1)
!        RETURN
!      END IF
      CALL rtsafe(evaluateImkM_wrapper, tL, tR, xacc, tZero, errorHere)

      IF (errorHere) THEN
        error = .TRUE.
        IF (Cverbose) CALL DBLEPR("ERROR: cannot solve", -1, tZero, 1)
        RETURN
      END IF
    ELSE
      ! Default to rtnewton for "easy" cases (e.g., initial zeros)
      tstart_update = (tL + tR) / 2.0_C_DOUBLE
  !    CALL rtnewton(i, evaluateImkM_wrapper, tstart_update, xacc, tZero)
      CALL rtsafe(evaluateImkM_wrapper, tL, tR, xacc, tZero, errorHere)
      
      IF (errorHere) THEN
        error = .TRUE.
        RETURN
      END IF
    END IF
    
  
  END SUBROUTINE findExactZeros
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  ! Define the wrapper subroutine
  SUBROUTINE evaluateImkM_wrapper(x, f, df, error)
    ! Evaluate  Im k(t) - m * pi (CDF) - pi/2  or  Im k(t) - m * pi (CDF)
    
    USE Calcs_Imag, ONLY: evaluateImkM  ! Explicitly tell the wrapper where to find the math
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: x
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: f, df
    LOGICAL(C_BOOL), INTENT(INOUT)    :: error
       
    CALL evaluateImkM(x, f, df, m_shared, error) ! Pass the captured 'm' value

  END SUBROUTINE evaluateImkM_wrapper

  
END MODULE Calcs_K

