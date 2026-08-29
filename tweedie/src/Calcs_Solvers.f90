MODULE Calcs_Solvers

  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE R_interfaces

  IMPLICIT NONE
  PRIVATE
  
  ! Make the public subroutines available to other files using the module
  PUBLIC :: rtsafe, funcd_signature
  
  ! --- Interface Template ---
  ! This defines the signature that rtnewton and rtsafe rely on.
  INTERFACE
    SUBROUTINE funcd_signature(x, f, df, error)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
      
      IMPLICIT NONE
      
      REAL(KIND=C_DOUBLE), INTENT(IN)     :: x
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: f, df
      LOGICAL(C_BOOL), INTENT(INOUT)      :: error
    END SUBROUTINE funcd_signature
  END INTERFACE

CONTAINS

  
  SUBROUTINE rtsafe(funcd, x1, x2, xacc, root, errorHere)
    ! Adapted from NUMERICAL RECIPES Sect. 9.4.
    ! Uses a combination of Newton-Raphson and bisection, find the root of a function bracketed
    ! between x1 and x2. The root is refined until its accuracy is known within plus/minus xacc. 
    !
    ! funcd is a user-supplied subroutine which returns
    ! both the function value and the first derivative of the function.
    
    USE tweedie_params_mod, ONLY: Cverbose
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
    
    IMPLICIT NONE
    
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: x1, x2, xacc
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: root
    LOGICAL(C_BOOL), INTENT(OUT)      :: errorHere
    
    INTEGER(C_INT)                    :: MAXIT, j
    REAL(KIND=C_DOUBLE)               :: df, dx, dxold, f, fh, fl, temp, xh, xl
    
    PARAMETER (MAXIT=100)   ! Maximum allowed number of iterations.
    
    PROCEDURE(funcd_signature) :: funcd
    
    errorHere = .FALSE.
    
    CALL funcd(x1, fl, df, errorHere)
    IF (errorHere) THEN
      IF (Cverbose) CALL DBLEPR("ERROR: In rtsafe (x1)", -1, xl, 1)
      RETURN
    END IF

    CALL funcd(x2, fh, df, errorHere)
    IF (errorHere) THEN
      IF (Cverbose) CALL DBLEPR("ERROR: In rtsafe (x2)", -1, xl, 1)
      RETURN
    END IF

    IF ( (fl .GT. 0.0_C_DOUBLE) .AND. (fh .GT. 0.0_C_DOUBLE) &  
           .OR.                                                &
           (fl .LT. 0.0_C_DOUBLE) .AND. (fh .LT. 0.0_C_DOUBLE) ) THEN
      IF (Cverbose) CALL DBLEPR("ERROR (rtsafe): Root not bracketed: ", -1, x2, 1)
      errorHere = .TRUE.
      RETURN
    END IF

    IF (fl .EQ. 0.0_C_DOUBLE) THEN
      root = x1
      RETURN
    ELSE IF (fh .EQ. 0.0_C_DOUBLE) THEN
      root = x2
      RETURN
    ELSE IF (fl .LT. 0.0_C_DOUBLE) THEN 
    ! Orient the search so that f(xl) < 0.
      xl = x1
      xh = x2
    ELSE
      xh = x1
      xl = x2
    END IF

    root = 0.5_C_DOUBLE * (x1 + x2)   ! Initialize the guess for root,
    dxold = DABS(x2 - x1)             ! the “stepsize before last,”
    dx = dxold                        ! and the last step.
    CALL funcd(root, f, df, errorHere)
    IF (errorHere) THEN
      IF (Cverbose) CALL DBLEPR("ERROR (rtsafe): Computing: ", -1, root, 1)
      RETURN
    END IF

    DO j = 1, MAXIT ! Loop over allowed iterations.
      IF ( ( (root - xh) * df - f) * ( (root - xl) * df - f) .GT. 0.0_C_DOUBLE    & 
             .OR.                                                                   &
             DABS(2.0_C_DOUBLE * f) .GT. DABS(dxold * df) ) THEN
        ! Bisect if Newton out of range,
        ! or not decreasing fast enough.
        dxold = dx
        dx = 0.5_C_DOUBLE *(xh - xl)
        root = xl + dx
        
        IF (xl .EQ. root) RETURN  ! Change in root is negligible.
      ELSE  ! Newton step acceptable; take it.
        dxold = dx
        dx = f/df
        temp = root
        root = root - dx
        IF (temp .EQ. root) RETURN
      END IF
      
      IF (DABS(dx) .LT. xacc) RETURN  ! Convergence criterion.
      CALL funcd(root, f, df, errorHere)      ! The one new function evaluation per iteration.
      IF (errorHere) THEN
        IF (Cverbose) CALL DBLEPR("ERROR (rtsafe): Computing: ", -1, root, 1)
        RETURN
      END IF

      IF (f .LT. 0.0_C_DOUBLE) THEN   ! Maintain the bracket on the root.
        xl = root
      ELSE
        xh = root
      END IF
    END DO
    
    ! If we get here, things have not converged
    errorHere = .TRUE.
    IF (j .GE. MAXIT) THEN
      IF (Cverbose) THEN
        CALL INTPR( "ERROR (rtsafe): failed to convergence after iterations:", -1, MAXIT, 1)
        CALL DBLEPR("  Function value is:", -1, f, 1)
        CALL DBLEPR("  at t:", -1, root, 1)
        CALL DBLEPR("  and slope is:", -1, df, 1)
      END IF
    END IF
    
    RETURN
    
  END SUBROUTINE rtsafe

END MODULE Calcs_Solvers

