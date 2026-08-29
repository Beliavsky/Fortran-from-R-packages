
MODULE Calcs_Imag

  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE R_interfaces
  
  IMPLICIT NONE
  
  PUBLIC :: evaluateImk, evaluateImkM, evaluateImkd, evaluateImkdd, evaluateImkdZero
  
CONTAINS

  SUBROUTINE evaluateImk(t, Imk, errorHere) 
    ! Evaluate Im k(t)
    
    IMPLICIT NONE
    
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: Imk
    LOGICAL(C_BOOL), INTENT(OUT)        :: errorHere
    
    REAL(KIND=C_DOUBLE)   :: tanArg, omega, front, alpha

  
    ! Initialisation
    errorHere = .FALSE.
    Imk = 0.0_C_DOUBLE
    
    front = current_mu ** (2.0_C_DOUBLE - Cp) / ( current_phi * (2.0_C_DOUBLE - Cp))
    tanArg = (1.0_C_DOUBLE - Cp) * t * current_phi  / (current_mu ** (1.0_C_DOUBLE - Cp) )
    omega = DATAN( tanArg )
  
    IF ((omega .GT. 0.0_C_DOUBLE ) .OR.    &    
        (omega .LT. (-PI/2.0_C_DOUBLE)) ) THEN
      errorHere = .TRUE.
      IF (Cverbose) CALL DBLEPR("ERROR: Imk: omega out of bounds:", -1, omega, 1)
      RETURN
    END IF
    alpha = (2.0_C_DOUBLE - Cp)/(1.0_C_DOUBLE - Cp)
  
    Imk = front *   &
          DSIN(omega * alpha)/(DCOS(omega) ** alpha) - t * current_y
  
    RETURN
  
  END SUBROUTINE evaluateImk
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  SUBROUTINE evaluateImkM(t, f, df, m, error)
    ! Evaluate  Im k(t) - m*pi - pi/2  (PDF) or  Im k(t) - m*pi  (CDF) for 
    ! finding the zeros of the integrand, for the PDf and CDF
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: f, df
    INTEGER(C_INT), INTENT(IN)        :: m
    LOGICAL(C_BOOL)                   :: error

    REAL(KIND=C_DOUBLE)               :: Imk_val
    LOGICAL(C_BOOL)                   :: errorHere
  
    
    ! Initialise
    f = 0.0_C_DOUBLE
    df = 0.0_C_DOUBLE
    errorHere = .FALSE.
    
    CALL evaluateImk(t, Imk_val, errorHere)
    IF (errorHere) THEN
      IF (Cverbose) CALL DBLEPR("ERROR: error evaluating evaluate Im k(t) =", -1, t, 1)
      error = .TRUE.
      RETURN
    END IF
  
    ! The expression depends on whether we are working with the PDF or the CDF.
    ! The PDF uses cos Im k(t) in the integrand; the CDF has sin Im k(t) in the integrand.
    ! Thus, the PDF has integrand zeros at Im k(t) = pi/2 + m pi/y;
    !       the CDF has integrand zeros at Im k(t) =        m pi/y.
    IF (Cpdf) THEN
      f = Imk_val - REAL(m, KIND=C_DOUBLE) * PI - PI/2.0_C_DOUBLE
    ELSE
      f = Imk_val - REAL(m, KIND=C_DOUBLE) * PI
    END IF
    CALL evaluateImkd(t, df, errorHere)
    IF (errorHere) THEN
      error = .TRUE.
      IF (Cverbose) CALL DBLEPR("ERROR: Im k(t) - m:", -1, f, 1)
      RETURN
    END IF
    
  END SUBROUTINE evaluateImkM

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE evaluateImkd(t, Imkd, errorHere)
    ! Evaluate Im k'(t)
    
    IMPLICIT NONE
    
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Imkd  ! The result of the calculation
    LOGICAL(C_BOOL), INTENT(OUT)      :: errorHere
    REAL(KIND=C_DOUBLE) :: omega, pindex

  
    ! Initlaise
    Imkd = 0.0_C_DOUBLE
    errorHere = .FALSE.
    
    pindex = 1.0_C_DOUBLE / (1.0_C_DOUBLE - Cp)
    omega = DATAN( ( (1.0_C_DOUBLE - Cp) * t * current_phi) / (current_mu ** (1.0_C_DOUBLE - Cp) ) )

    IF ((omega .GT. 0.0_C_DOUBLE ) .OR.    &    
        (omega .LT. (-PI/2.0_C_DOUBLE)) ) THEN
      IF (Cverbose) CALL DBLEPR("ERROR: Imkd: omega out of bounds:", -1, omega, 1)
      errorHere = .TRUE.
      
      RETURN
    END IF
    
    Imkd = current_mu * (DCOS(omega * pindex) / (DCOS(omega) ** pindex)) - current_y
  
  END SUBROUTINE evaluateImkd


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE evaluateImkdd(t, Imkdd, errorHere)
    ! Evaluate Im k''(t)
    
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Imkdd
    LOGICAL(C_BOOL), INTENT(INOUT)    :: errorHere
    
    REAL(KIND=C_DOUBLE)    :: front, omega, pindex

    
    ! Initialise
    Imkdd = 0.0_C_DOUBLE
    errorHere = .FALSE.
    
    pindex = Cp / (1.0_C_DOUBLE - Cp)
    front = -current_phi * current_mu ** (Cp/(1.0_C_DOUBLE - Cp))
    omega = DATAN( ( (1.0_C_DOUBLE - Cp) * t * current_phi) / (current_mu ** (1.0_C_DOUBLE - Cp) ) )
    
    IF ((omega .GT. 0.0_C_DOUBLE ) .OR.    &    
        (omega .LT. (-PI/2.0_C_DOUBLE)) ) THEN
      IF (Cverbose) CALL DBLEPR("ERROR: Imkdd: omega out of bounds:", -1, omega, 1)
      errorHere = .TRUE.
      RETURN
    END IF

  
    Imkdd = front * (DSIN(omega * pindex) / (DCOS(omega) ** pindex) )


  END SUBROUTINE evaluateImkdd
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE evaluateImkdZero(t, f, df, error) 
  ! Evaluate Im k'(t)  and  Im k''(t)  for solving for Kmax (i.e., Im k'(t) = 0)
    
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: f, df
    LOGICAL(C_BOOL), INTENT(INOUT)      :: error
  
    REAL(KIND=C_DOUBLE)   :: Imkd, Imkdd
    LOGICAL(C_BOOL)       :: errorHere
  
    ! Initialise
    f = 0.0_C_DOUBLE
    df = 0.0_C_DOUBLE
    
    CALL evaluateImkd( t, Imkd, errorHere)
    IF (errorHere) THEN
      error = .TRUE.
    END IF
    
    CALL evaluateImkdd(t, Imkdd, error)
    IF (errorHere) THEN
      error = .TRUE.
    END IF
    
    f  = Imkd
    df = Imkdd
  END SUBROUTINE evaluateImkdZero



END MODULE Calcs_Imag

