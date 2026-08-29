SUBROUTINE accelerate(xvec, wvec, nzeros, Mmatrix, Nmatrix, West)
  ! Accelerate the convergence of the infinite integral., the tail area
  
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  
  IMPLICIT NONE
  
  INTEGER, PARAMETER                  :: VEC_SIZE = 202
  INTEGER(C_INT), INTENT(IN)          :: nzeros
  REAL(KIND=C_DOUBLE), INTENT(IN)     :: xvec(VEC_SIZE), wvec(VEC_SIZE)
  REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: Mmatrix(2, VEC_SIZE), Nmatrix(2, VEC_SIZE)
  REAL(KIND=C_DOUBLE), INTENT(OUT)    :: West
  INTEGER(C_INT)                  :: p, q, l_nzeros, maxSize
  REAL(KIND=C_DOUBLE)             :: denom, psi_new, FF_current, sumw
  REAL(KIND=C_DOUBLE)             :: tinyDenom, scale_denom
  REAL(KIND=C_DOUBLE)             :: inv_x_i, inv_x_l
  REAL(KIND=C_DOUBLE)             :: s, maxinv, invx
  REAL(KIND=C_DOUBLE), PARAMETER  :: SAFETY_SCALE = 1.0D-12
  REAL(KIND=C_DOUBLE), PARAMETER  :: HUGE_LIMIT   = 1.0D300
  REAL(KIND=C_DOUBLE), ALLOCATABLE :: xscaled(:)

  maxSize  = 200
  l_nzeros = MIN(nzeros, maxSize)
  ALLOCATE(xscaled(l_nzeros))

  maxinv = 0.0E0_C_DOUBLE
  DO p = 1, l_nzeros
     invx = DABS(1.0E0_C_DOUBLE / xvec(p))
     IF (invx .GT. maxinv) maxinv = invx
  END DO
  s = MAX(1.0E0_C_DOUBLE, 1.0E0_C_DOUBLE / maxinv)
  DO p = 1, l_nzeros
     xscaled(p) = xvec(p) / s
  END DO

  psi_new = wvec(l_nzeros)
  sumw = 0.0_C_DOUBLE
  DO p = 1, l_nzeros
     sumw = sumw + wvec(p)
  END DO
  FF_current = sumw

  IF (DABS(psi_new) .LT. 1.0E-31_C_DOUBLE) THEN
     West = FF_current
     DO q = 1, l_nzeros
        Mmatrix(1, q) = Mmatrix(2, q)
        Nmatrix(1, q) = Nmatrix(2, q)
     END DO
     DEALLOCATE(xscaled)
     RETURN
  END IF

  Mmatrix(2, 1) = FF_current / psi_new
  Nmatrix(2, 1) = 1.0E0_C_DOUBLE / psi_new
  tinyDenom = SAFETY_SCALE
  inv_x_l = 1.0E0_C_DOUBLE / xscaled(l_nzeros)

  DO p = 2, l_nzeros
     inv_x_i = 1.0E0_C_DOUBLE / xscaled(l_nzeros + 1 - p)
     denom   = inv_x_i - inv_x_l
     scale_denom = MAX( DABS(inv_x_i), DABS(inv_x_l), 1.0E0_C_DOUBLE )
     tinyDenom   = SAFETY_SCALE * scale_denom
     IF ( (.NOT.(denom .EQ. denom)) .OR. (DABS(denom) .LT. tinyDenom) ) THEN
        West = FF_current
        DO q = 1, l_nzeros
           Mmatrix(1, q) = Mmatrix(2, q)
           Nmatrix(1, q) = Nmatrix(2, q)
        END DO
        DEALLOCATE(xscaled)
        RETURN
     END IF

     Mmatrix(2, p) = (Mmatrix(1, p - 1) - Mmatrix(2, p - 1)) / denom
     Nmatrix(2, p) = (Nmatrix(1, p - 1) - Nmatrix(2, p - 1)) / denom

     IF (.NOT.(Mmatrix(2, p) .EQ. Mmatrix(2, p)) .OR. &
         .NOT.(Nmatrix(2, p) .EQ. Nmatrix(2, p))) THEN
        West = FF_current
        DO q = 1, l_nzeros
           Mmatrix(1, q) = Mmatrix(2, q)
           Nmatrix(1, q) = Nmatrix(2, q)
        END DO
        DEALLOCATE(xscaled)
        RETURN
     END IF

     IF ( (DABS(Mmatrix(2, p)) .GT. HUGE_LIMIT) .OR. &
          (DABS(Nmatrix(2, p)) .GT. HUGE_LIMIT) ) THEN
        West = FF_current
        DO q = 1, l_nzeros
           Mmatrix(1, q) = Mmatrix(2, q)
           Nmatrix(1, q) = Nmatrix(2, q)
        END DO
        DEALLOCATE(xscaled)
        RETURN
     END IF
  END DO

  IF (l_nzeros > 1) THEN
    IF (.NOT.(Nmatrix(2, l_nzeros) .EQ. Nmatrix(2, l_nzeros))) THEN
       West = Mmatrix(2, l_nzeros - 1) / Nmatrix(2, l_nzeros - 1)
    ELSE IF (DABS(Nmatrix(2, l_nzeros)) .LT. tinyDenom .OR. &
             DABS(Mmatrix(2, l_nzeros)) .GT. 1.0E300_C_DOUBLE) THEN
       West = Mmatrix(2, l_nzeros - 1) / Nmatrix(2, l_nzeros - 1)
    ELSE
       West = Mmatrix(2, l_nzeros) / Nmatrix(2, l_nzeros)
    END IF
  ELSE
    West = FF_current
  END IF

  DO q = 1, l_nzeros
     Mmatrix(1, q) = Mmatrix(2, q)
     Nmatrix(1, q) = Nmatrix(2, q)
  END DO
  DEALLOCATE(xscaled)

  RETURN
END SUBROUTINE accelerate
