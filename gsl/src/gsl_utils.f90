module gsl_utils
  use iso_c_binding, only: c_double, c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use gsl_special, only: elljac
  implicit none
  private

  public :: strictify, gsl_poly
  public :: gsl_sn, gsl_cn, gsl_dn, gsl_ns, gsl_nc, gsl_nd
  public :: gsl_sc, gsl_sd, gsl_cs, gsl_cd, gsl_ds, gsl_dc

contains

  subroutine strictify(x, status)
    real(c_double), intent(inout) :: x(:)
    integer(c_int), intent(in) :: status(:)
    integer :: i
    if (size(x) /= size(status)) error stop 'strictify: size mismatch'
    do i = 1, size(x)
      if (status(i) /= 0_c_int) x(i) = ieee_value(x(i), ieee_quiet_nan)
    end do
  end subroutine strictify

  pure function gsl_poly(c, x) result(y)
    real(c_double), intent(in) :: c(:)
    real(c_double), intent(in) :: x(:)
    real(c_double) :: y(size(x))
    integer :: i, j
    if (size(c) == 0) then
      y = 0.0_c_double
      return
    end if
    do i = 1, size(x)
      y(i) = c(size(c))
      do j = size(c) - 1, 1, -1
        y(i) = y(i) * x(i) + c(j)
      end do
    end do
  end function gsl_poly

  function gsl_sn(z, m) result(v)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double) :: v
    complex(c_double) :: sn, cn, dn
    call jacobi_complex(z, m, sn, cn, dn)
    v = sn
  end function gsl_sn

  function gsl_cn(z, m) result(v)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double) :: v
    complex(c_double) :: sn, cn, dn
    call jacobi_complex(z, m, sn, cn, dn)
    v = cn
  end function gsl_cn

  function gsl_dn(z, m) result(v)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double) :: v
    complex(c_double) :: sn, cn, dn
    call jacobi_complex(z, m, sn, cn, dn)
    v = dn
  end function gsl_dn

  function gsl_ns(z, m) result(v)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double) :: v
    v = 1.0_c_double / gsl_sn(z, m)
  end function gsl_ns

  function gsl_nc(z, m) result(v)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double) :: v
    v = 1.0_c_double / gsl_cn(z, m)
  end function gsl_nc

  function gsl_nd(z, m) result(v)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double) :: v
    v = 1.0_c_double / gsl_dn(z, m)
  end function gsl_nd

  function gsl_sc(z, m) result(v)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double) :: v
    v = gsl_sn(z, m) / gsl_cn(z, m)
  end function gsl_sc

  function gsl_sd(z, m) result(v)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double) :: v
    v = gsl_sn(z, m) / gsl_dn(z, m)
  end function gsl_sd

  function gsl_cs(z, m) result(v)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double) :: v
    v = gsl_cn(z, m) / gsl_sn(z, m)
  end function gsl_cs

  function gsl_cd(z, m) result(v)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double) :: v
    v = gsl_cn(z, m) / gsl_dn(z, m)
  end function gsl_cd

  function gsl_ds(z, m) result(v)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double) :: v
    v = gsl_dn(z, m) / gsl_sn(z, m)
  end function gsl_ds

  function gsl_dc(z, m) result(v)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double) :: v
    v = gsl_dn(z, m) / gsl_cn(z, m)
  end function gsl_dc

  subroutine jacobi_complex(z, m, sn, cn, dn)
    complex(c_double), intent(in) :: z
    real(c_double), intent(in) :: m
    complex(c_double), intent(out) :: sn, cn, dn
    real(c_double), target :: ur(1), ui(1), mr(1), mi(1)
    real(c_double), target :: sr(1), cr(1), dr(1), si(1), ci(1), di(1)
    integer(c_int), target :: status_r(1), status_i(1)
    real(c_double) :: den

    ur(1) = real(z, c_double)
    ui(1) = aimag(z)
    mr(1) = m
    mi(1) = 1.0_c_double - m
    call elljac(ur, mr, sr, cr, dr, status_r)
    if (abs(ui(1)) <= tiny(1.0_c_double)) then
      sn = cmplx(sr(1), 0.0_c_double, c_double)
      cn = cmplx(cr(1), 0.0_c_double, c_double)
      dn = cmplx(dr(1), 0.0_c_double, c_double)
      return
    end if
    call elljac(ui, mi, si, ci, di, status_i)
    den = ci(1)**2 + m * sr(1)**2 * si(1)**2
    sn = cmplx(sr(1) * di(1), cr(1) * dr(1) * si(1) * ci(1), c_double) / den
    cn = cmplx(cr(1) * ci(1), -sr(1) * dr(1) * si(1) * di(1), c_double) / den
    dn = cmplx(dr(1) * ci(1) * di(1), -m * sr(1) * cr(1) * si(1), c_double) / den
  end subroutine jacobi_complex

end module gsl_utils
