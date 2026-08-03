! SPDX-License-Identifier: GPL-3.0-only
module indgen_moebius
  use indgen_kinds, only : dp
  use indgen_special, only : normal_quantile
  use indgen_types, only : indgen_success, indgen_invalid_argument, indgen_numerical_error
  implicit none
  private

  public :: moebius_full_stat

contains

  subroutine cdf_at_observations(x, fn_le, fn_lt, mass)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: fn_le(:), fn_lt(:), mass(:)
    integer :: n, i

    n = size(x)
    do i = 1, n
      fn_le(i) = real(count(x <= x(i)),dp)/real(n,dp)
      fn_lt(i) = real(count(x < x(i)),dp)/real(n,dp)
      mass(i) = fn_le(i)-fn_lt(i)
    end do
  end subroutine cdf_at_observations

  elemental function lg_integral(u) result(v)
    real(dp), intent(in) :: u
    real(dp) :: v, z

    if (u > 0.0_dp .and. u < 1.0_dp) then
      z = normal_quantile(u)
      v = exp(-0.5_dp*z*z)/sqrt(2.0_dp*acos(-1.0_dp))
    else
      v = 0.0_dp
    end if
  end function lg_integral

  elemental function le_integral(u) result(v)
    real(dp), intent(in) :: u
    real(dp) :: v

    if (u > 0.0_dp) then
      v = u-u*log(u)
    else
      v = 0.0_dp
    end if
  end function le_integral

  subroutine moebius_scores(x, score_s, score_g, score_e, sd_s, sd_g, sd_e)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: score_s(:), score_g(:), score_e(:)
    real(dp), intent(out) :: sd_s, sd_g, sd_e
    real(dp), allocatable :: f1(:), f0(:), mass(:)
    integer :: n

    n = size(x)
    allocate(f1(n),f0(n),mass(n))
    call cdf_at_observations(x,f1,f0,mass)
    score_s = 0.5_dp*(f1+f0-1.0_dp)
    score_g = (lg_integral(f1)-lg_integral(f0))/mass
    score_e = (le_integral(f1)-le_integral(f0))/mass-1.0_dp
    sd_s = sqrt(sum((score_s-sum(score_s)/real(n,dp))**2)/real(n,dp))
    sd_g = sqrt(sum((score_g-sum(score_g)/real(n,dp))**2)/real(n,dp))
    sd_e = sqrt(sum((score_e-sum(score_e)/real(n,dp))**2)/real(n,dp))
  end subroutine moebius_scores

  subroutine moebius_full_stat(x, spearman, vdw, savage, status)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: spearman, vdw, savage
    integer, intent(out) :: status
    real(dp), allocatable :: ss(:,:), sg(:,:), se(:,:)
    real(dp), allocatable :: sds(:), sdg(:), sde(:)
    real(dp) :: ps, pg, pe, ns, ng, ne
    integer :: n, d, i, j

    n = size(x,1)
    d = size(x,2)
    spearman = 0.0_dp
    vdw = 0.0_dp
    savage = 0.0_dp
    if (n < 2 .or. d < 2) then
      status = indgen_invalid_argument
      return
    end if

    allocate(ss(n,d),sg(n,d),se(n,d),sds(d),sdg(d),sde(d))
    do j = 1, d
      call moebius_scores(x(:,j),ss(:,j),sg(:,j),se(:,j),sds(j),sdg(j),sde(j))
    end do
    ns = product(sds)
    ng = product(sdg)
    ne = product(sde)
    if (min(ns,ng,ne) <= tiny(1.0_dp)) then
      status = indgen_numerical_error
      return
    end if

    do i = 1, n
      ps = 1.0_dp
      pg = 1.0_dp
      pe = 1.0_dp
      do j = 1, d
        ps = ps*ss(i,j)
        pg = pg*sg(i,j)
        pe = pe*se(i,j)
      end do
      spearman = spearman+ps
      vdw = vdw+pg
      savage = savage+pe
    end do
    spearman = spearman/real(n,dp)/ns
    vdw = vdw/real(n,dp)/ng
    savage = savage/real(n,dp)/ne
    status = indgen_success
  end subroutine moebius_full_stat

end module indgen_moebius
