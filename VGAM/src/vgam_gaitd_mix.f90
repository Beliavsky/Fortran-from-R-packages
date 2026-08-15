! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_gaitd_mix
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use vgam_kinds, only : dp
   use vgam_distributions, only : dpois_v, dnbinom_v
   use vgam_gaitd, only : gaitd_distribution_t
   implicit none
   private

   public :: gaitd_mix_poisson, gaitd_mix_negative_binomial

contains

   subroutine gaitd_mix_poisson(lambda_p, max_support, dist, a_mix, i_mix, d_mix, &
      pobs_mix, pstr_mix, pdip_mix, lambda_a, lambda_i, lambda_d, truncate, &
      a_mlm, pobs_mlm, i_mlm, pstr_mlm, d_mlm, pdip_mlm)
      real(dp), intent(in) :: lambda_p
      integer, intent(in) :: max_support
      type(gaitd_distribution_t), intent(out) :: dist
      integer, intent(in), optional :: a_mix(:), i_mix(:), d_mix(:), truncate(:)
      integer, intent(in), optional :: a_mlm(:), i_mlm(:), d_mlm(:)
      real(dp), intent(in), optional :: pobs_mix, pstr_mix, pdip_mix
      real(dp), intent(in), optional :: lambda_a, lambda_i, lambda_d
      real(dp), intent(in), optional :: pobs_mlm(:), pstr_mlm(:), pdip_mlm(:)
      real(dp) :: la, li, ld
      real(dp), allocatable :: parent(:), outer_a(:), outer_i(:), outer_d(:)
      integer :: k
      if (lambda_p < 0.0_dp .or. max_support < 0) then
         dist%status = 101
         return
      end if
      la = lambda_p; if (present(lambda_a)) la = lambda_a
      li = lambda_p; if (present(lambda_i)) li = lambda_i
      ld = lambda_p; if (present(lambda_d)) ld = lambda_d
      if (la < 0.0_dp .or. li < 0.0_dp .or. ld < 0.0_dp) then
         dist%status = 102
         return
      end if
      allocate(parent(max_support + 1), outer_a(max_support + 1), &
         outer_i(max_support + 1), outer_d(max_support + 1))
      do k = 0, max_support
         parent(k + 1) = dpois_v(k, lambda_p)
         outer_a(k + 1) = dpois_v(k, la)
         outer_i(k + 1) = dpois_v(k, li)
         outer_d(k + 1) = dpois_v(k, ld)
      end do
      call build_mix_distribution(parent, outer_a, outer_i, outer_d, 0, dist, &
         a_mix, i_mix, d_mix, pobs_mix, pstr_mix, pdip_mix, truncate, &
         a_mlm, pobs_mlm, i_mlm, pstr_mlm, d_mlm, pdip_mlm)
   end subroutine gaitd_mix_poisson

   subroutine gaitd_mix_negative_binomial(mu_p, size_p, max_support, dist, a_mix, i_mix, d_mix, &
      pobs_mix, pstr_mix, pdip_mix, mu_a, size_a, mu_i, size_i, mu_d, size_d, truncate, &
      a_mlm, pobs_mlm, i_mlm, pstr_mlm, d_mlm, pdip_mlm)
      real(dp), intent(in) :: mu_p, size_p
      integer, intent(in) :: max_support
      type(gaitd_distribution_t), intent(out) :: dist
      integer, intent(in), optional :: a_mix(:), i_mix(:), d_mix(:), truncate(:)
      integer, intent(in), optional :: a_mlm(:), i_mlm(:), d_mlm(:)
      real(dp), intent(in), optional :: pobs_mix, pstr_mix, pdip_mix
      real(dp), intent(in), optional :: mu_a, size_a, mu_i, size_i, mu_d, size_d
      real(dp), intent(in), optional :: pobs_mlm(:), pstr_mlm(:), pdip_mlm(:)
      real(dp) :: ma, mi, md, sa, si, sd, pp
      real(dp), allocatable :: parent(:), outer_a(:), outer_i(:), outer_d(:)
      integer :: k
      if (mu_p < 0.0_dp .or. size_p <= 0.0_dp .or. max_support < 0) then
         dist%status = 111
         return
      end if
      ma = mu_p; mi = mu_p; md = mu_p
      sa = size_p; si = size_p; sd = size_p
      if (present(mu_a)) ma = mu_a
      if (present(mu_i)) mi = mu_i
      if (present(mu_d)) md = mu_d
      if (present(size_a)) sa = size_a
      if (present(size_i)) si = size_i
      if (present(size_d)) sd = size_d
      if (min(ma, mi, md) < 0.0_dp .or. min(sa, si, sd) <= 0.0_dp) then
         dist%status = 112
         return
      end if
      allocate(parent(max_support + 1), outer_a(max_support + 1), &
         outer_i(max_support + 1), outer_d(max_support + 1))
      do k = 0, max_support
         pp = size_p/(size_p + mu_p)
         parent(k + 1) = dnbinom_v(k, size_p, pp)
         pp = sa/(sa + ma); outer_a(k + 1) = dnbinom_v(k, sa, pp)
         pp = si/(si + mi); outer_i(k + 1) = dnbinom_v(k, si, pp)
         pp = sd/(sd + md); outer_d(k + 1) = dnbinom_v(k, sd, pp)
      end do
      call build_mix_distribution(parent, outer_a, outer_i, outer_d, 0, dist, &
         a_mix, i_mix, d_mix, pobs_mix, pstr_mix, pdip_mix, truncate, &
         a_mlm, pobs_mlm, i_mlm, pstr_mlm, d_mlm, pdip_mlm)
   end subroutine gaitd_mix_negative_binomial

   subroutine build_mix_distribution(parent, outer_a, outer_i, outer_d, min_support, dist, &
      a_mix, i_mix, d_mix, pobs_mix, pstr_mix, pdip_mix, truncate, &
      a_mlm, pobs_mlm, i_mlm, pstr_mlm, d_mlm, pdip_mlm)
      real(dp), intent(in) :: parent(:), outer_a(:), outer_i(:), outer_d(:)
      integer, intent(in) :: min_support
      type(gaitd_distribution_t), intent(out) :: dist
      integer, intent(in), optional :: a_mix(:), i_mix(:), d_mix(:), truncate(:)
      integer, intent(in), optional :: a_mlm(:), i_mlm(:), d_mlm(:)
      real(dp), intent(in), optional :: pobs_mix, pstr_mix, pdip_mix
      real(dp), intent(in), optional :: pobs_mlm(:), pstr_mlm(:), pdip_mlm(:)
      integer, allocatable :: amix(:), imix(:), dmix(:), tr(:), amlm(:), imlm(:), dmlm(:)
      real(dp), allocatable :: pa(:), pi(:), pd(:), work(:), wa(:), wi(:), wd(:)
      real(dp) :: pobsm, pstrm, pdipm, sumt, suma, suma_prob, sumi_prob, sumd_prob
      real(dp) :: cdfmax, denom, tmp6, delta
      integer :: k, j, idx

      call copy_int(a_mix, amix); call copy_int(i_mix, imix); call copy_int(d_mix, dmix)
      call copy_int(truncate, tr); call copy_int(a_mlm, amlm); call copy_int(i_mlm, imlm)
      call copy_int(d_mlm, dmlm); call copy_real(pobs_mlm, pa); call copy_real(pstr_mlm, pi)
      call copy_real(pdip_mlm, pd)
      pobsm = 0.0_dp; pstrm = 0.0_dp; pdipm = 0.0_dp
      if (present(pobs_mix)) pobsm = pobs_mix
      if (present(pstr_mix)) pstrm = pstr_mix
      if (present(pdip_mix)) pdipm = pdip_mix
      if (size(parent) <= 0 .or. size(outer_a) /= size(parent) .or. size(outer_i) /= size(parent) .or. &
          size(outer_d) /= size(parent) .or. min(pobsm, pstrm, pdipm) < 0.0_dp .or. &
          max(pobsm, pstrm, pdipm) >= 1.0_dp) then
         dist%status = 120
         return
      end if
      if (size(amlm) /= size(pa) .or. size(imlm) /= size(pi) .or. size(dmlm) /= size(pd) .or. &
          any(pa < 0.0_dp) .or. any(pi < 0.0_dp) .or. any(pd < 0.0_dp) .or. &
          sum(pa) > 1.0_dp .or. sum(pi) > 1.0_dp .or. sum(pd) > 1.0_dp) then
         dist%status = 121
         return
      end if
      if (.not. valid_points([amix, imix, dmix, tr, amlm, imlm, dmlm], min_support, size(parent))) then
         dist%status = 122
         return
      end if
      if (has_duplicates([amix, imix, dmix, tr, amlm, imlm, dmlm])) then
         dist%status = 123
         return
      end if
      if (size(amix) == 0 .and. pobsm /= 0.0_dp .or. size(amix) > 0 .and. pobsm == 0.0_dp) then
         dist%status = 124
         return
      end if
      if (size(imix) == 0 .and. pstrm /= 0.0_dp .or. size(imix) > 0 .and. pstrm == 0.0_dp) then
         dist%status = 125
         return
      end if
      if (size(dmix) == 0 .and. pdipm /= 0.0_dp .or. size(dmix) > 0 .and. pdipm == 0.0_dp) then
         dist%status = 126
         return
      end if

      cdfmax = sum(parent)
      sumt = point_mass(parent, tr, min_support)
      suma = point_mass(parent, amlm, min_support) + point_mass(parent, amix, min_support)
      suma_prob = sum(pa); sumi_prob = sum(pi); sumd_prob = sum(pd)
      denom = cdfmax - sumt - suma
      tmp6 = 1.0_dp - suma_prob - sumi_prob - pobsm - pstrm + sumd_prob + pdipm
      if (denom <= tiny(1.0_dp) .or. tmp6 <= 0.0_dp) then
         dist%status = 127
         return
      end if
      delta = tmp6/denom
      allocate(work(size(parent))); work = delta*parent
      do j = 1, size(tr)
         work(point_index(tr(j), min_support)) = 0.0_dp
      end do
      do j = 1, size(amlm)
         work(point_index(amlm(j), min_support)) = pa(j)
      end do
      do j = 1, size(amix)
         work(point_index(amix(j), min_support)) = 0.0_dp
      end do
      do j = 1, size(imlm)
         idx = point_index(imlm(j), min_support); work(idx) = work(idx) + pi(j)
      end do
      do j = 1, size(dmlm)
         idx = point_index(dmlm(j), min_support); work(idx) = work(idx) - pd(j)
      end do

      call restricted_weights(outer_a, amix, min_support, wa)
      call restricted_weights(outer_i, imix, min_support, wi)
      call restricted_weights(outer_d, dmix, min_support, wd)
      if ((size(amix) > 0 .and. size(wa) == 0) .or. (size(imix) > 0 .and. size(wi) == 0) .or. &
          (size(dmix) > 0 .and. size(wd) == 0)) then
         dist%status = 128
         return
      end if
      do j = 1, size(amix)
         idx = point_index(amix(j), min_support); work(idx) = work(idx) + pobsm*wa(j)
      end do
      do j = 1, size(imix)
         idx = point_index(imix(j), min_support); work(idx) = work(idx) + pstrm*wi(j)
      end do
      do j = 1, size(dmix)
         idx = point_index(dmix(j), min_support); work(idx) = work(idx) - pdipm*wd(j)
      end do
      if (any(work < -1.0e-12_dp) .or. any(.not. ieee_is_finite(work))) then
         dist%status = 129
         return
      end if
      where (work < 0.0_dp) work = 0.0_dp
      if (abs(sum(work) - 1.0_dp) > 5.0e-10_dp) then
         dist%status = 130
         return
      end if
      call finalize_distribution(work, min_support, dist)
   end subroutine build_mix_distribution

   subroutine finalize_distribution(work, min_support, dist)
      real(dp), intent(in) :: work(:)
      integer, intent(in) :: min_support
      type(gaitd_distribution_t), intent(out) :: dist
      real(dp) :: x, ex2
      integer :: k
      dist%min_support = min_support
      dist%pmf = work
      allocate(dist%cdf(size(work)))
      dist%cdf(1) = work(1)
      do k = 2, size(work)
         dist%cdf(k) = min(1.0_dp, dist%cdf(k - 1) + work(k))
      end do
      dist%cdf(size(work)) = 1.0_dp
      dist%mean = 0.0_dp; ex2 = 0.0_dp
      do k = 1, size(work)
         x = real(min_support + k - 1, dp)
         dist%mean = dist%mean + x*work(k)
         ex2 = ex2 + x*x*work(k)
      end do
      dist%variance = max(0.0_dp, ex2 - dist%mean*dist%mean)
      dist%status = 0
   end subroutine finalize_distribution

   subroutine restricted_weights(pmf, points, min_support, w)
      real(dp), intent(in) :: pmf(:)
      integer, intent(in) :: points(:), min_support
      real(dp), allocatable, intent(out) :: w(:)
      real(dp) :: total
      integer :: j
      allocate(w(size(points)))
      if (size(points) == 0) return
      do j = 1, size(points)
         w(j) = pmf(point_index(points(j), min_support))
      end do
      total = sum(w)
      if (total <= tiny(1.0_dp)) then
         deallocate(w); allocate(w(0))
      else
         w = w/total
      end if
   end subroutine restricted_weights

   real(dp) function point_mass(pmf, points, min_support) result(s)
      real(dp), intent(in) :: pmf(:)
      integer, intent(in) :: points(:), min_support
      integer :: j
      s = 0.0_dp
      do j = 1, size(points)
         s = s + pmf(point_index(points(j), min_support))
      end do
   end function point_mass

   integer function point_index(point, min_support) result(idx)
      integer, intent(in) :: point, min_support
      idx = point - min_support + 1
   end function point_index

   logical function valid_points(points, min_support, n) result(ok)
      integer, intent(in) :: points(:), min_support, n
      ok = all(points >= min_support) .and. all(points <= min_support + n - 1)
   end function valid_points

   logical function has_duplicates(points) result(dup)
      integer, intent(in) :: points(:)
      integer :: i, j
      dup = .false.
      do i = 1, size(points) - 1
         do j = i + 1, size(points)
            if (points(i) == points(j)) then
               dup = .true.
               return
            end if
         end do
      end do
   end function has_duplicates

   subroutine copy_int(input, output)
      integer, intent(in), optional :: input(:)
      integer, allocatable, intent(out) :: output(:)
      if (present(input)) then
         output = input
      else
         allocate(output(0))
      end if
   end subroutine copy_int

   subroutine copy_real(input, output)
      real(dp), intent(in), optional :: input(:)
      real(dp), allocatable, intent(out) :: output(:)
      if (present(input)) then
         output = input
      else
         allocate(output(0))
      end if
   end subroutine copy_real

end module vgam_gaitd_mix
