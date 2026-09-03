program test_corstructures
   use geepack
   implicit none
   integer, parameter :: ng = 12
   integer, parameter :: s = 4
   integer, parameter :: n = ng * s
   real(dp) :: y(n)
   real(dp) :: x(n, 2)
   integer :: cs(ng)
   integer :: waves(n)
   real(dp) :: cor_param(n)
   real(dp) :: b0(2)
   real(dp), allocatable :: zuser(:, :)
   real(dp), allocatable :: zfixed(:, :)
   type(gee_spec) :: se
   type(gee_spec) :: su
   type(gee_spec) :: sa
   type(gee_spec) :: sf
   type(gee_spec) :: sun
   type(gee_result) :: fe
   type(gee_result) :: fu
   type(gee_result) :: fa
   type(gee_result) :: fa_irregular
   type(gee_result) :: ff
   type(gee_result) :: fun
   integer :: g
   integer :: j
   integer :: pos
   integer :: row
   integer :: k

   cs = s
   pos = 0
   do g = 1, ng
      do j = 1, s
         pos = pos + 1
         waves(pos) = j
         select case (j)
         case (1)
            cor_param(pos) = 0.0_dp
         case (2)
            cor_param(pos) = 1.0_dp
         case (3)
            cor_param(pos) = 3.0_dp
         case default
            cor_param(pos) = 6.0_dp
         end select
         x(pos, 1) = 1.0_dp
         x(pos, 2) = real(j - 1, dp)
         y(pos) = 0.8_dp + 0.35_dp * x(pos, 2) + 0.12_dp * real(g - 6, dp) + &
            0.03_dp * real(mod(2 * g + j, 5) - 2, dp)
      end do
   end do
   b0 = [0.8_dp, 0.3_dp]
   se%corstr = COR_EXCHANGEABLE
   se%scale_fixed = .true.
   se%tolerance = 1.0e-8_dp
   allocate(se%mean_links(1), se%variance_codes(1), se%scale_links(1))
   se%mean_links = LINK_IDENTITY
   se%variance_codes = VAR_GAUSSIAN
   se%scale_links = LINK_IDENTITY
   call fit_geese(y, x, cs, b0, se, fe, waves=waves)
   if (fe%error /= GEE_OK) error stop 1

   su = se
   su%corstr = COR_USERDEFINED
   allocate(zuser(ng * s * (s - 1) / 2, 1))
   zuser = 1.0_dp
   call fit_geese(y, x, cs, b0, su, fu, waves=waves, zcor=zuser)
   if (fu%error /= GEE_OK) error stop 2
   if (maxval(abs(fu%beta - fe%beta)) > 1.0e-8_dp) error stop 3
   if (abs(fu%alpha(1) - fe%alpha(1)) > 1.0e-8_dp) error stop 4

   sa = se
   sa%corstr = COR_AR1
   sa%max_iterations = 100
   call fit_geese(y, x, cs, b0, sa, fa, waves=waves)
   if (fa%error /= GEE_OK .or. abs(fa%alpha(1)) >= 1.0_dp) error stop 5
   call fit_geese(y, x, cs, b0, sa, fa_irregular, waves=waves, cor_param=cor_param)
   if (fa_irregular%error /= GEE_OK .or. abs(fa_irregular%alpha(1)) >= 1.0_dp) error stop 10
   if (abs(fa_irregular%alpha(1) - fa%alpha(1)) < 1.0e-6_dp) error stop 11

   sf = se
   sf%corstr = COR_FIXED
   allocate(zfixed(size(zuser, 1), 1))
   row = 0
   do g = 1, ng
      do j = 1, s - 1
         do k = j + 1, s
            row = row + 1
            zfixed(row, 1) = 0.3_dp ** (k - j)
         end do
      end do
   end do
   call fit_geese(y, x, cs, b0, sf, ff, waves=waves, zcor=zfixed)
   if (ff%error /= GEE_OK) error stop 6
   if (size(ff%alpha) /= 1 .or. abs(ff%alpha(1) - 1.0_dp) > 1.0e-14_dp) error stop 7
   sun = se
   sun%corstr = COR_UNSTRUCTURED
   call fit_geese(y, x, cs, b0, sun, fun, waves=waves)
   if (fun%error /= GEE_OK .or. size(fun%alpha) /= 6) error stop 8
   if (any(abs(fun%alpha) > 1.0_dp)) error stop 9
   print *, 'test_corstructures: PASS'
end program test_corstructures
