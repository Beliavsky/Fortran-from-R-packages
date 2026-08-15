program test_core
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
   use rmutil
   implicit none
   real(dp) :: ans, eps, infv
   real(dp), allocatable :: nodes(:), weights(:), y(:,:), rk(:), aligned(:), c(:,:)
   integer :: m, info

   call check(abs(normal_quantile(normal_cdf(0.7_dp))-0.7_dp) < 2.0e-12_dp, "normal round trip")

   ans = integrate_romberg(fsin,0.0_dp,pi,1.0e-10_dp)
   call check(abs(ans-2.0_dp) < 1.0e-9_dp, "Romberg sin")

   ans = integrate_2d(fxy,0.0_dp,1.0_dp,0.0_dp,1.0_dp)
   call check(abs(ans-0.25_dp) < 1.0e-12_dp, "2D integration")

   m = 101
   info = 4
   eps = 1.0e-10_dp
   call toms614_integrate(fsquare,0.0_dp,1.0_dp,1.0_dp,m,0.0_dp,eps,info,ans)
   call check(abs(ans-1.0_dp/3.0_dp) < 1.0e-7_dp, "TOMS614 finite integral")

   call gauss_hermite(12,nodes,weights,status=info)
   call check(abs(sum(weights)-1.0_dp) < 2.0e-15_dp, "Gauss-Hermite weights")
   call check(abs(sum(weights*nodes*nodes)-1.0_dp) < 2.0e-12_dp, "Gauss-Hermite variance")

   rk = runge_kutta(fode,1.0_dp,[0.0_dp,0.1_dp,0.2_dp,0.3_dp,0.4_dp,0.5_dp])
   call check(abs(rk(6)-exp(0.5_dp)) < 1.0e-5_dp, "RK4")

   y = lin_diff_eqn(reshape([1.0_dp,0.0_dp,0.0_dp,-2.0_dp],[2,2]), &
      [1.0_dp,2.0_dp],[0.0_dp,0.5_dp])
   call check(abs(y(2,1)-exp(0.5_dp)) < 2.0e-12_dp, "matrix exponential 1")
   call check(abs(y(2,2)-2.0_dp*exp(-1.0_dp)) < 2.0e-12_dp, "matrix exponential 2")

   aligned = gettvc([1.0_dp,2.0_dp,4.0_dp,1.5_dp,3.0_dp], &
      [0.5_dp,2.0_dp,3.0_dp,1.0_dp,2.0_dp], [10.0_dp,20.0_dp,30.0_dp,4.0_dp,5.0_dp], &
      [3,2],[3,2],.true.)
   call check(maxval(abs(aligned-[10.0_dp,20.0_dp,30.0_dp,4.0_dp,5.0_dp])) < 1.0e-14_dp, "gettvc ties")

   c = contrast_mean(3)
   call check(maxval(abs(c-reshape([1.0_dp,0.0_dp,-1.0_dp,0.0_dp,1.0_dp,-1.0_dp],[3,2]))) < &
      1.0e-14_dp, "contr.mean")

   infv = ieee_value(0.0_dp,ieee_positive_inf)
   ans = integrate_romberg(fnormal,-infv,infv,1.0e-8_dp)
   call check(abs(ans-1.0_dp) < 2.0e-6_dp, "Romberg infinite normal")

   print *, "test_core: PASS"
contains
   function fsin(x) result(v)
      real(dp),intent(in)::x; real(dp)::v
      v=sin(x)
   end function fsin
   function fxy(x,y) result(v)
      real(dp),intent(in)::x,y; real(dp)::v
      v=x*y
   end function fxy
   function fsquare(x) result(v)
      real(dp),intent(in)::x; real(dp)::v
      v=x*x
   end function fsquare
   function fode(y,x) result(v)
      real(dp),intent(in)::y,x; real(dp)::v
      v=y+0.0_dp*x
   end function fode
   function fnormal(x) result(v)
      real(dp),intent(in)::x; real(dp)::v
      v=normal_pdf(x)
   end function fnormal
   subroutine check(ok,msg)
      logical,intent(in)::ok; character(*),intent(in)::msg
      if(.not.ok)then
         print *, "FAIL: ", trim(msg)
         error stop 1
      end if
   end subroutine check
end program test_core
