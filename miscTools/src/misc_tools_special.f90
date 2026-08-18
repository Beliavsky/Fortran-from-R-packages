module misc_tools_special
   use misc_tools_kinds, only : dp
   implicit none
   private

   public :: regularized_beta, student_t_two_sided_p

contains

   real(dp) function beta_cf(a,b,x,status) result(h)
      real(dp), intent(in) :: a,b,x
      integer, intent(out), optional :: status
      integer, parameter :: maxit = 400
      real(dp), parameter :: eps = 8.0_dp*epsilon(1.0_dp)
      real(dp), parameter :: fpmin = tiny(1.0_dp)/epsilon(1.0_dp)
      real(dp) :: qab,qap,qam,c,d,aa,del
      integer :: m,m2

      if (present(status)) status = 0
      qab = a+b
      qap = a+1.0_dp
      qam = a-1.0_dp
      c = 1.0_dp
      d = 1.0_dp-qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d

      do m = 1, maxit
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x / &
              ((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c

         aa = -(a+real(m,dp))*(qab+real(m,dp))*x / &
              ((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= eps) return
      end do
      if (present(status)) status = 1
   end function beta_cf

   real(dp) function regularized_beta(x,a,b,status) result(p)
      real(dp), intent(in) :: x,a,b
      integer, intent(out), optional :: status
      real(dp) :: logbt,bt,cf
      integer :: istat

      if (present(status)) status = 0
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         p = 0.0_dp
         if (present(status)) status = 2
         return
      end if
      if (x <= 0.0_dp) then
         p = 0.0_dp
         return
      end if
      if (x >= 1.0_dp) then
         p = 1.0_dp
         return
      end if

      logbt = log_gamma(a+b)-log_gamma(a)-log_gamma(b) + &
              a*log(x)+b*log(1.0_dp-x)
      if (logbt < log(tiny(1.0_dp))) then
         bt = 0.0_dp
      else
         bt = exp(logbt)
      end if

      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
         cf = beta_cf(a,b,x,istat)
         p = bt*cf/a
      else
         cf = beta_cf(b,a,1.0_dp-x,istat)
         p = 1.0_dp-bt*cf/b
      end if
      p = min(1.0_dp,max(0.0_dp,p))
      if (present(status)) status = istat
   end function regularized_beta

   real(dp) function student_t_two_sided_p(t,df) result(p)
      real(dp), intent(in) :: t,df
      real(dp) :: x
      if (df <= 0.0_dp) then
         p = -1.0_dp
      else
         x = df/(df+t*t)
         p = regularized_beta(x,0.5_dp*df,0.5_dp)
      end if
   end function student_t_two_sided_p

end module misc_tools_special
