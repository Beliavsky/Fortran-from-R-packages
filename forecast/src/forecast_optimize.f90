module forecast_optimize
   use forecast_kinds, only : dp
   implicit none
   private
   public :: pattern_search, golden_minimize
   abstract interface
      function objective_vec(x) result(f)
         import dp
         real(dp),intent(in) :: x(:)
         real(dp) :: f
      end function
      function objective_scalar(x) result(f)
         import dp
         real(dp),intent(in) :: x
         real(dp) :: f
      end function
   end interface
contains
   subroutine pattern_search(fun,x,lower,upper,maxit,tol,fval)
      procedure(objective_vec) :: fun
      real(dp),intent(inout) :: x(:)
      real(dp),intent(in) :: lower(:),upper(:)
      integer,intent(in),optional :: maxit
      real(dp),intent(in),optional :: tol
      real(dp),intent(out),optional :: fval
      real(dp),allocatable :: step(:),trial(:)
      real(dp) :: f,ft,eps
      integer :: it,mi,j,sgn
      mi=600
      if(present(maxit)) mi=maxit
      eps=1.0e-6_dp
      if(present(tol)) eps=tol
      allocate(step(size(x)),trial(size(x)))
      x=max(lower,min(upper,x))
      step=max(0.05_dp*(upper-lower),1.0e-4_dp)
      f=fun(x)
      do it=1,mi
         do j=1,size(x)
            do sgn=-1,1,2
               trial=x
               trial(j)=max(lower(j),min(upper(j),x(j)+real(sgn,dp)*step(j)))
               ft=fun(trial)
               if(ft<f) then
               x=trial
               f=ft
               end if
            end do
         end do
         step=step*0.72_dp
         if(maxval(step)<=eps) exit
      end do
      if(present(fval)) fval=f
   end subroutine

   function golden_minimize(fun,a,b,tol) result(xmin)
      procedure(objective_scalar) :: fun
      real(dp),intent(in) :: a,b
      real(dp),intent(in),optional :: tol
      real(dp) :: xmin,l,r,x1,x2,f1,f2,eps,gr
      eps=1.0e-7_dp
      if(present(tol)) eps=tol
      gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
      l=a
      r=b
      x1=r-gr*(r-l)
      x2=l+gr*(r-l)
      f1=fun(x1)
      f2=fun(x2)
      do while(r-l>eps*(1.0_dp+abs(l)+abs(r)))
         if(f1<f2) then
            r=x2
            x2=x1
            f2=f1
            x1=r-gr*(r-l)
            f1=fun(x1)
         else
            l=x1
            x1=x2
            f1=f2
            x2=l+gr*(r-l)
            f2=fun(x2)
         end if
      end do
      xmin=0.5_dp*(l+r)
   end function
end module forecast_optimize
