program basic
   use discrete_weibull
   implicit none
   type(dweibull_fit_result) :: fit
   integer(i64) :: x(20)

   call seed_rng(12345)
   call rdweibull(x,0.6_dp,0.8_dp,.false.)
   fit = estdweibull(x,"ML")

   print '(a,2f12.6)', "type-I ML q,beta: ", fit%pars
   print '(a,f12.6)', "mean at fit:       ", Edweibull(fit%pars(1),fit%pars(2))

   call rdweibull3(x,0.3_dp,0.75_dp)
   fit = estdweibull3(x,"ML")
   print '(a,2f12.6)', "type-III ML c,beta:", fit%pars

contains
   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: n,i
      integer, allocatable :: s(:)
      call random_seed(size=n)
      allocate(s(n))
      do i=1,n
         s(i)=mod(seed+104729*i,2147483646)+1
      end do
      call random_seed(put=s)
   end subroutine seed_rng
end program basic
