module derivative_context
   use gkwdist_kinds, only : dp
   use gkwdist_core, only : family_nll
   implicit none
   integer :: current_family=1
   real(dp),allocatable :: current_data(:)
contains
   function objective(x) result(v)
      real(dp),intent(in)::x(:)
      real(dp)::v
      v=family_nll(current_family,x,current_data)
   end function objective
end module derivative_context

program test_derivatives
   use gkwdist
   use gkwdist_core, only : family_derivatives
   use derivative_context
   use numderiv, only : grad,hessian,deriv_options,hessian_options
   implicit none
   integer :: fails,fam,n,status
   real(dp),allocatable :: par(:),ga(:),ha(:,:),gn(:),hn(:,:)
   real(dp) :: nll,gtol,htol
   type(deriv_options) :: hopts
   fails=0
   current_data=[0.12_dp,0.19_dp,0.27_dp,0.35_dp,0.48_dp,0.61_dp,0.73_dp,0.86_dp]
   hopts=hessian_options()
   hopts%r=5
   do fam=1,7
      select case(fam)
      case(1); par=[1.7_dp,2.2_dp,1.3_dp,0.7_dp,1.25_dp]
      case(2); par=[1.7_dp,2.2_dp,1.3_dp,0.7_dp]
      case(3); par=[1.7_dp,2.2_dp,0.7_dp,1.25_dp]
      case(4); par=[1.7_dp,2.2_dp,1.25_dp]
      case(5); par=[1.3_dp,0.7_dp,1.25_dp]
      case(6); par=[1.7_dp,2.2_dp]
      case(7); par=[1.3_dp,0.7_dp]
      end select
      n=size(par); allocate(ga(n),ha(n,n),gn(n))
      current_family=fam
      call family_derivatives(fam,par,current_data,nll,ga,ha)
      call grad(objective,par,gn,status=status)
      call hessian(objective,par,hn,options=hopts,status=status)
      gtol=3.0e-6_dp; htol=3.0e-4_dp
      if(maxval(abs(ga-gn))>gtol) then
         print '(a,i0,es12.4)','gradient mismatch family ',fam,maxval(abs(ga-gn)); fails=fails+1
      end if
      if(maxval(abs(ha-hn))>htol) then
         print '(a,i0,es12.4)','hessian mismatch family ',fam,maxval(abs(ha-hn)); fails=fails+1
      end if
      deallocate(par,ga,ha,gn,hn)
   end do
   if(fails==0) then
      print '(a)','test_derivatives: PASS'
   else
      print '(a,i0)','test_derivatives: FAIL ',fails
      error stop 1
   end if
end program test_derivatives
