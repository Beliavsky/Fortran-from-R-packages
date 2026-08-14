program test_solvers
   use isotone
   implicit none
   integer :: a(2,2), i
   real(dp) :: z(4), y(4), w(4), wm(4,4)
   type(isotone_solver_result) :: r
   a=reshape([1,3,2,4],[2,2])
   z=[2.0_dp,2.0_dp,3.0_dp,3.0_dp]
   y=[1.0_dp,3.0_dp,2.0_dp,5.0_dp]
   w=[1.0_dp,2.0_dp,1.5_dp,0.5_dp]
   call ls_solver(z,a,y,w,r); call good(r,'ls')
   call d_solver(z,a,y,w,r); call good(r,'l1')
   call p_solver(z,a,y,w,0.3_dp,0.7_dp,r); call good(r,'quantile')
   wm=0.0_dp; do i=1,4;wm(i,i)=w(i);end do
   call lf_solver(z,a,y,wm,r);call good(r,'gls')
   call s_solver([1.0_dp,1.0_dp,2.0_dp,2.0_dp],a,y,r);call good(r,'poisson')
   call o_solver(z,a,y,w,1.5_dp,r);call good(r,'lp')
   call a_solver(z,a,y,w,2.0_dp,1.0_dp,r);call good(r,'asyls')
   call e_solver(z,a,y,w,1.0e-4_dp,r);call good(r,'l1eps')
   call h_solver(z,a,y,w,0.5_dp,r);call good(r,'huber')
   call i_solver(z,a,y,w,0.5_dp,0.8_dp,r);call good(r,'silf')
   call m_solver(z,a,y,w,r);call good(r,'chebyshev')
   call f_solver(z,a,y,custom_loss,r);call good(r,'custom')
   print *, 'test_solvers: PASS'
contains
   subroutine good(r,name)
      type(isotone_solver_result),intent(in)::r
      character(*),intent(in)::name
      if(r%status/=ISO_SUCCESS) then
         print *,trim(name),' failed'
         error stop 1
      end if
      if(abs(r%x(1)-r%x(2))>1.0e-8_dp .or. abs(r%x(3)-r%x(4))>1.0e-8_dp) then
         print *,trim(name),' violated active equalities'
         error stop 1
      end if
   end subroutine good
   subroutine custom_loss(x,f,g)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::f,g(:)
      f=sum((x-y)**2);g=2.0_dp*(x-y)
   end subroutine custom_loss
end program
