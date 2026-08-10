! SPDX-License-Identifier: GPL-3.0-only
program test_scs
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use scs_kinds, only : dp
   use scs_types
   use scs_sparse, only : dense_to_csc, dense_upper_to_csc
   use scs_solver, only : scs
   use scs_cones, only : project_primal_cone
   implicit none
   call test_example
   call test_unbounded
   call test_socp
   call test_qp
   call test_psd
   call test_cone_projections
   call test_acceleration
   print '(a)', 'All SCS Fortran tests passed.'
contains
   subroutine assert_close(x,y,tol,msg)
      real(dp),intent(in)::x(:),y(:),tol
      character(*),intent(in)::msg
      if(size(x)/=size(y) .or. any(.not.ieee_is_finite(x)) .or. any(.not.ieee_is_finite(y)) .or. &
         any(abs(x-y)>tol))then
         print *,trim(msg),x,y
         error stop 1
      end if
   end subroutine

   subroutine test_example
      type(scs_data)::d
      type(scs_cone)::k
      type(scs_settings)::st
      type(scs_solution)::sol
      type(scs_info)::info
      real(dp)::a(2,1)
      a(:,1)=[1.0_dp,1.0_dp]
      d%m=2;d%n=1;call dense_to_csc(a,d%A);d%b=[1.0_dp,1.0_dp];d%c=[1.0_dp]
      k%z=2;st%eps_rel=1e-3_dp;st%eps_abs=1e-3_dp;st%max_iters=50
      call scs(d,k,st,sol,info)
      if(info%status_val/=scs_solved)then;print *,trim(info%status),info%status_val;error stop 'example status';end if
      call assert_close(sol%x,[1.0_dp],1e-6_dp,'example x')
      call assert_close(sol%y,[-0.5_dp,-0.5_dp],1e-5_dp,'example y')
      if(trim(info%lin_sys_solver)/='native-sparse-qdldl-natural') error stop 'linear solver name'
      if(info%kkt_nnz<=0 .or. info%symbolic_analyses/=1 .or. info%factorizations<1) error stop 'linear solver stats'
   end subroutine

   subroutine test_unbounded
      type(scs_data)::d;type(scs_cone)::k;type(scs_settings)::st;type(scs_solution)::sol;type(scs_info)::info
      real(dp)::a(2,1)
      a(:,1)=[0.5_dp,2.0_dp];d%m=2;d%n=1;call dense_to_csc(a,d%A);d%b=[3.0_dp,1.0_dp];d%c=[1.0_dp];k%l=2
      st%max_iters=5000;call scs(d,k,st,sol,info)
      if (info%status_val /= scs_unbounded) then
         print *,trim(info%status),info%status_val,info%res_unbdd_a,info%res_unbdd_p
         error stop 'unbounded status'
      end if
   end subroutine

   subroutine test_socp
      type(scs_data)::d;type(scs_cone)::k;type(scs_settings)::st;type(scs_solution)::sol;type(scs_info)::info
      real(dp)::a(4,3)
      a=0.0_dp;a(1,1)=1.0_dp;a(2,1)=-1.0_dp;a(3,2)=-1.0_dp;a(4,3)=-1.0_dp
      d%m=4;d%n=3;call dense_to_csc(a,d%A);d%b=[sqrt(2.0_dp),0.0_dp,0.0_dp,0.0_dp];d%c=[1.0_dp,1.0_dp,1.0_dp]
      k%z=1;allocate(k%q(1),source=3);st%max_iters=100000;call scs(d,k,st,sol,info)
      if(info%status_val/=scs_solved)then;print *,trim(info%status),info%status_val;error stop 'socp status';end if
      call assert_close(sol%x,[sqrt(2.0_dp),-1.0_dp,-1.0_dp],1e-4_dp,'socp x')
   end subroutine

   subroutine test_qp
      type(scs_data)::d;type(scs_cone)::k;type(scs_settings)::st;type(scs_solution)::sol;type(scs_info)::info
      real(dp)::a(2,1),p(1,1)
      a(:,1)=[-1.0_dp,1.0_dp];p(1,1)=2.0_dp
      d%m=2;d%n=1;call dense_to_csc(a,d%A);call dense_upper_to_csc(p,d%P);d%has_p=.true.
      d%b=[0.0_dp,1.0_dp];d%c=[-1.0_dp];k%l=2
      st%max_iters=10000;call scs(d,k,st,sol,info)
      if(info%status_val/=scs_solved)then;print *,trim(info%status),info%status_val;error stop 'qp status';end if
      call assert_close(sol%x,[0.5_dp],2e-4_dp,'qp x')
   end subroutine

   subroutine test_psd
      type(scs_data) :: d
      type(scs_cone) :: k
      type(scs_settings) :: st
      type(scs_solution) :: sol
      type(scs_info) :: info
      real(dp) :: a(9,3)
      integer :: ii(25), jj(25), t
      real(dp) :: vv(25)
      a=0.0_dp
      ii=[1,2,3,4,5,7,8,9,1,2,3,5,6,7,8,9,1,2,3,4,5,6,7,8,9]
      jj=[1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3]
      vv=[-7.0_dp,-15.556349186104_dp,3.0_dp,-21.0_dp,-15.556349186104_dp,10.0_dp, &
          11.3137084989848_dp,5.0_dp,7.0_dp,-25.4558441227157_dp,8.0_dp, &
          14.142135623731_dp,22.6274169979695_dp,-10.0_dp,-14.142135623731_dp,3.0_dp, &
          -2.0_dp,-11.3137084989848_dp,1.0_dp,-5.0_dp,2.82842712474619_dp, &
          -24.0416305603426_dp,-6.0_dp,11.3137084989848_dp,6.0_dp]
      do t=1,size(vv)
         a(ii(t),jj(t))=vv(t)
      end do
      d%m=9; d%n=3; call dense_to_csc(a,d%A)
      d%b=[33.0_dp,-12.7279220613579_dp,26.0_dp,14.0_dp,12.7279220613579_dp, &
           56.5685424949238_dp,91.0_dp,14.142135623731_dp,15.0_dp]
      d%c=[1.0_dp,-1.0_dp,1.0_dp]
      allocate(k%s(2)); k%s=[2,3]
      st%max_iters=100000
      call scs(d,k,st,sol,info)
      if(info%status_val/=scs_solved)then
         print *,trim(info%status),info%status_val
         error stop 'psd status'
      end if
      call assert_close(sol%x,[-0.36778629730617_dp,1.89834935806707_dp,-0.88746074956525_dp], &
                        1.0e-4_dp,'psd x')
      call assert_close(sol%y,[0.00396134156379081_dp,-0.00613634113889057_dp, &
                        0.00475276897567086_dp,0.0558030988006113_dp,-0.00340904764830297_dp, &
                        0.0342440719001259_dp,0.000104130479854425_dp,-0.00147926212203912_dp, &
                        0.0105070907711366_dp],1.0e-4_dp,'psd y')
      if(abs(info%pobj-(-3.1535964049385_dp))>1.0e-4_dp) error stop 'psd pobj'
      if(abs(info%dobj-(-3.1535519733301_dp))>1.0e-4_dp) error stop 'psd dobj'
   end subroutine test_psd

   subroutine test_cone_projections
      type(scs_cone) :: k
      real(dp) :: x(3), y(3)

      k%ep=1
      x=[1.0_dp,-0.5_dp,0.1_dp]
      call project_primal_cone(x,k)
      y=x; call project_primal_cone(y,k)
      call assert_close(y,x,5.0e-10_dp,'exp projection idempotence')

      k=scs_cone(); k%ed=1
      x=[-1.0_dp,0.25_dp,-0.5_dp]
      call project_primal_cone(x,k)
      y=x; call project_primal_cone(y,k)
      call assert_close(y,x,5.0e-10_dp,'dual exp projection idempotence')

      k=scs_cone(); allocate(k%p(1)); k%p=[0.5_dp]
      x=[-1.0_dp,2.0_dp,3.0_dp]
      call project_primal_cone(x,k)
      y=x; call project_primal_cone(y,k)
      call assert_close(y,x,5.0e-9_dp,'power projection idempotence')

      k=scs_cone(); allocate(k%p(1)); k%p=[-0.5_dp]
      x=[1.0_dp,-2.0_dp,2.5_dp]
      call project_primal_cone(x,k)
      y=x; call project_primal_cone(y,k)
      call assert_close(y,x,5.0e-9_dp,'dual power projection idempotence')

      k=scs_cone(); k%bsize=3; allocate(k%bl(2),k%bu(2))
      k%bl=[-1.0_dp,0.0_dp]; k%bu=[1.0_dp,2.0_dp]
      x=[0.5_dp,5.0_dp,-3.0_dp]
      call project_primal_cone(x,k)
      if(x(1)<-1.0e-12_dp .or. x(2)<-x(1)-1.0e-10_dp .or. x(2)>x(1)+1.0e-10_dp .or. &
         x(3)<-1.0e-10_dp .or. x(3)>2.0_dp*x(1)+1.0e-10_dp) error stop 'box projection feasibility'
      y=x; call project_primal_cone(y,k)
      call assert_close(y,x,5.0e-10_dp,'box projection idempotence')
   end subroutine test_cone_projections

   subroutine test_acceleration
      type(scs_data)::d
      type(scs_cone)::k
      type(scs_settings)::st
      type(scs_solution)::sol
      type(scs_info)::info
      real(dp)::a(4,3)
      a=0.0_dp; a(1,1)=1.0_dp; a(2,1)=-1.0_dp; a(3,2)=-1.0_dp; a(4,3)=-1.0_dp
      d%m=4; d%n=3; call dense_to_csc(a,d%A)
      d%b=[sqrt(2.0_dp),0.0_dp,0.0_dp,0.0_dp]; d%c=[1.0_dp,1.0_dp,1.0_dp]
      k%z=1; allocate(k%q(1),source=3)
      st%max_iters=100000; st%acceleration_lookback=5; st%acceleration_interval=5
      call scs(d,k,st,sol,info)
      if(info%status_val/=scs_solved)then
         print *,trim(info%status),info%status_val
         error stop 'acceleration status'
      end if
      call assert_close(sol%x,[sqrt(2.0_dp),-1.0_dp,-1.0_dp],1.0e-4_dp,'accelerated socp x')
   end subroutine test_acceleration
end program test_scs
