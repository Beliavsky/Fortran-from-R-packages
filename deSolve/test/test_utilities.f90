program test_utilities
  use desolve_kinds, only : dp
  use desolve_types, only : ode_result
  use desolve_utilities, only : forcing_table,event_table,history_buffer,sparsity_pattern, &
      brent_root,iterate_map,sparsity_2d,sparsity_3d,hermite_value,hermite_deriv,clean_event_times
  implicit none
  type(forcing_table)::ft
  type(event_table)::et
  type(history_buffer)::hb
  type(sparsity_pattern)::sp
  type(ode_result)::it
  real(dp),allocatable::v(:),ct(:)
  real(dp)::y(2),root
  integer::status
  ft%time=[0.0_dp,1.0_dp,2.0_dp];allocate(ft%value(2,3));ft%value(1,:)=[0.0_dp,2.0_dp,4.0_dp]
  ft%value(2,:)=[1.0_dp,3.0_dp,5.0_dp];v=ft%eval(0.25_dp)
  if(maxval(abs(v-[0.5_dp,1.5_dp]))>1e-14_dp)error stop 'forcing'
  allocate(et%time(3),et%state(3),et%value(3),et%method(3));et%time=1.0_dp;et%state=[1,2,1]
  et%value=[3.0_dp,2.0_dp,2.0_dp];et%method=[1,2,3];y=[1.0_dp,4.0_dp];call et%apply(1.0_dp,y)
  if(maxval(abs(y-[6.0_dp,6.0_dp]))>1e-14_dp)error stop 'events'
  call hb%init(1,10);call hb%append(0.0_dp,[0.0_dp],[1.0_dp]);call hb%append(1.0_dp,[1.0_dp],[1.0_dp]);v=hb%lag_value(0.3_dp)
  if(abs(v(1)-0.3_dp)>1e-14_dp)error stop 'history'
  if(abs(hermite_value(0.0_dp,1.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp,0.4_dp)-0.4_dp)>1e-14_dp)error stop 'hermite'
  if(abs(hermite_deriv(0.0_dp,1.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp,0.4_dp)-1.0_dp)>1e-14_dp)error stop 'dhermite'
  root=brent_root(f,0.0_dp,2.0_dp,tol=1e-12_dp,status=status)
  if(status/=0.or.abs(root-sqrt(2.0_dp))>1e-10_dp)error stop 'brent'
  it=iterate_map(map,[1.0_dp],[0.0_dp,1.0_dp,2.0_dp],nsteps=2)
  if(abs(it%y(1,3)-16.0_dp)>1e-14_dp)error stop 'iteration'
  sp=sparsity_2d(1,2,2,.false.,.false.);if(size(sp%row_ptr)/=5.or.size(sp%col_ind)/=12)error stop 'sparsity2d'
  sp=sparsity_3d(1,2,2,2,.false.,.false.,.false.);if(size(sp%row_ptr)/=9)error stop 'sparsity3d'
  ct=clean_event_times([0.0_dp,1.0_dp,2.0_dp],[0.5_dp,1.0_dp]);if(size(ct)/=4)error stop 'clean events'
  print *,'test_utilities: PASS'
contains
  function f(x)result(vv)
    real(dp),intent(in)::x;real(dp)::vv;vv=x*x-2.0_dp
  end function f
  subroutine map(t,x,xn)
    real(dp),intent(in)::t,x(:);real(dp),intent(out)::xn(:);xn=2.0_dp*x;if(t< -huge(1.0_dp))stop
  end subroutine map
end program test_utilities
