program branin_mbo
  use mlrmbo
  implicit none
  type(mbo_space) :: space
  type(mbo_control) :: control
  type(mbo_result) :: result

  call init_space(space,[mbo_real,mbo_real],[-5.0_dp,0.0_dp],[10.0_dp,15.0_dp])
  call init_control(control,1,[.true.])
  control%n_init=10
  control%max_iter=6
  control%infill_criterion=crit_ei
  control%focus_restarts=2
  control%focus_iterations=3
  control%focus_points=120
  control%km_multistart=1
  control%km_pop_size=10
  control%km_max_iter=100
  control%seed=20260813_i8
  call seed_intrinsic(2026)
  call mbo(space,branin,control,result)
  write(*,'(a,2f12.6)') 'best x: ',result%best_x
  write(*,'(a,f14.6)') 'best y: ',result%best_y(1)
contains
  subroutine branin(x,y)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::y(:)
    real(dp)::a,b,c,r,s,t
    a=1.0_dp;b=5.1_dp/(4.0_dp*pi_dp*pi_dp);c=5.0_dp/pi_dp
    r=6.0_dp;s=10.0_dp;t=1.0_dp/(8.0_dp*pi_dp)
    y(1)=a*(x(2)-b*x(1)*x(1)+c*x(1)-r)**2+s*(1.0_dp-t)*cos(x(1))+s
  end subroutine branin
  subroutine seed_intrinsic(a)
    integer,intent(in)::a
    integer::n,i
    integer,allocatable::z(:)
    call random_seed(size=n);allocate(z(n));do i=1,n;z(i)=a+29*i;end do;call random_seed(put=z)
  end subroutine seed_intrinsic
end program branin_mbo
