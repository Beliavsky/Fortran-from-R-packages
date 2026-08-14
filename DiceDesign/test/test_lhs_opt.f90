program test_lhs_opt
  use iso_fortran_env, only : int64
  use dicedesign, only : dp, lhs_design, discrepancy_value, phi_p, lhs_optimization_result, &
    discrep_ese_lhs, maximin_ese_lhs, discrep_sa_lhs
  implicit none
  real(dp), allocatable :: x(:,:)
  real(dp) :: before
  type(lhs_optimization_result) :: r1,r2,r3

  call lhs_design(10,3,x,randomized=.false.,seed=2_int64)
  before=discrepancy_value(x,'C2')
  call discrep_ese_lhs(x,r1,inner_iterations=20,candidates=12,outer_iterations=2,criterion='C2',seed=44_int64)
  if (discrepancy_value(r1%design,'C2')>before+1e-13_dp) error stop 'discrep ESE did not improve best design'

  before=phi_p(x,20.0_dp)
  call maximin_ese_lhs(x,r2,inner_iterations=20,candidates=12,outer_iterations=2,p=20.0_dp,seed=45_int64)
  if (phi_p(r2%design,20.0_dp)>before+1e-13_dp) error stop 'maximin ESE did not improve best design'

  call discrep_sa_lhs(x,r3,t0=0.02_dp,cooling=0.95_dp,iterations=40,criterion='W2',profile='GEOM',seed=46_int64)
  call check_same_columns(x,r3%design)
  if (r3%steps/=40) error stop 'SA step count'

  print *, 'test_lhs_opt: PASS'

contains
  subroutine check_same_columns(a,b)
    real(dp), intent(in) :: a(:,:),b(:,:)
    real(dp), allocatable :: aa(:),bb(:)
    integer :: j
    allocate(aa(size(a,1)),bb(size(a,1)))
    do j=1,size(a,2)
      aa=a(:,j);bb=b(:,j)
      call sort_vec(aa);call sort_vec(bb)
      if (maxval(abs(aa-bb))>1e-14_dp) error stop 'LHS column values changed'
    end do
  end subroutine check_same_columns
  subroutine sort_vec(v)
    real(dp), intent(inout) :: v(:)
    integer :: i,j
    real(dp)::t
    do i=2,size(v)
      t=v(i);j=i-1
      do while(j>=1)
        if(v(j)<=t) exit
        v(j+1)=v(j);j=j-1
      end do
      v(j+1)=t
    end do
  end subroutine sort_vec
end program test_lhs_opt
