program test_pareto_reference
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rmoo
  implicit none
  real(dp) :: f(6,2), cd(6), gdv
  integer :: rank(6)
  real(dp), allocatable :: ref(:,:)
  integer, allocatable :: comp(:,:)
  real(dp) :: ideal(2), worst(2), wf(2), wp(2), nadir(2)
  real(dp) :: extreme(2,2), smin(2)
  real(dp) :: x(2,2), y(2,2), d(2,2)

  f = reshape([ &
    1.0_dp,4.0_dp, &
    2.0_dp,3.0_dp, &
    3.0_dp,2.0_dp, &
    4.0_dp,1.0_dp, &
    3.0_dp,4.0_dp, &
    5.0_dp,5.0_dp], shape(f), order=[2,1])
  call non_dominated_sort(f,rank)
  call assert_true(all(rank == [1,1,1,1,2,3]), "non-dominated ranks")

  call crowding_distance(f,rank,cd)
  call assert_true(.not.ieee_is_finite(cd(1)), "crowding endpoint 1")
  call assert_true(.not.ieee_is_finite(cd(4)), "crowding endpoint 4")
  call assert_true(cd(2) > 0.0_dp .and. cd(3) > 0.0_dp, "crowding interior")

  call get_fixed_rowsum_integer_matrix(3,2,comp)
  call assert_true(size(comp,1)==6 .and. size(comp,2)==3, "composition shape")
  call assert_true(all(sum(comp,dim=2)==2), "composition sums")
  call generate_reference_points(3,2,ref)
  call assert_true(size(ref,1)==6, "reference point count")
  call assert_close(maxval(abs(sum(ref,dim=2)-1.0_dp)),0.0_dp,1.0e-14_dp, &
    "reference hyperplane")

  gdv = generational_distance(f(1:4,:),f(1:4,:))
  call assert_close(gdv,0.0_dp,1.0e-14_dp,"GD identity")
  call assert_close(igd(f(1:4,:),f(1:4,:)),0.0_dp,1.0e-14_dp,"IGD identity")

  ideal = 0.0_dp
  worst = 2.0_dp
  wf = 1.0_dp
  wp = 2.0_dp
  extreme = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])
  call get_nadir_point(extreme,ideal,worst,wf,wp,nadir)
  call assert_close(maxval(abs(nadir-1.0_dp)),0.0_dp,1.0e-12_dp,"nadir")

  smin = huge(1.0_dp)
  extreme = 0.0_dp
  call perform_scalarizing(reshape([1.0_dp,4.0_dp,4.0_dp,1.0_dp],[2,2],order=[2,1]), &
    [0.0_dp,0.0_dp],smin,extreme)
  call assert_true(all(ieee_is_finite(smin)),"scalarizing finite")

  x = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2],order=[2,1])
  y = x
  call compute_perpendicular_distance(x,y,d)
  call assert_close(d(1,1),0.0_dp,1.0e-14_dp,"perpendicular diagonal")
  call assert_close(d(1,2),1.0_dp,1.0e-14_dp,"perpendicular orthogonal")

  print *, "test_pareto_reference: PASS"
contains
  subroutine assert_true(cond,msg)
    logical,intent(in)::cond
    character(*),intent(in)::msg
    if(.not.cond)then
      write(*,*) "FAIL: ",trim(msg)
      error stop 1
    end if
  end subroutine
  subroutine assert_close(xv,yv,tol,msg)
    real(dp),intent(in)::xv,yv,tol
    character(*),intent(in)::msg
    if(abs(xv-yv)>tol)then
      write(*,*) "FAIL: ",trim(msg),xv,yv
      error stop 1
    end if
  end subroutine
end program test_pareto_reference
