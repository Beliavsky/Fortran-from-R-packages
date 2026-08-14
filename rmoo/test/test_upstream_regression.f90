program test_upstream_regression
  use rmoo
  implicit none
  type(rmoo_real_result) :: r2,r3,r1
  real(dp) :: sug(10,3),lo(3),hi(3)
  real(dp),allocatable :: ref(:,:)
  real(dp) :: x(10,3),y(1,3),d(10,1)

  sug=1.0_dp;lo=0.0_dp;hi=1.0_dp
  call rmoo_optimize_real(identity3,lo,hi,3,10,1,r2, &
    algorithm=ALG_NSGA2,pcrossover=0.0_dp,pmutation=0.0_dp, &
    suggestions=sug,seed=1)
  call assert_close(maxval(abs(r2%fitness-sug)),0.0_dp,1.0e-14_dp, &
    "upstream NSGA-II suggestions")

  call generate_reference_points(3,3,ref)
  call rmoo_optimize_real(identity3,lo,hi,3,10,1,r3, &
    algorithm=ALG_NSGA3,reference_dirs=ref,pcrossover=0.0_dp, &
    pmutation=0.0_dp,suggestions=sug,seed=1)
  call assert_true(all(shape(r3%fitness)==[10,3]),"upstream NSGA-III dimensions")

  call rmoo_optimize_real(identity3,lo,hi,3,10,1,r1, &
    algorithm=ALG_NSGA1,pcrossover=0.0_dp,pmutation=0.0_dp, &
    suggestions=sug,seed=1)
  call assert_close(maxval(abs(r1%fitness-sug)),0.0_dp,1.0e-14_dp, &
    "upstream NSGA-I suggestions")

  x=3.0_dp;y=1.0_dp
  call compute_perpendicular_distance(x,y,d)
  call assert_close(maxval(abs(d)),0.0_dp,1.0e-14_dp, &
    "upstream perpendicular-distance regression")

  call assert_close(generational_distance(ref,ref),0.0_dp,1.0e-14_dp, &
    "upstream GD regression")
  print *,"test_upstream_regression: PASS"
contains
  subroutine identity3(xv,f)
    real(dp),intent(in)::xv(:)
    real(dp),intent(out)::f(:)
    f=xv
  end subroutine identity3
  subroutine assert_true(cond,msg)
    logical,intent(in)::cond
    character(*),intent(in)::msg
    if(.not.cond)then
      write(*,*)"FAIL: ",trim(msg)
      error stop 1
    end if
  end subroutine
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol
    character(*),intent(in)::msg
    if(abs(a-b)>tol)then
      write(*,*)"FAIL: ",trim(msg),a,b
      error stop 1
    end if
  end subroutine
end program test_upstream_regression
