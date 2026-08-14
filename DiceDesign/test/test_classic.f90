program test_classic
  use iso_fortran_env, only : int64
  use dicedesign, only : dp, fact_design, lhs_design, olh_design, nolh_design, nolhdr_design, &
    runif_faure, faureprime_design, scale_design, unscale_design
  implicit none
  real(dp), allocatable :: x(:,:), y(:,:), z(:,:)
  real(dp) :: expected_olh(9,3), expected_faure(5,3)
  integer :: p

  call fact_design(2,[2,3],x)
  if (size(x,1)/=6 .or. size(x,2)/=2) error stop 'fact_design shape'
  call assert_close(x(1,1),0.0_dp,1e-15_dp,'fact 11')
  call assert_close(x(4,1),1.0_dp,1e-15_dp,'fact 41')
  call assert_close(x(2,2),0.5_dp,1e-15_dp,'fact 22')

  call lhs_design(8,3,y,randomized=.false.,seed=7_int64)
  call check_centered_lhs(y)

  expected_olh(1,:)=[0.625_dp,0.750_dp,0.875_dp]
  expected_olh(2,:)=[0.750_dp,0.375_dp,0.000_dp]
  expected_olh(3,:)=[0.875_dp,1.000_dp,0.375_dp]
  expected_olh(4,:)=[1.000_dp,0.125_dp,0.750_dp]
  expected_olh(5,:)=[0.500_dp,0.500_dp,0.500_dp]
  expected_olh(6,:)=[0.375_dp,0.250_dp,0.125_dp]
  expected_olh(7,:)=[0.250_dp,0.625_dp,1.000_dp]
  expected_olh(8,:)=[0.125_dp,0.000_dp,0.625_dp]
  expected_olh(9,:)=[0.000_dp,0.875_dp,0.250_dp]
  call olh_design(3,x)
  call assert_matrix(x,expected_olh,1e-14_dp,'olh')

  call nolh_design(7,x)
  call assert_close(x(1,1),0.3125_dp,1e-15_dp,'nolh 11')
  call assert_close(x(1,2),1.0_dp,1e-15_dp,'nolh 12')
  call assert_close(x(1,7),0.5625_dp,1e-15_dp,'nolh 17')
  if (size(x,1)/=17) error stop 'nolh row count'

  call nolhdr_design(8,x)
  call assert_close(x(1,1),0.0625_dp,1e-15_dp,'nolhdr 11')
  call assert_close(x(1,6),0.96875_dp,1e-15_dp,'nolhdr 16')
  if (size(x,1)/=33) error stop 'nolhdr row count'

  expected_faure(1,:)=[1.0_dp/3.0_dp,1.0_dp/3.0_dp,1.0_dp/3.0_dp]
  expected_faure(2,:)=[2.0_dp/3.0_dp,2.0_dp/3.0_dp,2.0_dp/3.0_dp]
  expected_faure(3,:)=[1.0_dp/9.0_dp,4.0_dp/9.0_dp,7.0_dp/9.0_dp]
  expected_faure(4,:)=[4.0_dp/9.0_dp,7.0_dp/9.0_dp,1.0_dp/9.0_dp]
  expected_faure(5,:)=[7.0_dp/9.0_dp,1.0_dp/9.0_dp,4.0_dp/9.0_dp]
  call runif_faure(5,3,x)
  call assert_matrix(x,expected_faure,2e-15_dp,'faure')

  call faureprime_design(3,2,x,p)
  if (p/=3 .or. size(x,1)/=8) error stop 'faureprime metadata'
  if (minval(x)<=0.0_dp .or. maxval(x)>=1.0_dp) error stop 'faureprime default open range'

  allocate(z(3,2))
  z(1,:)=[2.0_dp,10.0_dp]
  z(2,:)=[4.0_dp,20.0_dp]
  z(3,:)=[6.0_dp,30.0_dp]
  call scale_design(z,x,lower=[2.0_dp,10.0_dp],upper=[6.0_dp,30.0_dp])
  call unscale_design(x,y,lower=[2.0_dp,10.0_dp],upper=[6.0_dp,30.0_dp])
  call assert_matrix(y,z,1e-14_dp,'scale roundtrip')

  print *, 'test_classic: PASS'

contains
  subroutine assert_close(a,b,tol,name)
    real(dp), intent(in) :: a,b,tol
    character(len=*), intent(in) :: name
    if (abs(a-b)>tol) then
      print *, trim(name), a, b, abs(a-b)
      error stop 'assert_close failed'
    end if
  end subroutine assert_close

  subroutine assert_matrix(a,b,tol,name)
    real(dp), intent(in) :: a(:,:),b(:,:),tol
    character(len=*), intent(in) :: name
    if (any(shape(a)/=shape(b)) .or. maxval(abs(a-b))>tol) then
      print *, trim(name), maxval(abs(a-b))
      error stop 'assert_matrix failed'
    end if
  end subroutine assert_matrix

  subroutine check_centered_lhs(a)
    real(dp), intent(in) :: a(:,:)
    integer :: n,j,i,k
    logical, allocatable :: seen(:)
    n=size(a,1)
    allocate(seen(n))
    do j=1,size(a,2)
      seen=.false.
      do i=1,n
        k=nint(a(i,j)*real(n,dp)+0.5_dp)
        if (k<1 .or. k>n) error stop 'lhs stratum out of range'
        if (seen(k)) error stop 'lhs repeated stratum'
        seen(k)=.true.
      end do
    end do
  end subroutine check_centered_lhs
end program test_classic
