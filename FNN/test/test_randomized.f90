program test_randomized
  use fnn
  implicit none
  integer :: rep,n,p,k,m,nseed,i
  integer, allocatable :: seed(:)
  real(dp), allocatable, target :: x(:,:),q(:,:)
  type(knn_result) :: ref,kd,cover
  real(dp), parameter :: tol=5.0e-12_dp

  call random_seed(size=nseed); allocate(seed(nseed)); seed=24681357; call random_seed(put=seed)
  do rep=1,200
    n=8+mod(7*rep,33); p=1+mod(5*rep,7); k=1+mod(3*rep,min(8,n-1)); m=1+mod(11*rep,9)
    allocate(x(n,p),q(m,p)); call random_number(x); call random_number(q)
    do i=1,n
      x(i,:)=x(i,:)+real(i,dp)*1.0e-8_dp
    end do
    ref=get_knn(x,k,"brute"); kd=get_knn(x,k,"kd_tree"); cover=get_knn(x,k,"cover_tree")
    call check_equal(ref,kd,"self kd")
    call check_equal(ref,cover,"self cover")
    ref=get_knnx(x,q,k,"brute"); kd=get_knnx(x,q,k,"kd_tree"); cover=get_knnx(x,q,k,"cover_tree")
    call check_equal(ref,kd,"query kd")
    call check_equal(ref,cover,"query cover")
    deallocate(x,q)
  end do
  print *, "test_randomized: PASS (200 data sets)"
contains
  subroutine check_equal(a,b,label)
    type(knn_result), intent(in) :: a,b
    character(len=*), intent(in) :: label
    if(.not.all(a%index==b%index)) then
      print *, "index mismatch: ",trim(label)
      error stop 1
    end if
    if(maxval(abs(a%distance-b%distance))>tol) then
      print *, "distance mismatch: ",trim(label),maxval(abs(a%distance-b%distance))
      error stop 1
    end if
  end subroutine check_equal
end program test_randomized
