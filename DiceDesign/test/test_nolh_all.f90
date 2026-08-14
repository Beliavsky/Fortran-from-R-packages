program test_nolh_all
  use dicedesign, only : dp, nolh_design, nolhdr_design
  implicit none
  real(dp), allocatable :: x(:,:)
  integer :: d, j, expected_n

  do d=2,29
    call nolh_design(d,x)
    if (d<=7) then
      expected_n=17
    else if (d<=11) then
      expected_n=33
    else if (d<=16) then
      expected_n=65
    else if (d<=22) then
      expected_n=129
    else
      expected_n=257
    end if
    if (size(x,1)/=expected_n .or. size(x,2)/=d) error stop 'nolh shape'
    do j=1,d
      call check_latin_column(x(:,j))
    end do
  end do

  do d=2,29
    call nolhdr_design(d,x)
    if (d<=7) then
      expected_n=17
    else if (d<=11) then
      expected_n=33
    else if (d<=16) then
      expected_n=65
    else if (d<=22) then
      expected_n=129
    else
      expected_n=257
    end if
    if (size(x,1)/=expected_n .or. size(x,2)/=d) error stop 'nolhdr shape'
    do j=1,d
      call check_latin_column(x(:,j))
    end do
  end do
  print *, 'test_nolh_all: PASS'
contains
  subroutine check_latin_column(v)
    real(dp), intent(in) :: v(:)
    real(dp), allocatable :: s(:)
    integer :: k,n
    n=size(v)
    allocate(s(n));s=v
    call sort_vec(s)
    do k=1,n
      if (abs(s(k)-real(k-1,dp)/real(n-1,dp))>2e-14_dp) error stop 'not Latin levels'
    end do
  end subroutine check_latin_column
  subroutine sort_vec(v)
    real(dp), intent(inout)::v(:)
    integer::i,j
    real(dp)::key
    do i=2,size(v)
      key=v(i);j=i-1
      do while(j>=1)
        if(v(j)<=key) exit
        v(j+1)=v(j);j=j-1
      end do
      v(j+1)=key
    end do
  end subroutine sort_vec
end program test_nolh_all
