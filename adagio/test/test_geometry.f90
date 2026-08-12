program test_geometry
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
  use adagio
  implicit none
  type(maxsub_result) :: r1
  type(maxsub2d_result) :: r2
  type(maxempty_result) :: re
  type(count_result) :: rc
  type(hamiltonian_result) :: rh
  real(dp) :: a(4,4), z(4), lo(4), up(4), x(4), xr(4)
  integer, allocatable :: io(:)

  r1 = maxsub([-2._dp, 3._dp, 4._dp, -10._dp, 5._dp])
  call check(abs(r1%sum-7._dp)<1e-12_dp .and. r1%first==2 .and. r1%last==3, 'maxsub')

  a = transpose(reshape([0._dp,-2._dp,-7._dp,0._dp, &
                         9._dp,2._dp,-6._dp,2._dp, &
                        -4._dp,1._dp,-4._dp,1._dp, &
                        -1._dp,8._dp,0._dp,2._dp], [4,4]))
  r2 = maxsub2d(a)
  call check(abs(r2%sum-15._dp)<1e-12_dp, 'maxsub2d sum')
  call check(all(r2%inds == [2,4,1,2]), 'maxsub2d indices')

  re = maxempty([0.25_dp,0.75_dp], [0.5_dp,0.5_dp])
  call check(abs(re%area-0.5_dp)<1e-12_dp, 'maxempty')

  rc = count_values([3._dp,1._dp,3._dp,2._dp,1._dp])
  call check(maxval(abs(rc%values-[1._dp,2._dp,3._dp])) < 1e-15_dp, 'count values')
  call check(all(rc%counts == [2,1,2]), 'count counts')

  io = occurs([2._dp,3._dp], [1._dp,2._dp,3._dp,2._dp,3._dp,4._dp])
  call check(all(io == [2,4]), 'occurs')

  lo = [-2._dp, ieee_value(0._dp,ieee_negative_inf), 0._dp, ieee_value(0._dp,ieee_negative_inf)]
  up = [ 3._dp, ieee_value(0._dp,ieee_positive_inf), ieee_value(0._dp,ieee_positive_inf), 4._dp]
  x  = [ 0.2_dp, -1.5_dp, 2.0_dp, 1.0_dp]
  z = transfinite_forward(x, lo, up)
  xr = transfinite_inverse(z, lo, up)
  call check(maxval(abs(xr-x)) < 1e-11_dp, 'transfinite round trip')

  rh = hamiltonian([1,2, 2,3, 3,4, 4,1], start=1, cycle=.true.)
  call check(rh%found .and. size(rh%path)==4, 'hamiltonian cycle')

  print *, 'test_geometry: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print *, 'FAIL: ',trim(msg)
      error stop 1
    end if
  end subroutine
end program test_geometry
