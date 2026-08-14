program test_criteria
  use dicedesign, only : dp, mindist, coverage, mesh_ratio, phi_p, discrepancy_all, mst_result, mst_criteria
  implicit none
  real(dp) :: x(5,2), disc(7)
  type(mst_result) :: mst

  x(1,:) = [0.1_dp,0.2_dp]
  x(2,:) = [0.7_dp,0.3_dp]
  x(3,:) = [0.4_dp,0.9_dp]
  x(4,:) = [0.9_dp,0.8_dp]
  x(5,:) = [0.2_dp,0.6_dp]

  call assert_close(mindist(x),0.36055512754639896_dp,1.0e-13_dp,'mindist')
  call assert_close(coverage(x),0.17123200398554297_dp,1.0e-13_dp,'coverage')
  call assert_close(mesh_ratio(x),1.4935759876113537_dp,1.0e-13_dp,'mesh_ratio')
  call assert_close(phi_p(x,7.0_dp),2.9730634402966736_dp,1.0e-12_dp,'phi_p')

  call discrepancy_all(x,disc)
  call assert_close(disc(1),0.14738083698741605_dp,1.0e-12_dp,'C2')
  call assert_close(disc(2),0.04658802898217998_dp,1.0e-12_dp,'L2')
  call assert_close(disc(3),0.09067034306271882_dp,1.0e-12_dp,'L2star')
  call assert_close(disc(4),0.15775860603395817_dp,1.0e-12_dp,'M2')
  call assert_close(disc(5),0.39438278078255196_dp,1.0e-12_dp,'S2')
  call assert_close(disc(6),0.18784627284623603_dp,1.0e-12_dp,'W2')
  call assert_close(disc(7),0.18964191402628985_dp,1.0e-12_dp,'Mix2')

  call mst_criteria(x,mst)
  call assert_close(mst%mean_length,0.45532103054522344_dp,1.0e-13_dp,'mst mean')
  call assert_close(mst%sd_length,0.08312848804299428_dp,1.0e-13_dp,'mst sd')
  if (size(mst%edges,2) /= 4) error stop 'mst edge count'

  print *, 'test_criteria: PASS'

contains
  subroutine assert_close(a,b,tol,name)
    real(dp), intent(in) :: a,b,tol
    character(len=*), intent(in) :: name
    if (abs(a-b)>tol) then
      print *, trim(name), a, b, abs(a-b)
      error stop 'assert_close failed'
    end if
  end subroutine assert_close
end program test_criteria
