! SPDX-License-Identifier: GPL-3.0-only
program test_multivariate_mds
  use mass
  use test_support
  implicit none
  real(dp) :: table(3,3), points0(4,2), distance(4,4), fitted_distance(4,4)
  real(dp), allocatable :: points(:,:), target(:), fitted(:), disparities(:)
  integer :: codes(8,2), levels(2), i, j, status
  type(correspondence_result) :: correspondence
  type(mca_result) :: mca
  type(mds_result) :: sammon_result, iso_result

  table=reshape([30.0_dp,5.0_dp,2.0_dp,4.0_dp,25.0_dp,3.0_dp, &
    1.0_dp,4.0_dp,28.0_dp],[3,3])
  call correspondence_analysis(table,2,correspondence)
  call assert_true(correspondence%status==mass_success,'correspondence status')
  call assert_true(all(correspondence%correlations>=0.0_dp),'correspondence correlations')

  levels=[2,2]
  codes=reshape([1,1,1,2,2,2,2,1, 1,1,2,1,2,2,1,2],[8,2])
  call mca_fit(codes,levels,2,mca)
  call assert_true(mca%status==mass_success,'MCA status')
  call assert_all_finite(reshape(mca%row_scores,[size(mca%row_scores)]),'MCA finite')

  points0=reshape([0.0_dp,1.0_dp,1.0_dp,0.0_dp, &
                   0.0_dp,0.0_dp,1.0_dp,1.0_dp],[4,2])
  do i=1,4
    do j=1,4
      distance(i,j)=sqrt(sum((points0(i,:)-points0(j,:))**2))
    end do
  end do
  call classical_mds(distance,2,points,status)
  call assert_true(status==mass_success,'classical MDS status')
  do i=1,4
    do j=1,4
      fitted_distance(i,j)=sqrt(sum((points(i,:)-points(j,:))**2))
    end do
  end do
  call assert_true(maxval(abs(fitted_distance-distance))<1.0e-7_dp,'classical MDS distances')

  call sammon(distance,sammon_result,dimensions=2,maxit=100)
  call assert_true(sammon_result%status==mass_success,'Sammon status')
  call assert_true(sammon_result%stress<1.0e-6_dp,'Sammon stress')
  call iso_mds(distance,iso_result,dimensions=2,maxit=100)
  call assert_true(iso_result%stress>=0.0_dp,'isoMDS finite stress')
  call shepard(distance,points,target,fitted,disparities,status)
  call assert_true(status==mass_success .and. size(target)==6,'Shepard output')
  write(*,'(a)') 'test_multivariate_mds: PASS'
end program test_multivariate_mds
