program test_multivariate
  use extra_distr
  use test_support
  implicit none
  integer :: i,j
  real(dp) :: s, expected
  real(dp) :: x2(2),a2(2),prob2(2),means(2),sds(2),weights(2)
  integer :: xi2(2), nc2(2)
  real(dp), allocatable :: rd(:,:), rb(:,:), rm(:)
  integer, allocatable :: ri(:,:)

  call assert_close(dbvnorm(0.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.0_dp), &
       1.0_dp/(2.0_dp*acos(-1.0_dp)),1.0e-12_dp,'bivariate normal independent density')

  s=0.0_dp
  do i=0,25
    do j=0,25
      s=s+dbvpois(i,j,1.2_dp,1.7_dp,0.4_dp)
    end do
  end do
  call assert_close(s,1.0_dp,2.0e-10_dp,'bivariate Poisson normalization')

  x2=[0.3_dp,0.7_dp]; a2=[1.0_dp,1.0_dp]
  call assert_close(ddirichlet(x2,a2),1.0_dp,1.0e-12_dp,'uniform Dirichlet density')
  xi2=[1,1]
  call assert_close(ddirmnom(xi2,2,a2),1.0_dp/3.0_dp,1.0e-12_dp,'Dirichlet-multinomial fixture')

  means=[-1.0_dp,2.0_dp]; sds=[1.0_dp,0.5_dp]; weights=[0.25_dp,0.75_dp]
  expected=weights(1)*exp(-0.5_dp)/(sqrt(2.0_dp*acos(-1.0_dp))) + &
           weights(2)*exp(-8.0_dp)/(0.5_dp*sqrt(2.0_dp*acos(-1.0_dp)))
  call assert_close(dmixnorm(0.0_dp,means,sds,weights),expected,1.0e-12_dp,'normal mixture density')
  call assert_close(pmixnorm(0.0_dp,means,sds,weights,lower_tail=.false.), &
       1.0_dp-pmixnorm(0.0_dp,means,sds,weights),1.0e-12_dp,'normal mixture upper tail')

  means=[1.5_dp,5.0_dp]; weights=[0.4_dp,0.6_dp]
  s=0.0_dp; do i=0,100; s=s+dmixpois(i,means,weights); end do
  call assert_close(s,1.0_dp,1.0e-12_dp,'Poisson mixture normalization')

  prob2=[0.25_dp,0.75_dp]
  s=0.0_dp
  do i=0,2
    xi2=[i,2-i]
    s=s+dmnom(xi2,2,prob2)
  end do
  call assert_close(s,1.0_dp,1.0e-12_dp,'multinomial normalization')
  xi2=[1,1]
  call assert_close(dmnom(xi2,2,prob2),0.375_dp,1.0e-12_dp,'multinomial fixture')

  nc2=[2,3]
  s=0.0_dp
  do i=0,2
    xi2=[i,2-i]
    s=s+dmvhyper(xi2,nc2,2)
  end do
  call assert_close(s,1.0_dp,1.0e-12_dp,'multivariate hypergeometric normalization')

  call seed_rng(20260804)
  rd=rdirichlet(200,a2)
  call assert_true(all(rd>=0.0_dp),'Dirichlet random support')
  do i=1,size(rd,1)
    call assert_close(sum(rd(i,:)),1.0_dp,2.0e-14_dp,'Dirichlet random simplex')
  end do

  rb=rbvnorm(4000,0.0_dp,0.0_dp,1.0_dp,2.0_dp,0.5_dp)
  call assert_true(abs(sum(rb(:,1))/real(size(rb,1),dp))<0.08_dp,'bivariate normal random mean')
  call assert_true(abs(sum(rb(:,1)*rb(:,2))/real(size(rb,1),dp)-1.0_dp)<0.15_dp,'bivariate normal random covariance')

  ri=rmnom(100,7,prob2)
  call assert_true(all(sum(ri,dim=2)==7),'multinomial random totals')
  ri=rmvhyper(100,nc2,2)
  call assert_true(all(sum(ri,dim=2)==2),'multivariate hypergeometric random totals')

  rm=rmixnorm(20,[-1.0_dp,1.0_dp],[1.0_dp,1.0_dp],[0.5_dp,0.5_dp])
  call assert_int(size(rm),20,'mixture random length')

  call finish_tests('multivariate and mixtures')
end program test_multivariate
