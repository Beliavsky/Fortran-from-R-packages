program test_outlier_calibration
  use compositions
  implicit none
  real(dp) :: x(18,3),qe,qm
  real(dp), allocatable :: d(:)
  logical, allocatable :: flag(:)
  type(outlier_classification_result) :: c
  integer :: i
  do i=1,17
    x(i,1)=0.30_dp+0.002_dp*real(mod(i,5)-2,dp)
    x(i,2)=0.30_dp+0.002_dp*real(mod(2*i,7)-3,dp)
    x(i,3)=1.0_dp-x(i,1)-x(i,2)
  end do
  x(18,:)=[0.90_dp,0.05_dp,0.05_dp]
  d=acomp_mahalanobis(x,robust=.false.)
  if(d(18)<=maxval(d(1:17))) error stop 'outlier distance ordering'
  qe=q_empirical_mahalanobis(0.95_dp,18,2,49,robust=.false.,seed=99)
  qm=q_max_mahalanobis(0.95_dp,18,2,49,robust=.false.,seed=99)
  if(qm<qe) error stop 'max calibration should exceed typical calibration'
  flag=is_mahalanobis_outlier(x,alpha=0.05_dp,replicates=79,corrected=.true.,robust=.false.,seed=77)
  if(.not.flag(18)) error stop 'extreme composition not flagged'
  c=outlier_classifier_best(x,alpha=0.05_dp,replicates=49,robust=.false.,seed=88)
  if(.not.c%is_outlier(18)) error stop 'classifier failed'
  print *, 'test_outlier_calibration: PASS'
end program test_outlier_calibration
