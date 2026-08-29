program linear_gmm_example
use gmm
implicit none
real(dp) :: y(8), x(8,2), z(8,3)
type(linear_gmm_result_t) :: fit
integer :: i

y=[1.2_dp,1.9_dp,2.7_dp,3.8_dp,4.1_dp,5.2_dp,5.9_dp,7.1_dp]
do i=1,8
   x(i,:)=[1.0_dp,real(i-1,dp)]
   z(i,:)=[1.0_dp,real(i-1,dp),merge(1.0_dp,-1.0_dp,mod(i,2)==1)]
end do
call linear_gmm_fit(y,x,z,fit,method=LINEAR_TWO_STEP,covariance=COV_MDS)
print '(a,2f14.8)', 'coefficients: ',fit%coefficients
print '(a,f14.8)', 'J statistic:  ',fit%j_stat
print '(a,f14.8)', 'J p-value:    ',fit%j_pvalue
end program linear_gmm_example
