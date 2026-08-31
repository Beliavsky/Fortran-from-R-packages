program test_fft_spatial_process
use fields, only: dp,circulant_embedding_2d,circulant_setup_2d,circulant_sample_2d, &
                  spatial_profile_grid_result,spatial_profile_grid
implicit none
type(circulant_embedding_2d)::emb
type(spatial_profile_grid_result)::grid
real(dp),allocatable::draws(:,:,:)
real(dp)::x(8,1),y(8),ranges(3),lams(3)
integer::i,info
call circulant_setup_2d(8,8,0.1_dp,0.1_dp,emb,model='exponential',a_range=0.35_dp,info=info)
call check(info==0,'circulant setup')
draws=circulant_sample_2d(emb,4,info)
call check(info==0 .and. all(draws==draws),'circulant draw')
call check(abs(sum(draws)/real(size(draws),dp))<1.5_dp,'draw mean sanity')
do i=1,8;x(i,1)=real(i-1,dp)/7;y(i)=sin(4*x(i,1));end do
ranges=[0.15_dp,0.4_dp,1.0_dp];lams=[0.01_dp,0.1_dp,1.0_dp]
grid=spatial_profile_grid(x,y,ranges,lams,model='exponential',criterion='reml')
call check(grid%best_range>=1 .and. grid%best_lambda>=1,'profile grid best')
call check(all(grid%objective<huge(1.0_dp)),'profile finite')
print *,'test_fft_spatial_process: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok;character(len=*),intent(in)::msg
if(.not.ok)then;print *,'FAIL: ',msg;error stop 1;end if
end subroutine
end program
