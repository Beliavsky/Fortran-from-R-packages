! SPDX-License-Identifier: GPL-3.0-only
program test_statistics
    use rsdc, only: dp, rsdc_model, rsdc_control, rsdc_bands_result, rsdc_bootstrap_result
    use rsdc, only: rsdc_estimate, rsdc_corr_bands, rsdc_parametric_bootstrap
    use rsdc, only: rsdc_scores, rsdc_robust_vcov, rsdc_const
    implicit none
    integer, parameter :: t=90,k=2
    real(dp) :: y(t,k)
    real(dp), allocatable :: scores(:,:), v(:,:)
    type(rsdc_model) :: fit
    type(rsdc_control) :: ctl, boot_ctl
    type(rsdc_bands_result) :: bands
    type(rsdc_bootstrap_result) :: boot
    logical :: ok
    integer :: i
    do i=1,t
        y(i,1)=sin(0.19_dp*i)+0.2_dp*cos(0.07_dp*i)
        y(i,2)=0.45_dp*y(i,1)+0.8_dp*cos(0.11_dp*i)
    end do
    y(:,1)=y(:,1)/sqrt(sum(y(:,1)**2)/real(t,dp))
    y(:,2)=y(:,2)/sqrt(sum(y(:,2)**2)/real(t,dp))
    ctl%population_size=24; ctl%max_global_iterations=45; ctl%max_local_iterations=35
    ctl%compute_vcov=.true.; ctl%seed=22
    call rsdc_estimate(rsdc_const,y,1,fit,control=ctl,ok=ok)
    call check(ok .and. allocated(fit%vcov),'vcov estimation')
    call rsdc_scores(fit,y,scores,ok=ok)
    call check(ok .and. all(shape(scores)==[t,1]),'scores')
    call rsdc_robust_vcov(fit,y,'opg',v,ok=ok)
    call check(ok .and. v(1,1)>0.0_dp,'opg')
    call rsdc_corr_bands(fit,y,20,bands,level=0.9_dp,seed=4,ok=ok)
    call check(ok .and. bands%n_used>=2,'bands')
    boot_ctl=ctl; boot_ctl%compute_vcov=.false.; boot_ctl%max_global_iterations=18
    boot_ctl%max_local_iterations=15; boot_ctl%population_size=12
    call rsdc_parametric_bootstrap(fit,y,3,boot,seed=8,control=boot_ctl,ok=ok)
    call check(ok .and. boot%n_success>=2 .and. allocated(boot%covariance),'bootstrap')
    call check(all(shape(boot%confidence_interval)==[1,2]),'bootstrap intervals')
    print '(a)', 'test_statistics: PASS'
contains
    subroutine check(condition,message)
        logical,intent(in)::condition
        character(len=*),intent(in)::message
        if(.not.condition) error stop message
    end subroutine check
end program test_statistics
