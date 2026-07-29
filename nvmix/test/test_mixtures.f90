! SPDX-License-Identifier: GPL-3.0-or-later
program test_mixtures
  use nvmix
  implicit none
  type(nvmix_model) :: model
  type(sample_result) :: draw
  type(integration_control) :: ctrl
  real(dp), allocatable :: corr(:,:)
  real(dp) :: loc(2),scale(2,2),means(2),vars(2),sample_corr,v
  integer :: groups(2),families(2),n
  loc=0.0_dp; scale=reshape([1.0_dp,0.5_dp,0.5_dp,1.0_dp],[2,2])
  model=make_nvmix_model(loc,scale,mix_inverse_gamma,8.0_dp)
  draw=nvmix_random_sample(60000,model,12345_i8)
  call assert_true(draw%ok,'Student simulation')
  n=size(draw%x,1)
  means=sum(draw%x,dim=1)/real(n,dp)
  vars=sum((draw%x-spread(means,1,n))**2,dim=1)/real(n-1,dp)
  call assert_close(means(1),0.0_dp,0.03_dp,'sample mean')
  call assert_close(vars(1),8.0_dp/6.0_dp,0.05_dp,'sample variance')
  groups=[1,2]; families=mix_inverse_gamma
  model=make_grouped_model(loc,scale,groups,families,[6.0_dp,12.0_dp])
  corr=corgnvmix(model)
  draw=nvmix_random_sample(80000,model,2222_i8)
  means=sum(draw%x,dim=1)/real(size(draw%x,1),dp)
  sample_corr=sum((draw%x(:,1)-means(1))*(draw%x(:,2)-means(2)))
  sample_corr=sample_corr/sqrt(sum((draw%x(:,1)-means(1))**2)*sum((draw%x(:,2)-means(2))**2))
  call assert_close(sample_corr,corr(1,2),0.025_dp,'grouped correlation')
  ctrl%samples=4096
  v=nvmix_pdf([0.0_dp,0.0_dp],model,ctrl)
  call assert_true(v>0.0_dp .and. v<1.0_dp,'grouped density')
  call assert_close(dgammamix(3.0_dp,4,mix_constant,1.0_dp),chi_square_pdf(3.0_dp,4.0_dp),1.0e-13_dp,'gamma mix constant')
  call assert_close(pgammamix(4.0_dp,4,mix_inverse_gamma,9.0_dp),f_cdf(1.0_dp,4.0_dp,9.0_dp),1.0e-13_dp,'gamma mix t')
  print '(a)','test_mixtures: PASS'
contains
  subroutine assert_close(a,b,tol,label)
    real(dp), intent(in) :: a,b,tol
    character(*), intent(in) :: label
    if(abs(a-b)>tol*max(1.0_dp,abs(b)))then
      write(*,'(a,3es24.15)')trim(label)//' mismatch: ',a,b,abs(a-b); error stop 1
    end if
  end subroutine
  subroutine assert_true(ok,label)
    logical, intent(in) :: ok
    character(*), intent(in) :: label
    if(.not.ok)then; write(*,'(a)')trim(label)//' failed'; error stop 1; end if
  end subroutine
end program
