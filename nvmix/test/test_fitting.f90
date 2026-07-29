! SPDX-License-Identifier: GPL-3.0-or-later
program test_fitting
  use nvmix
  implicit none
  type(sample_result) :: draw
  type(fit_result) :: fit
  type(qq_result) :: qq
  type(nvmix_model) :: model
  real(dp) :: loc(2),scale(2,2)
  loc=[0.4_dp,-0.2_dp]; scale=reshape([1.0_dp,0.25_dp,0.25_dp,0.7_dp],[2,2])
  draw=rnorm_mv(5000,loc,scale,2026_i8)
  fit=fit_norm(draw%x)
  call assert_true(fit%ok,'normal fit')
  call assert_close(fit%loc(1),loc(1),0.04_dp,'normal location')
  call assert_close(fit%scale(1,1),scale(1,1),0.05_dp,'normal scale')
  draw=rstudent_mv(3500,8.0_dp,loc,scale,9281_i8)
  fit=fit_student(draw%x,max_iterations=40)
  call assert_true(fit%ok,'Student fit')
  call assert_true(fit%mixing_parameter(1)>4.0_dp .and. fit%mixing_parameter(1)<20.0_dp,'Student df')
  model=make_nvmix_model(fit%loc,fit%scale,mix_inverse_gamma,fit%mixing_parameter(1))
  qq=qqplot_maha(draw%x(1:200,:),model)
  call assert_true(qq%ok .and. all(qq%observed(2:)>=qq%observed(:199)),'QQ data')
  print '(a)','test_fitting: PASS'
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
