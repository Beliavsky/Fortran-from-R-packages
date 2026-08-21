program test_spectral_parity
  use tsa, only : dp, spectral_estimate, spec_pgram, spec_ar, &
    modified_daniell_weights, spectral_kernel_df, spectral_kernel_bandwidth
  implicit none
  type(spectral_estimate) :: sp
  real(dp), allocatable :: w(:)
  real(dp) :: x(10), z(12), xm(10,2)
  integer :: status

  x = [1.0_dp,2.0_dp,4.0_dp,7.0_dp,11.0_dp,16.0_dp,22.0_dp,29.0_dp,37.0_dp,46.0_dp]
  sp = spec_pgram(x,spans=[3,5],taper=0.1_dp,pad=1,fast=.true.,demean=.false.,detrend=.true.)
  call check(sp%status == 0,'spec.pgram status')
  call check(sp%n_used == 20 .and. sp%orig_n == 10,'spec.pgram lengths')
  call close(sp%frequency(1),0.05_dp,1.0e-12_dp,'spec.pgram freq')
  call close(sp%spectrum(1),22.68626435_dp,2.0e-9_dp,'spec.pgram spec1')
  call close(sp%spectrum(2),21.23401364_dp,2.0e-9_dp,'spec.pgram spec2')
  call close(sp%degrees_freedom,4.680073126142595_dp,2.0e-11_dp,'spec.pgram df')
  call close(sp%bandwidth,0.07216878364870322_dp,2.0e-11_dp,'spec.pgram bandwidth')

  xm(:,1) = x
  xm(:,2) = [2.0_dp,1.0_dp,3.0_dp,5.0_dp,8.0_dp,13.0_dp,21.0_dp,34.0_dp,55.0_dp,89.0_dp]
  sp = spec_pgram(xm,spans=[3],taper=0.1_dp,pad=0,fast=.false.,demean=.false.,detrend=.true.)
  call check(sp%status == 0 .and. sp%n_series == 2,'multivariate spec.pgram status')
  call check(size(sp%spectrum_matrix,2) == 2 .and. size(sp%coherence,2) == 1, &
    'multivariate spectral shapes')
  call close(sp%spectrum_matrix(1,1),31.6020778_dp,3.0e-8_dp,'multivariate spec1')
  call close(sp%spectrum_matrix(1,2),379.642046_dp,3.0e-8_dp,'multivariate spec2')
  call close(sp%coherence(1,1),0.93698874_dp,3.0e-8_dp,'multivariate coherence')
  call close(sp%phase(1,1),0.20656071_dp,3.0e-8_dp,'multivariate phase')
  call close(sp%degrees_freedom,4.777574649603899_dp,2.0e-11_dp,'multivariate df')

  call modified_daniell_weights([3],w,status)
  call check(status == 0 .and. size(w) == 3,'kernel status')
  call close(w(1),0.25_dp,1.0e-14_dp,'kernel edge')
  call close(w(2),0.5_dp,1.0e-14_dp,'kernel center')
  call close(spectral_kernel_df(w),16.0_dp/3.0_dp,1.0e-14_dp,'kernel df')
  call close(spectral_kernel_bandwidth(w),sqrt(7.0_dp/12.0_dp),1.0e-14_dp,'kernel bandwidth')

  z = [1.0_dp,2.0_dp,1.0_dp,3.0_dp,5.0_dp,4.0_dp,6.0_dp,8.0_dp,7.0_dp,9.0_dp,10.0_dp,8.0_dp]
  sp = spec_ar(z,n_freq=5,order=2,demean=.true.,frequency=1.0_dp)
  call check(sp%status == 0 .and. sp%order == 2,'spec.ar status/order')
  call close(sp%frequency(2),0.125_dp,1.0e-14_dp,'spec.ar freq')
  call close(sp%spectrum(1),67.96245914770164_dp,2.0e-9_dp,'spec.ar zero')
  call close(sp%spectrum(2),13.03601415352587_dp,2.0e-9_dp,'spec.ar 0.125')
  call close(sp%spectrum(5),1.419642086945521_dp,2.0e-9_dp,'spec.ar nyquist')

  print '(a)', 'test_spectral_parity: PASS'
contains
  subroutine close(a,b,t,msg)
    real(dp), intent(in) :: a,b,t
    character(len=*), intent(in) :: msg
    if (abs(a-b) > t*max(1.0_dp,abs(b))) then
      print '(a,2es24.14)', trim(msg)//' FAIL: ',a,b
      error stop 1
    end if
  end subroutine close
  subroutine check(ok,msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg
    if (.not. ok) then
      print '(a)', trim(msg)//' FAIL'
      error stop 1
    end if
  end subroutine check
end program test_spectral_parity
