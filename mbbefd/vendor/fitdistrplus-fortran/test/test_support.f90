module test_support
  use fitdistrplus_kinds, only : dp
  implicit none
  private
  public :: assert_close, assert_true, laplace_logpdf_test, laplace_cdf_test
contains
  subroutine assert_close(actual,expected,tolerance,message)
    real(dp),intent(in)::actual,expected,tolerance
    character(len=*),intent(in)::message
    if(abs(actual-expected)>tolerance*max(1.0_dp,abs(expected)))then
      write(*,'(a,2(1x,es24.16))')trim(message),actual,expected
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition)then
      write(*,'(a)')trim(message)
      error stop 1
    end if
  end subroutine assert_true
  function laplace_logpdf_test(value,par) result(log_density)
    real(dp),intent(in)::value,par(:)
    real(dp)::log_density
    if(par(2)<=0.0_dp)then
      log_density=-huge(1.0_dp)
    else
      log_density=-log(2.0_dp*par(2))-abs(value-par(1))/par(2)
    end if
  end function laplace_logpdf_test

  function laplace_cdf_test(value,par) result(probability)
    real(dp),intent(in)::value,par(:)
    real(dp)::probability,z
    z=(value-par(1))/par(2)
    if(z<0.0_dp)then
      probability=0.5_dp*exp(z)
    else
      probability=1.0_dp-0.5_dp*exp(-z)
    end if
  end function laplace_cdf_test
end module test_support
