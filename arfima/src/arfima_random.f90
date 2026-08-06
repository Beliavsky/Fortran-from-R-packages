module arfima_random
  use arfima_kinds, only : dp, pi_dp
  implicit none
  private
  public :: set_random_seed, random_normal_vector
contains
  subroutine set_random_seed(seed)
    integer,intent(in)::seed
    integer::n,i
    integer,allocatable::put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i=1,n
      put(i)=mod(abs(seed)+104729*i,2147483646)+1
    end do
    call random_seed(put=put)
  end subroutine set_random_seed

  subroutine random_normal_vector(x)
    real(dp),intent(out)::x(:)
    integer::i
    real(dp)::u1,u2,r
    i=1
    do while(i<=size(x))
      call random_number(u1); call random_number(u2)
      u1=max(u1,tiny(1.0_dp))
      r=sqrt(-2.0_dp*log(u1))
      x(i)=r*cos(2.0_dp*pi_dp*u2)
      if(i+1<=size(x)) x(i+1)=r*sin(2.0_dp*pi_dp*u2)
      i=i+2
    end do
  end subroutine random_normal_vector
end module arfima_random
