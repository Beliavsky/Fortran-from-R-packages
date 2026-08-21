! SPDX-License-Identifier: GPL-2.0-or-later
module mnb_simulation
  use mnb_kinds, only : dp
  use mnb_math, only : gamma_rng, poisson_rng
  implicit none
  private
  public :: simulate_mnb
contains
  subroutine simulate_mnb(n,mi,x,par,y,offset,random_effect)
    integer,intent(in)::n,mi
    real(dp),intent(in)::x(:,:),par(:)
    integer,intent(out)::y(:)
    real(dp),intent(in),optional::offset(:)
    real(dp),intent(out),optional::random_effect(:)
    real(dp),allocatable::eta(:),u(:)
    integer::i,j,k
    if(size(y)/=n*mi)error stop 'simulate_mnb: wrong y size'
    if(par(1)<=0.0_dp)error stop 'simulate_mnb: phi must be positive'
    allocate(eta(n*mi),u(n));eta=matmul(x,par(2:));if(present(offset))eta=eta+offset
    do i=1,n
      u(i)=gamma_rng(par(1))/par(1)
      do j=1,mi
        k=(i-1)*mi+j;y(k)=poisson_rng(exp(eta(k))*u(i))
      end do
    end do
    if(present(random_effect))random_effect=u
  end subroutine simulate_mnb
end module mnb_simulation
