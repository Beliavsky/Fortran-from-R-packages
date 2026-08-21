! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_counts
  use compositions_kinds, only: dp
  use compositions_geometry, only: closure
  use bayesm_rng, only: rand_poisson, rand_uniform
  implicit none
  private
  public :: rpois_composition, rmultinom_composition, count_totals, count_proportions
contains
  function rpois_composition(n,p,lambda) result(x)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(:),lambda
    integer :: x(n,size(p))
    real(dp) :: prob(size(p))
    integer :: i,j
    prob=closure(p)
    do i=1,n; do j=1,size(p); x(i,j)=rand_poisson(lambda*prob(j)); end do; end do
  end function rpois_composition

  function rmultinom_composition(n,p,total) result(x)
    integer, intent(in) :: n,total
    real(dp), intent(in) :: p(:)
    integer :: x(n,size(p))
    real(dp) :: prob(size(p)),u,cum
    integer :: i,j,k
    prob=closure(p); x=0
    do i=1,n
      do k=1,total
        u=rand_uniform(); cum=0.0_dp
        do j=1,size(p)
          cum=cum+prob(j)
          if(u<=cum) then; x(i,j)=x(i,j)+1; exit; end if
        end do
      end do
    end do
  end function rmultinom_composition

  function count_totals(x) result(t)
    integer, intent(in) :: x(:,:)
    integer :: t(size(x,1))
    t=sum(x,dim=2)
  end function count_totals

  function count_proportions(x) result(p)
    integer, intent(in) :: x(:,:)
    real(dp) :: p(size(x,1),size(x,2))
    integer :: i,t
    do i=1,size(x,1)
      t=sum(x(i,:)); if(t>0) then; p(i,:)=real(x(i,:),dp)/real(t,dp); else; p(i,:)=0.0_dp; end if
    end do
  end function count_proportions
end module compositions_counts
