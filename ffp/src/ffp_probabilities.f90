module ffp_probabilities
  use ffp_kinds, only : dp
  implicit none
  private
  public :: normalize_probabilities, crisp_probabilities, exp_decay_probabilities
  public :: normal_kernel_probabilities, half_life, effective_scenarios, relative_entropy
contains
  subroutine normalize_probabilities(p, info)
    real(dp), intent(inout) :: p(:)
    integer, intent(out), optional :: info
    real(dp) :: s
    if (present(info)) info=0
    p=max(p,1.0e-300_dp); s=sum(p)
    if (.not.(s>0.0_dp)) then
      if (present(info)) info=1
      return
    end if
    p=p/s
  end subroutine normalize_probabilities

  subroutine crisp_probabilities(condition,p)
    logical, intent(in) :: condition(:)
    real(dp), intent(out) :: p(:)
    p=merge(1.0_dp,1.0e-30_dp,condition)
    call normalize_probabilities(p)
  end subroutine crisp_probabilities

  subroutine exp_decay_probabilities(n,lambda,p)
    integer, intent(in) :: n
    real(dp), intent(in) :: lambda
    real(dp), intent(out) :: p(n)
    integer :: i
    do i=1,n
      p(i)=exp(-lambda*real(n-i,dp))
    end do
    call normalize_probabilities(p)
  end subroutine exp_decay_probabilities

  subroutine normal_kernel_probabilities(x,mean,sigma,p,info)
    real(dp), intent(in) :: x(:,:),mean(:),sigma(:,:)
    real(dp), intent(out) :: p(:)
    integer, intent(out) :: info
    real(dp), allocatable :: a(:,:),d(:),z(:)
    real(dp) :: q
    integer :: i,n,k,stat
    n=size(x,1); k=size(x,2); info=0
    if (size(mean)/=k .or. size(sigma,1)/=k .or. size(sigma,2)/=k .or. size(p)/=n) then
      info=-1; return
    end if
    allocate(a(k,k),d(k),z(k)); a=sigma
    call cholesky_lower(a,stat)
    if (stat/=0) then; info=stat; return; end if
    do i=1,n
      d=x(i,:)-mean
      call forward_substitute(a,d,z)
      q=dot_product(z,z)
      p(i)=exp(-0.5_dp*min(q,1400.0_dp))
    end do
    call normalize_probabilities(p,info)
  contains
    subroutine cholesky_lower(a,info)
      real(dp),intent(inout)::a(:,:); integer,intent(out)::info
      integer::ii,jj,kk,nn; real(dp)::s
      nn=size(a,1); info=0
      do jj=1,nn
        do ii=jj,nn
          s=a(ii,jj)
          do kk=1,jj-1; s=s-a(ii,kk)*a(jj,kk); end do
          if (ii==jj) then
            if (s<=0.0_dp) then; info=jj; return; end if
            a(jj,jj)=sqrt(s)
          else
            a(ii,jj)=s/a(jj,jj)
          end if
        end do
        if (jj<nn) a(jj,jj+1:nn)=0.0_dp
      end do
    end subroutine
    subroutine forward_substitute(l,b,y)
      real(dp),intent(in)::l(:,:),b(:); real(dp),intent(out)::y(:)
      integer::ii
      do ii=1,size(b)
        y(ii)=(b(ii)-dot_product(l(ii,1:ii-1),y(1:ii-1)))/l(ii,ii)
      end do
    end subroutine
  end subroutine normal_kernel_probabilities

  pure real(dp) function half_life(lambda) result(days)
    real(dp),intent(in)::lambda
    days=anint(log(2.0_dp)/lambda)
  end function

  pure real(dp) function effective_scenarios(p) result(ens)
    real(dp),intent(in)::p(:)
    ens=exp(-sum(p*log(max(p,1.0e-300_dp))))
  end function

  pure real(dp) function relative_entropy(prior,posterior) result(value)
    real(dp),intent(in)::prior(:),posterior(:)
    value=sum(posterior*(log(max(posterior,1.0e-300_dp))-log(max(prior,1.0e-300_dp))))
  end function
end module ffp_probabilities
