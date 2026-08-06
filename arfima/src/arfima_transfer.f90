module arfima_transfer
  use arfima_kinds, only : dp
  use arfima_status, only : arfima_ok, arfima_invalid_input
  use arfima_types, only : transfer_spec, arfima_error, set_error
  implicit none
  private
  public :: apply_transfer_function, apply_static_regression

contains

  subroutine apply_static_regression(y,x,beta,mean_value,residual,error)
    real(dp),intent(in)::y(:),x(:,:),beta(:),mean_value
    real(dp),allocatable,intent(out)::residual(:)
    type(arfima_error),intent(out)::error
    call set_error(error,arfima_ok,'')
    if(size(x,1)/=size(y) .or. size(x,2)/=size(beta)) then
      allocate(residual(0)); call set_error(error,arfima_invalid_input,'x and beta dimensions do not match y'); return
    end if
    allocate(residual(size(y)))
    residual=y-matmul(x,beta)-mean_value
  end subroutine apply_static_regression

  subroutine apply_transfer_function(y,transfer,mean_value,residual,effect,error)
    real(dp),intent(in)::y(:),mean_value
    type(transfer_spec),intent(in)::transfer
    real(dp),allocatable,intent(out)::residual(:)
    real(dp),allocatable,intent(out),optional::effect(:)
    type(arfima_error),intent(out)::error
    integer::n,nx,k,i,j,u,rr,ss,rk,sk,bk
    real(dp),allocatable::total(:),state(:)

    call set_error(error,arfima_ok,'')
    n=size(y)
    if(.not.allocated(transfer%x) .or. .not.allocated(transfer%r) .or. .not.allocated(transfer%s) .or. &
       .not.allocated(transfer%b) .or. .not.allocated(transfer%delta) .or. .not.allocated(transfer%omega)) then
      allocate(residual(0)); call set_error(error,arfima_invalid_input,'transfer specification is incomplete'); return
    end if
    nx=size(transfer%r)
    if(size(transfer%s)/=nx .or. size(transfer%b)/=nx .or. size(transfer%x,1)/=n .or. size(transfer%x,2)/=nx .or. &
       sum(transfer%r)/=size(transfer%delta) .or. sum(transfer%s)/=size(transfer%omega)) then
      allocate(residual(0)); call set_error(error,arfima_invalid_input,'transfer dimensions are inconsistent'); return
    end if
    if(any(transfer%r<0) .or. any(transfer%s<1) .or. any(transfer%b<0)) then
      allocate(residual(0)); call set_error(error,arfima_invalid_input,'transfer orders must satisfy r>=0, s>=1, b>=0'); return
    end if
    allocate(total(n)); total=0.0_dp
    rr=0; ss=0
    do k=1,nx
      rk=transfer%r(k); sk=transfer%s(k); bk=transfer%b(k)
      allocate(state(n)); state=-mean_value
      u=max(rk,sk+bk)-1
      do i=max(1,u+1),n
        do j=1,rk
          state(i)=state(i)+transfer%delta(rr+j)*state(i-j)
        end do
        if(i-bk>=1) state(i)=state(i)+transfer%omega(ss+1)*transfer%x(i-bk,k)
        do j=2,sk
          if(i-(j-1)-bk>=1) state(i)=state(i)-transfer%omega(ss+j)*transfer%x(i-(j-1)-bk,k)
        end do
      end do
      total=total+state
      rr=rr+rk; ss=ss+sk
      deallocate(state)
    end do
    allocate(residual(n)); residual=y-total
    if(present(effect)) then
      allocate(effect(n)); effect=total
    end if
  end subroutine apply_transfer_function

end module arfima_transfer
