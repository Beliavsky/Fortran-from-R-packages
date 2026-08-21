module relsurv_data
  use relsurv_kinds, only : dp
  implicit none
  private
  type, public :: split_data_type
    integer, allocatable :: source(:), episode(:), event(:)
    real(dp), allocatable :: start(:), stop(:)
  end type split_data_type
  public :: survsplit_counting, inverse_time_monotone
contains
  subroutine survsplit_counting(start,stop,event,cut,out,zero)
    real(dp),intent(in),optional::start(:)
    real(dp),intent(in)::stop(:),cut(:)
    integer,intent(in)::event(:)
    type(split_data_type),intent(out)::out
    real(dp),intent(in),optional::zero
    real(dp)::z,a,b
    integer::n,i,j,k,nout
    logical::keep
    n=size(stop); z=0.0_dp; if(present(zero))z=zero
    nout=0
    do i=1,n
      a=z; if(present(start))a=start(i)
      do j=0,size(cut)
        if(j==0) then; b=merge(cut(1),stop(i),size(cut)>0)
        else if(j<size(cut)) then; b=cut(j+1)
        else; b=stop(i)
        end if
        b=min(b,stop(i))
        keep=a<b
        if(keep) nout=nout+1
        if(stop(i)<=b) exit
        a=b
      end do
    end do
    allocate(out%source(nout),out%episode(nout),out%event(nout),out%start(nout),out%stop(nout))
    k=0
    do i=1,n
      a=z; if(present(start))a=start(i)
      do j=0,size(cut)
        if(size(cut)==0) then
          b=stop(i)
        else if(j<size(cut)) then
          b=cut(j+1)
        else
          b=stop(i)
        end if
        b=min(b,stop(i)); keep=a<b
        if(keep) then
          k=k+1; out%source(k)=i; out%episode(k)=j
          out%start(k)=a; out%stop(k)=b
          out%event(k)=merge(event(i),0,stop(i)<=b+10.0_dp*epsilon(1.0_dp))
        end if
        if(stop(i)<=b) exit
        a=b
      end do
    end do
  end subroutine survsplit_counting

  function inverse_time_monotone(target,t,prob) result(x)
    real(dp),intent(in)::target,t(:),prob(:)
    real(dp)::x,w
    integer::i
    if(target<=prob(1)) then; x=t(1); return; end if
    if(target>=prob(size(prob))) then; x=t(size(t)); return; end if
    do i=2,size(prob)
      if(target<=prob(i)) then
        if(prob(i)==prob(i-1)) then
          x=t(i)
        else
          w=(target-prob(i-1))/(prob(i)-prob(i-1)); x=t(i-1)+w*(t(i)-t(i-1))
        end if
        return
      end if
    end do
    x=t(size(t))
  end function inverse_time_monotone
end module relsurv_data
