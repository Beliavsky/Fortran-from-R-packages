! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_probability
  use gpareto_kinds, only : dp
  use gpareto_math, only : normal_cdf
  use gpareto_pareto, only : nondominated_points
  implicit none
  private
  public :: probability_nondomination
contains
  subroutine probability_nondomination(front,mean,sd,pn)
    real(dp),intent(in)::front(:,:),mean(:,:),sd(:,:)
    real(dp),allocatable,intent(out)::pn(:)
    real(dp),allocatable::pf(:,:),phix(:,:),phiy(:,:),phiz(:,:),sub(:,:)
    integer,allocatable::ord(:),idx(:)
    integer::nq,np,m,i,j,k,ns
    real(dp),allocatable::pi2(:),pz(:)
    nq=size(mean,1)
    m=size(mean,2)
    if(size(sd,1)/=nq.or.size(sd,2)/=m.or.size(front,2)/=m) error stop 'probability_nondomination: dimensions'
    if(m/=2.and.m/=3) error stop 'probability_nondomination: only 2 or 3 objectives'
    call nondominated_points(front,pf)
    np=size(pf,1)
    allocate(ord(np))
    ord=[(i,i=1,np)]
    if(m==2) then
      call sort_index(pf(:,1),ord,.false.)
      pf=pf(ord,:)
    else
      call sort_index(pf(:,3),ord,.true.)
      pf=pf(ord,:)
    end if
    deallocate(ord)
    allocate(phix(np,nq),phiy(np,nq))
    do i=1,np
    do j=1,nq
      phix(i,j)=cdf_or_step(pf(i,1),mean(j,1),sd(j,1))
      phiy(i,j)=cdf_or_step(pf(i,2),mean(j,2),sd(j,2))
    end do
    end do
    allocate(pn(nq))
    pn=0.0_dp
    if(m==2) then
      pn=phix(1,:)
      do i=2,np
      pn=pn+(phix(i,:)-phix(i-1,:))*phiy(i-1,:)
      end do
      pn=pn+(1.0_dp-phix(np,:))*phiy(np,:)
      return
    end if
    allocate(phiz(np,nq))
    do i=1,np
    do j=1,nq
      phiz(i,j)=cdf_or_step(pf(i,3),mean(j,3),sd(j,3))
    end do
    end do
    do i=1,np
      if(i==1) then
        pz=1.0_dp-phiz(i,:)
      else
        pz=phiz(i-1,:)-phiz(i,:)
      end if
      call nondominated_indices_local(pf(i:np,1:2),idx)
      ns=size(idx)
      allocate(sub(ns,2))
      do k=1,ns
      sub(k,:)=pf(i+idx(k)-1,1:2)
      end do
      allocate(ord(ns))
      ord=[(k,k=1,ns)]
      call sort_index(sub(:,1),ord,.false.)
      sub=sub(ord,:)
      allocate(pi2(nq))
      do j=1,nq
        pi2(j)=cdf_or_step(sub(1,1),mean(j,1),sd(j,1))
        do k=2,ns
          pi2(j)=pi2(j)+(cdf_or_step(sub(k,1),mean(j,1),sd(j,1))- &
            cdf_or_step(sub(k-1,1),mean(j,1),sd(j,1)))*cdf_or_step(sub(k-1,2),mean(j,2),sd(j,2))
        end do
        pi2(j)=pi2(j)+(1.0_dp-cdf_or_step(sub(ns,1),mean(j,1),sd(j,1)))* &
          cdf_or_step(sub(ns,2),mean(j,2),sd(j,2))
      end do
      pn=pn+pi2*pz
      deallocate(sub,ord,pi2,idx,pz)
    end do
    pn=pn+phiz(np,:)
  end subroutine probability_nondomination

  pure real(dp) function cdf_or_step(a,m,s) result(p)
    real(dp),intent(in)::a,m,s
    if(s<=sqrt(tiny(1.0_dp))) then
      if(a>m) then
      p=1.0_dp
      else
      p=0.0_dp
      end if
    else
      p=normal_cdf((a-m)/s)
    end if
  end function cdf_or_step

  subroutine sort_index(x,idx,decreasing)
    real(dp),intent(in)::x(:)
    integer,intent(inout)::idx(:)
    logical,intent(in)::decreasing
    integer::i,j,t
    do i=2,size(idx)
      t=idx(i)
      j=i-1
      do while(j>=1)
        if((.not.decreasing.and.x(idx(j))<=x(t)).or.(decreasing.and.x(idx(j))>=x(t))) exit
        idx(j+1)=idx(j)
        j=j-1
      end do
      idx(j+1)=t
    end do
  end subroutine sort_index

  subroutine nondominated_indices_local(p,idx)
    real(dp),intent(in)::p(:,:)
    integer,allocatable,intent(out)::idx(:)
    logical,allocatable::keep(:)
    integer::i,j,k
    allocate(keep(size(p,1)))
    keep=.true.
    do i=1,size(p,1)
    do j=1,size(p,1)
      if(i==j)cycle
      if(all(p(j,:)<=p(i,:)).and.any(p(j,:)<p(i,:)))then
      keep(i)=.false.
      exit
      end if
    end do
    end do
    allocate(idx(count(keep)))
    k=0
    do i=1,size(keep)
    if(keep(i))then
    k=k+1
    idx(k)=i
    end if
    end do
  end subroutine nondominated_indices_local
end module gpareto_probability
