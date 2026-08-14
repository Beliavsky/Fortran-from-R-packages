! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_impute
  use mclust_kinds, only : dp, pi_dp
  use mclust_types, only : mclust_fit
  use mclust_math, only : logsumexp
  implicit none
  private
  public :: impute_data

contains

  subroutine impute_data(fit,x,missing,ximp,status)
    type(mclust_fit),intent(in)::fit
    real(dp),intent(in)::x(:,:)
    logical,intent(in)::missing(:,:)
    real(dp),allocatable,intent(out)::ximp(:,:)
    integer,intent(out),optional::status
    integer,allocatable::obs(:),mis(:)
    real(dp),allocatable::lp(:),cm(:),s_oo(:,:),b(:),delta(:),mucond(:)
    integer::i,j,k,no,nm,info,d,g,jo,jm
    real(dp)::ld,q,ls

    d=size(x,2); g=fit%g
    if(size(x,1)/=size(missing,1) .or. d/=size(missing,2) .or. d/=fit%d) then
      allocate(ximp(0,0)); if(present(status)) status=-1; return
    end if
    allocate(ximp(size(x,1),d)); ximp=x
    do i=1,size(x,1)
      nm=count(missing(i,:)); if(nm==0) cycle
      no=d-nm; allocate(obs(no),mis(nm),lp(g),cm(nm),mucond(nm))
      jo=0; jm=0
      do j=1,d
        if(missing(i,j)) then; jm=jm+1; mis(jm)=j
        else; jo=jo+1; obs(jo)=j; end if
      end do
      lp=0.0_dp; cm=0.0_dp
      do k=1,g
        if(no==0) then
          lp(k)=log(max(fit%pro(k),tiny(1.0_dp)))
          mucond=fit%mean(mis,k)
        else
          allocate(s_oo(no,no),b(no),delta(no))
          do j=1,no
            s_oo(j,:)=fit%sigma(obs(j),obs,k)
          end do
          call chol_lower(s_oo,info)
          if(info/=0) then; if(present(status)) status=10+k; return; end if
          delta=x(i,obs)-fit%mean(obs,k); b=delta
          call forward_solve(s_oo,b)
          q=dot_product(b,b)
          ld=2.0_dp*sum(log([(s_oo(j,j),j=1,no)]))
          lp(k)=log(max(fit%pro(k),tiny(1.0_dp)))-0.5_dp*(no*log(2.0_dp*pi_dp)+ld+q)
          call solve_spd_chol(s_oo,delta)
          do jm=1,nm
            mucond(jm)=fit%mean(mis(jm),k)+dot_product(fit%sigma(mis(jm),obs,k),delta)
          end do
          deallocate(s_oo,b,delta)
        end if
        ls=lp(k) ! temporary store does not preserve conditional mean per k
        ! accumulate later after posterior is known: stash in ximp impossible; recompute below
      end do
      ls=logsumexp(lp); lp=exp(lp-ls)
      cm=0.0_dp
      do k=1,g
        if(no==0) then
          mucond=fit%mean(mis,k)
        else
          allocate(s_oo(no,no),delta(no))
          do j=1,no; s_oo(j,:)=fit%sigma(obs(j),obs,k); end do
          call chol_lower(s_oo,info)
          delta=x(i,obs)-fit%mean(obs,k); call solve_spd_chol(s_oo,delta)
          do jm=1,nm
            mucond(jm)=fit%mean(mis(jm),k)+dot_product(fit%sigma(mis(jm),obs,k),delta)
          end do
          deallocate(s_oo,delta)
        end if
        cm=cm+lp(k)*mucond
      end do
      ximp(i,mis)=cm
      deallocate(obs,mis,lp,cm,mucond)
    end do
    if(present(status)) status=0
  end subroutine impute_data

  subroutine chol_lower(a,info)
    real(dp),intent(inout)::a(:,:)
    integer,intent(out)::info
    integer::i,j,k
    real(dp)::s
    info=0
    do j=1,size(a,1)
      s=a(j,j); do k=1,j-1; s=s-a(j,k)**2; end do
      if(s<=0.0_dp) then; info=j; return; end if
      a(j,j)=sqrt(s)
      do i=j+1,size(a,1)
        s=a(i,j); do k=1,j-1; s=s-a(i,k)*a(j,k); end do
        a(i,j)=s/a(j,j)
      end do
      if(j<size(a,1)) a(j,j+1:)=0.0_dp
    end do
  end subroutine chol_lower
  subroutine forward_solve(l,b)
    real(dp),intent(in)::l(:,:); real(dp),intent(inout)::b(:); integer::i
    do i=1,size(b); if(i>1)b(i)=b(i)-dot_product(l(i,1:i-1),b(1:i-1)); b(i)=b(i)/l(i,i); end do
  end subroutine forward_solve
  subroutine back_solve_lt(l,b)
    real(dp),intent(in)::l(:,:); real(dp),intent(inout)::b(:); integer::i
    do i=size(b),1,-1; if(i<size(b))b(i)=b(i)-dot_product(l(i+1:,i),b(i+1:)); b(i)=b(i)/l(i,i); end do
  end subroutine back_solve_lt
  subroutine solve_spd_chol(l,b)
    real(dp),intent(in)::l(:,:); real(dp),intent(inout)::b(:)
    call forward_solve(l,b); call back_solve_lt(l,b)
  end subroutine solve_spd_chol
end module mclust_impute
