! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_utilities
  use mclust_kinds, only : dp
  implicit none
  private
  public :: unmap_classes, majority_vote, match_clusters
contains
  subroutine unmap_classes(classification,z,levels)
    integer,intent(in)::classification(:)
    real(dp),allocatable,intent(out)::z(:,:)
    integer,allocatable,intent(out),optional::levels(:)
    integer,allocatable::lev(:)
    integer::i,j,g
    call unique_sorted(classification,lev)
    g=size(lev)
    allocate(z(size(classification),g)); z=0.0d0
    do i=1,size(classification)
      do j=1,g
        if(classification(i)==lev(j)) then
          z(i,j)=1.0d0
          exit
        end if
      end do
    end do
    if(present(levels)) then
      allocate(levels(g)); levels=lev
    end if
  end subroutine unmap_classes

  subroutine majority_vote(votes,result)
    integer,intent(in)::votes(:,:)
    integer,intent(out)::result(:)
    integer,allocatable::lev(:),row(:)
    integer::i,j,k,best,best_count,c
    if(size(result)/=size(votes,1)) error stop 'majority_vote: result has wrong size'
    allocate(row(size(votes,2)))
    do i=1,size(votes,1)
      row=votes(i,:)
      call unique_sorted(row,lev)
      best=lev(1); best_count=-1
      do j=1,size(lev)
        c=0
        do k=1,size(row)
          if(row(k)==lev(j)) c=c+1
        end do
        if(c>best_count) then
          best_count=c; best=lev(j)
        end if
      end do
      result(i)=best
    end do
  end subroutine majority_vote

  subroutine match_clusters(reference,classification,mapped,mapping,status)
    integer,intent(in)::reference(:),classification(:)
    integer,intent(out)::mapped(:)
    integer,allocatable,intent(out),optional::mapping(:,:)
    integer,intent(out),optional::status
    integer,allocatable::rlev(:),clev(:),counts(:,:),assign(:)
    integer::i,j,k,g,ir,jc
    if(size(reference)/=size(classification) .or. size(mapped)/=size(reference)) then
      if(present(status)) status=-1
      if(size(mapped)>0) mapped=0
      return
    end if
    call unique_sorted(reference,rlev); call unique_sorted(classification,clev)
    g=max(size(rlev),size(clev)); allocate(counts(g,g)); counts=0
    do k=1,size(reference)
      ir=find_value(rlev,reference(k)); jc=find_value(clev,classification(k))
      if(ir>0 .and. jc>0) counts(ir,jc)=counts(ir,jc)+1
    end do
    call maximum_assignment(counts,assign)
    mapped=classification
    do j=1,size(clev)
      i=assign(j)
      if(i>=1 .and. i<=size(rlev)) then
        where(classification==clev(j)) mapped=rlev(i)
      end if
    end do
    if(present(mapping)) then
      allocate(mapping(size(clev),2))
      do j=1,size(clev)
        mapping(j,1)=clev(j)
        i=assign(j)
        if(i>=1 .and. i<=size(rlev)) then
          mapping(j,2)=rlev(i)
        else
          mapping(j,2)=clev(j)
        end if
      end do
    end if
    if(present(status)) status=0
  end subroutine match_clusters

  subroutine maximum_assignment(counts,assignment)
    integer,intent(in)::counts(:,:)
    integer,allocatable,intent(out)::assignment(:)
    real(dp),allocatable::u(:),v(:),minv(:)
    integer,allocatable::p(:),way(:)
    logical,allocatable::used(:)
    real(dp)::cur,delta,maxc
    integer::n,i,j,j0,j1,i0
    n=size(counts,1); maxc=real(maxval(counts),dp)
    allocate(u(0:n),v(0:n),p(0:n),way(0:n)); u=0.0d0; v=0.0d0; p=0; way=0
    do i=1,n
      p(0)=i
      allocate(minv(0:n),used(0:n)); minv=huge(1.0d0); used=.false.; j0=0
      do
        used(j0)=.true.; i0=p(j0); delta=huge(1.0d0); j1=0
        do j=1,n
          if(used(j)) cycle
          cur=(maxc-real(counts(i0,j),dp))-u(i0)-v(j)
          if(cur<minv(j)) then; minv(j)=cur; way(j)=j0; end if
          if(minv(j)<delta) then; delta=minv(j); j1=j; end if
        end do
        do j=0,n
          if(used(j)) then
            u(p(j))=u(p(j))+delta; v(j)=v(j)-delta
          else if(j>0) then
            minv(j)=minv(j)-delta
          end if
        end do
        j0=j1
        if(p(j0)==0) exit
      end do
      do
        j1=way(j0); p(j0)=p(j1); j0=j1
        if(j0==0) exit
      end do
      deallocate(minv,used)
    end do
    allocate(assignment(n)); assignment=0
    do j=1,n
      if(p(j)>0) assignment(j)=p(j)
    end do
  end subroutine maximum_assignment

  subroutine unique_sorted(x,u)
    integer,intent(in)::x(:)
    integer,allocatable,intent(out)::u(:)
    integer,allocatable::tmp(:)
    integer::i,j,key,n
    allocate(tmp(size(x))); tmp=x
    do i=2,size(tmp)
      key=tmp(i); j=i-1
      do while(j>=1)
        if(tmp(j)<=key) exit
        tmp(j+1)=tmp(j); j=j-1
      end do
      tmp(j+1)=key
    end do
    n=0
    do i=1,size(tmp)
      if(i==1) then
        n=n+1
      else if(tmp(i)/=tmp(i-1)) then
        n=n+1
      end if
    end do
    allocate(u(n)); n=0
    do i=1,size(tmp)
      if(i==1) then
        n=n+1; u(n)=tmp(i)
      else if(tmp(i)/=tmp(i-1)) then
        n=n+1; u(n)=tmp(i)
      end if
    end do
  end subroutine unique_sorted

  pure integer function find_value(x,value) result(pos)
    integer,intent(in)::x(:),value
    integer::i
    pos=0
    do i=1,size(x)
      if(x(i)==value) then; pos=i; return; end if
    end do
  end function find_value
end module mclust_utilities
