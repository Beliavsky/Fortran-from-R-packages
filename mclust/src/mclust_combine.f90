! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_combine
  use mclust_kinds, only : dp
  implicit none
  private
  type, public :: cluster_combination
    integer :: g=0
    integer,allocatable :: merge(:,:)
    real(dp),allocatable :: entropy(:)
  end type
  public :: clust_combi, entropy_z, apply_combination
contains
  pure real(dp) function entropy_z(z) result(e)
    real(dp),intent(in)::z(:,:)
    integer::i,j
    e=0.0_dp
    do i=1,size(z,1); do j=1,size(z,2); if(z(i,j)>0.0_dp)e=e-z(i,j)*log(z(i,j)); end do; end do
  end function entropy_z

  subroutine clust_combi(z,out)
    real(dp),intent(in)::z(:,:)
    type(cluster_combination),intent(out)::out
    real(dp),allocatable::cur(:,:),cand(:,:)
    integer,allocatable::labels(:)
    integer::g,s,a,b,ba,bb,j
    real(dp)::best,e
    g=size(z,2); out%g=g
    allocate(out%merge(max(0,g-1),2),out%entropy(max(1,g)),labels(g)); labels=[(j,j=1,g)]
    cur=z; out%entropy(1)=entropy_z(cur)
    do s=1,g-1
      best=huge(1.0_dp); ba=1; bb=2
      do a=1,size(cur,2)-1
        do b=a+1,size(cur,2)
          call merge_columns(cur,a,b,cand)
          e=entropy_z(cand)
          if(e<best) then; best=e; ba=a; bb=b; end if
          deallocate(cand)
        end do
      end do
      out%merge(s,:)=[labels(ba),labels(bb)]
      call merge_columns(cur,ba,bb,cand); call move_alloc(cand,cur)
      labels(ba)=-s
      labels(bb:size(labels)-1)=labels(bb+1:size(labels))
      labels=labels(:size(labels)-1)
      out%entropy(s+1)=best
    end do
  end subroutine clust_combi

  subroutine apply_combination(z,merge,n_groups,zout)
    real(dp),intent(in)::z(:,:)
    integer,intent(in)::merge(:,:),n_groups
    real(dp),allocatable,intent(out)::zout(:,:)
    real(dp),allocatable::cur(:,:),tmp(:,:)
    integer,allocatable::labels(:)
    integer::s,a,b,i
    cur=z; allocate(labels(size(z,2))); labels=[(i,i=1,size(z,2))]
    do s=1,size(z,2)-n_groups
      a=find_label(labels,merge(s,1)); b=find_label(labels,merge(s,2))
      if(a==0 .or. b==0) exit
      if(a>b) then; i=a; a=b; b=i; end if
      call merge_columns(cur,a,b,tmp); call move_alloc(tmp,cur)
      labels(a)=-s; labels(b:size(labels)-1)=labels(b+1:size(labels)); labels=labels(:size(labels)-1)
    end do
    call move_alloc(cur,zout)
  end subroutine apply_combination

  subroutine merge_columns(z,a,b,out)
    real(dp),intent(in)::z(:,:); integer,intent(in)::a,b
    real(dp),allocatable,intent(out)::out(:,:)
    integer::j,k
    allocate(out(size(z,1),size(z,2)-1)); k=0
    do j=1,size(z,2)
      if(j==b) cycle
      k=k+1
      if(j==a) then; out(:,k)=z(:,a)+z(:,b); else; out(:,k)=z(:,j); end if
    end do
  end subroutine merge_columns
  pure integer function find_label(labels,v) result(k)
    integer,intent(in)::labels(:),v; integer::i
    k=0; do i=1,size(labels); if(labels(i)==v) then; k=i; return; end if; end do
  end function find_label
end module mclust_combine
