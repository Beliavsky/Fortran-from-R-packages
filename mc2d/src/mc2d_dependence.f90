! SPDX-License-Identifier: GPL-2.0-or-later
!
! Rank-correlation induction using the Iman-Conover method used by mc2d::cornode.
module mc2d_dependence
  use mc2d_kinds, only : dp
  use mc2d_utils, only : argsort_dp, ranks_dp
  use mc2d_random, only : seed_random
  use mc2d_node, only : mcnode, mc_type_v, mc_type_u, mc_type_vu
  use mvtnorm_linalg, only : cholesky_lower
  implicit none
  private
  public :: cornode, cornode_matrix, cornode_nodes

  interface cornode
    module procedure cornode_matrix
    module procedure cornode_nodes
  end interface
contains
  subroutine random_permutation(n, p)
    integer, intent(in) :: n
    integer, intent(out) :: p(n)
    integer :: i, j, t
    real(dp) :: u
    p = [(i, i=1,n)]
    do i=n,2,-1
      call random_number(u)
      j = 1 + int(u*real(i,dp))
      if (j > i) j=i
      t=p(i); p(i)=p(j); p(j)=t
    end do
  end subroutine random_permutation

  function cornode_matrix(data, target, seed, rank_index) result(res)
    real(dp), intent(in) :: data(:,:), target(:,:)
    integer, intent(in), optional :: seed
    integer, allocatable, intent(out), optional :: rank_index(:,:)
    real(dp), allocatable :: res(:,:)
    real(dp), allocatable :: l(:,:), u(:,:), rr(:,:), ret(:,:), rk(:)
    integer, allocatable :: perm(:), xorder(:,:), rang(:,:), first_order(:)
    logical :: ok
    character(len=128) :: message
    integer :: n, p, i, j

    n=size(data,1); p=size(data,2)
    if (n < 2) error stop 'cornode: at least two rows are required'
    if (p < 2) error stop 'cornode: at least two columns are required'
    if (size(target,1)/=p .or. size(target,2)/=p) &
      error stop 'cornode: target has incorrect dimensions'
    if (present(seed)) call seed_random(seed)

    call cholesky_lower(target,l,ok,message)
    if (.not.ok) error stop 'cornode: target must be positive definite'
    u=transpose(l)
    allocate(rr(n,p),ret(n,p),rk(n),perm(n),xorder(n,p),rang(n,p))

    ! R builds random permutations of normal scores and then ranks them.
    ! Because the scores are strictly increasing, their ranks are exactly the
    ! integer permutation, so the normal-score values need not be materialized.
    do j=1,p
      call random_permutation(n,perm)
      rr(:,j)=real(perm,dp)
      call argsort_dp(data(:,j),xorder(:,j))
    end do
    ret=matmul(rr,u)
    do j=1,p
      call ranks_dp(ret(:,j),rk)
      do i=1,n
        rang(i,j)=xorder(nint(rk(i)),j)
      end do
    end do

    ! Match mc2d: restore the first variable's original row order.
    allocate(first_order(n))
    call argsort_int(rang(:,1),first_order)
    rang=rang(first_order,:)

    allocate(res(n,p))
    do j=1,p
      do i=1,n
        res(i,j)=data(rang(i,j),j)
      end do
    end do
    if (present(rank_index)) then
      allocate(rank_index(n,p)); rank_index=rang
    end if
  end function cornode_matrix

  function cornode_nodes(nodes, target, seed, outrank) result(res)
    type(mcnode), intent(in) :: nodes(:)
    real(dp), intent(in) :: target(:,:)
    integer, intent(in), optional :: seed
    logical, intent(in), optional :: outrank
    type(mcnode), allocatable :: res(:)
    real(dp), allocatable :: mat(:,:), sorted(:,:)
    integer :: p,nvar,nunc,nva,i,j,k,ir,nv_count
    logical :: all_u, rank_out

    p=size(nodes)
    if (p < 2) error stop 'cornode: at least two mcnodes are required'
    if (size(target,1)/=p .or. size(target,2)/=p) &
      error stop 'cornode: target has incorrect dimensions'
    all_u=all([(nodes(i)%node_type==mc_type_u,i=1,p)])
    if (.not.all_u) then
      if (.not.all([(nodes(i)%node_type==mc_type_v .or. &
                     nodes(i)%node_type==mc_type_vu,i=1,p)])) &
        error stop 'cornode: use all U nodes or a valid V/VU combination'
      nv_count=count([(nodes(i)%node_type==mc_type_v,i=1,p)])
      if (any([(nodes(i)%node_type==mc_type_vu,i=1,p)]) .and. nv_count>0) then
        if (nv_count/=1 .or. nodes(1)%node_type/=mc_type_v) &
          error stop 'cornode: with VU nodes, only the first node may be V'
      end if
    end if
    nva=nodes(1)%nvariates()
    if (any([(nodes(i)%nvariates()/=nva,i=1,p)])) &
      error stop 'cornode: inconsistent numbers of variates'
    if (present(seed)) call seed_random(seed)
    allocate(res(p)); res=nodes
    rank_out=.false.; if(present(outrank)) rank_out=outrank

    if (all_u) then
      nvar=maxval([(nodes(i)%nsu(),i=1,p)])
      if (any([(nodes(i)%nsu()/=nvar,i=1,p)])) &
        error stop 'cornode: inconsistent uncertainty dimensions'
      allocate(mat(nvar,p))
      do k=1,nva
        do i=1,p; mat(:,i)=nodes(i)%value(1,:,k); end do
        sorted=cornode_matrix(mat,target)
        do i=2,p; res(i)%value(1,:,k)=sorted(:,i); end do
      end do
      if (rank_out) then
        do k=1,nva; res(1)%value(1,:,k)=real([(j,j=1,nvar)],dp); end do
      end if
    else
      nvar=maxval([(nodes(i)%nsv(),i=1,p)])
      nunc=maxval([(nodes(i)%nsu(),i=1,p)])
      if (any([(nodes(i)%nsv()/=nvar,i=1,p)])) &
        error stop 'cornode: inconsistent variability dimensions'
      if (any([(nodes(i)%nsu()/=1 .and. nodes(i)%nsu()/=nunc,i=1,p)])) &
        error stop 'cornode: inconsistent uncertainty dimensions'
      allocate(mat(nvar,p))
      do k=1,nva
        do j=1,nunc
          do i=1,p
            ir=min(j,nodes(i)%nsu())
            mat(:,i)=nodes(i)%value(:,ir,k)
          end do
          sorted=cornode_matrix(mat,target)
          do i=2,p
            if (res(i)%nsu()==1) then
              if (j==1) res(i)%value(:,1,k)=sorted(:,i)
            else
              res(i)%value(:,j,k)=sorted(:,i)
            end if
          end do
        end do
      end do
      if (rank_out) then
        do k=1,nva
          res(1)%value(:,1,k)=real([(j,j=1,nvar)],dp)
        end do
      end if
    end if
  end function cornode_nodes

  subroutine argsort_int(x,idx)
    integer, intent(in) :: x(:)
    integer, intent(out) :: idx(size(x))
    integer :: i,j,t
    idx=[(i,i=1,size(x))]
    do i=2,size(x)
      t=idx(i); j=i-1
      do while(j>=1)
        if(x(idx(j))<=x(t)) exit
        idx(j+1)=idx(j); j=j-1
      end do
      idx(j+1)=t
    end do
  end subroutine argsort_int
end module mc2d_dependence
