! Problem construction and validation for the Rdsdp/DSDP5 translation.
! DSDP copyright/license: see licenses/DSDP-LICENSE.
module rdsdp_problem_mod
   use rdsdp_kinds, only : dp
   use rdsdp_types, only : dsdp_problem, dsdp_block, dsdp_sdp_block, dsdp_lp_block, &
      dsdp_data_dense, dsdp_data_sparse, dsdp_data_lowrank
   use rdsdp_data, only : compress_sdp_block, densify_sdp_block
   implicit none
   private
   public :: dsdp_from_sedumi, validate_problem, symmetrize_problem, compress_problem, densify_problem

contains

   subroutine validate_problem(prob,ok,message)
      type(dsdp_problem), intent(in) :: prob
      logical, intent(out) :: ok
      character(len=:), allocatable, intent(out) :: message
      integer :: k,i,storage
      ok = .false.
      if (prob%m <= 0) then
         message='number of constraints must be positive'; return
      end if
      if (.not.allocated(prob%b) .or. size(prob%b)/=prob%m) then
         message='b has wrong size'; return
      end if
      if (.not.allocated(prob%block) .or. size(prob%block)==0) then
         message='problem has no cone blocks'; return
      end if
      do k=1,size(prob%block)
         if (prob%block(k)%n<=0) then
            message='block size must be positive'; return
         end if
         select case(prob%block(k)%category)
         case(dsdp_sdp_block)
            if (.not.allocated(prob%block(k)%a_storage) .or. size(prob%block(k)%a_storage)/=prob%m) then
               ! Legacy dense problem objects remain accepted.
               if (.not.allocated(prob%block(k)%a)) then
                  message='SDP constraint storage not allocated'; return
               end if
            end if
            storage=prob%block(k)%c_storage
            if (.not.valid_data(prob%block(k),0,storage)) then
               message='SDP objective block data invalid'; return
            end if
            do i=1,prob%m
               if (allocated(prob%block(k)%a_storage)) then
                  storage=prob%block(k)%a_storage(i)
               else
                  storage=dsdp_data_dense
               end if
               if (.not.valid_data(prob%block(k),i,storage)) then
                  message='SDP constraint block data invalid'; return
               end if
            end do
         case(dsdp_lp_block)
            if (.not.allocated(prob%block(k)%cdiag) .or. .not.allocated(prob%block(k)%adiag)) then
               message='LP block data not allocated'; return
            end if
            if (size(prob%block(k)%cdiag)/=prob%block(k)%n) then
               message='LP objective block has wrong dimensions'; return
            end if
            if (size(prob%block(k)%adiag,1)/=prob%block(k)%n .or. size(prob%block(k)%adiag,2)/=prob%m) then
               message='LP constraint block has wrong dimensions'; return
            end if
         case default
            message='unknown block category'; return
         end select
      end do
      ok=.true.; message='ok'
   end subroutine validate_problem

   logical function valid_data(block,idx,storage) result(ok)
      type(dsdp_block), intent(in) :: block
      integer, intent(in) :: idx,storage
      ok=.false.
      select case(storage)
      case(dsdp_data_dense)
         if (idx==0) then
            ok=allocated(block%c)
            if (ok) ok=size(block%c,1)==block%n .and. size(block%c,2)==block%n
         else
            ok=allocated(block%a)
            if (ok) ok=size(block%a,1)==block%n .and. size(block%a,2)==block%n .and. size(block%a,3)>=idx
         end if
      case(dsdp_data_sparse)
         if (idx==0) then
            ok=block%c_sparse%n==block%n .and. block%c_sparse%nnz>=0
            if (ok .and. block%c_sparse%nnz>0) then
               ok=allocated(block%c_sparse%row) .and. allocated(block%c_sparse%col) .and. &
                  allocated(block%c_sparse%val)
            end if
         else
            ok=allocated(block%a_sparse)
            if (ok) ok=size(block%a_sparse)>=idx
            if (ok) ok=block%a_sparse(idx)%n==block%n .and. block%a_sparse(idx)%nnz>=0
         end if
      case(dsdp_data_lowrank)
         if (idx==0) then
            ok=block%c_lowrank%n==block%n .and. block%c_lowrank%rank>=0
            if (ok .and. block%c_lowrank%rank>0) ok=allocated(block%c_lowrank%coeff) .and. allocated(block%c_lowrank%vec)
         else
            ok=allocated(block%a_lowrank)
            if (ok) ok=size(block%a_lowrank)>=idx
            if (ok) ok=block%a_lowrank(idx)%n==block%n .and. block%a_lowrank(idx)%rank>=0
         end if
      end select
   end function valid_data

   subroutine symmetrize_problem(prob)
      type(dsdp_problem), intent(inout) :: prob
      integer :: k,j
      do k=1,size(prob%block)
         if (prob%block(k)%category==dsdp_sdp_block) then
            if (allocated(prob%block(k)%c) .and. prob%block(k)%c_storage==dsdp_data_dense) &
               prob%block(k)%c=0.5_dp*(prob%block(k)%c+transpose(prob%block(k)%c))
            if (allocated(prob%block(k)%a)) then
               do j=1,prob%m
                  if (.not.allocated(prob%block(k)%a_storage) .or. prob%block(k)%a_storage(j)==dsdp_data_dense) &
                     prob%block(k)%a(:,:,j)=0.5_dp*(prob%block(k)%a(:,:,j)+ &
                        transpose(prob%block(k)%a(:,:,j)))
               end do
            end if
         end if
      end do
   end subroutine symmetrize_problem

   subroutine compress_problem(prob,threshold,release_dense)
      type(dsdp_problem), intent(inout) :: prob
      real(dp), intent(in), optional :: threshold
      logical, intent(in), optional :: release_dense
      real(dp) :: t
      logical :: rel
      integer :: k
      t=0.20_dp; if (present(threshold)) t=max(0.0_dp,min(1.0_dp,threshold))
      rel=.true.; if (present(release_dense)) rel=release_dense
      do k=1,size(prob%block)
         if (prob%block(k)%category==dsdp_sdp_block) call compress_sdp_block(prob%block(k),prob%m,t,rel)
      end do
   end subroutine compress_problem

   subroutine densify_problem(prob)
      type(dsdp_problem), intent(inout) :: prob
      integer :: k
      do k=1,size(prob%block)
         if (prob%block(k)%category==dsdp_sdp_block) call densify_sdp_block(prob%block(k),prob%m)
      end do
   end subroutine densify_problem

   subroutine dsdp_from_sedumi(amat,b,c,l,s,prob)
      ! Construct the same mathematical problem accepted by Rdsdp::dsdp().
      real(dp), intent(in) :: amat(:,:), b(:), c(:)
      integer, intent(in) :: l
      integer, intent(in) :: s(:)
      type(dsdp_problem), intent(out) :: prob
      integer :: m,ncols,nb,k,i,j,q,off,p
      m=size(amat,1); ncols=size(amat,2)
      if (size(b)/=m) error stop 'dsdp_from_sedumi: b has wrong size'
      if (size(c)/=ncols) error stop 'dsdp_from_sedumi: C has wrong size'
      q=l
      do k=1,size(s); q=q+s(k)*s(k); end do
      if (q/=ncols) error stop 'dsdp_from_sedumi: cone dimensions do not match A/C'
      nb=size(s)+merge(1,0,l>0)
      prob%m=m
      allocate(prob%b(m),prob%block(nb)); prob%b=b
      off=0; q=0
      if (l>0) then
         q=q+1
         prob%block(q)%category=dsdp_lp_block
         prob%block(q)%n=l
         allocate(prob%block(q)%cdiag(l),prob%block(q)%adiag(l,m))
         prob%block(q)%cdiag=c(1:l)
         do i=1,m
            prob%block(q)%adiag(:,i)=amat(i,1:l)
         end do
         off=l
      end if
      do k=1,size(s)
         q=q+1
         prob%block(q)%category=dsdp_sdp_block
         prob%block(q)%n=s(k)
         prob%block(q)%c_storage=dsdp_data_dense
         allocate(prob%block(q)%a_storage(m)); prob%block(q)%a_storage=dsdp_data_dense
         allocate(prob%block(q)%c(s(k),s(k)),prob%block(q)%a(s(k),s(k),m))
         do i=1,s(k)
            do j=1,s(k)
               p=off+(i-1)*s(k)+j
               prob%block(q)%c(i,j)=c(p)
               prob%block(q)%a(i,j,:)=amat(:,p)
            end do
         end do
         off=off+s(k)*s(k)
      end do
      call symmetrize_problem(prob)
      ! Automatically compress sparse SeDuMi blocks while preserving dense matrices
      ! whose density makes BLAS more appropriate.
      call compress_problem(prob,0.20_dp,.true.)
   end subroutine dsdp_from_sedumi

end module rdsdp_problem_mod
