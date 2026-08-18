module mcmc_temper
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use mcmc_kinds, only : dp
   use mcmc_numerics, only : rand_uniform, rand_normal_vec
   use mcmc_metrop, only : mcmc_scale, output_callback, SCALE_MODE_CONSTANT, &
      SCALE_MODE_DIAGONAL, SCALE_MODE_FULL, scale_constant
   implicit none
   private

   abstract interface
      subroutine temper_density_callback(state,value,data)
         import dp
         real(dp), intent(in) :: state(:)
         real(dp), intent(out) :: value
         class(*), intent(in), optional :: data
      end subroutine temper_density_callback
      subroutine parallel_output_callback(state,value,data)
         import dp
         real(dp), intent(in) :: state(:,:)
         real(dp), intent(out) :: value(:)
         class(*), intent(in), optional :: data
      end subroutine parallel_output_callback
   end interface

   type, public :: temper_serial_result
      real(dp), allocatable :: batch(:,:)
      real(dp), allocatable :: out_batch(:,:)
      real(dp), allocatable :: ibatch(:,:)
      real(dp), allocatable :: acceptx(:)
      real(dp), allocatable :: accepti(:,:)
      real(dp), allocatable :: initial(:)
      real(dp), allocatable :: final(:)
      integer :: status=0
   end type temper_serial_result

   type, public :: temper_parallel_result
      real(dp), allocatable :: batch(:,:,:)
      real(dp), allocatable :: out_batch(:,:)
      real(dp), allocatable :: acceptx(:)
      real(dp), allocatable :: accepti(:,:)
      real(dp), allocatable :: initial(:,:)
      real(dp), allocatable :: final(:,:)
      integer :: status=0
   end type temper_parallel_result

   public :: temper_density_callback, parallel_output_callback
   public :: temper_serial, temper_parallel

contains

   logical function neighbors_valid(neighbors) result(ok)
      logical, intent(in) :: neighbors(:,:)
      integer :: i,k
      k=size(neighbors,1)
      ok=size(neighbors,2)==k
      if (.not. ok) return
      ok=all(neighbors .eqv. transpose(neighbors))
      if (.not. ok) return
      do i=1,k
         if (count(neighbors(i,:))==0) then
            ok=.false.
            return
         end if
      end do
   end function neighbors_valid

   subroutine call_lud(lud,state,value,data)
      procedure(temper_density_callback) :: lud
      real(dp), intent(in) :: state(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: data
      if (present(data)) then
         call lud(state,value,data)
      else
         call lud(state,value)
      end if
   end subroutine call_lud

   subroutine proposal_x(x,proposal,z,scale)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: proposal(size(x)),z(size(x))
      type(mcmc_scale), intent(in) :: scale
      integer :: i,j,d
      d=size(x)
      call rand_normal_vec(z)
      select case(scale%mode)
      case(SCALE_MODE_CONSTANT)
         proposal=x+scale%scalar*z
      case(SCALE_MODE_DIAGONAL)
         if (.not. allocated(scale%diagonal) .or. size(scale%diagonal)/=d) &
            error stop "temper: invalid diagonal scale"
         proposal=x+scale%diagonal*z
      case(SCALE_MODE_FULL)
         if (.not. allocated(scale%full) .or. size(scale%full,1)/=d .or. &
             size(scale%full,2)/=d) error stop "temper: invalid full scale"
         proposal=x
         do i=1,d
            do j=1,d
               proposal(j)=proposal(j)+scale%full(j,i)*z(i)
            end do
         end do
      case default
         error stop "temper: invalid scale mode"
      end select
   end subroutine proposal_x

   function choose_neighbor(neighbors,i) result(j)
      logical, intent(in) :: neighbors(:,:)
      integer, intent(in) :: i
      integer :: j,n,k,target
      n=count(neighbors(i,:))
      if (n<=0) then
         j=0; return
      end if
      target=1+int(rand_uniform()*real(n,dp))
      target=min(n,max(1,target))
      k=0; j=0
      do n=1,size(neighbors,2)
         if (neighbors(i,n)) then
            k=k+1
            if (k==target) then
               j=n; return
            end if
         end if
      end do
   end function choose_neighbor

   function scale_for(scales,i) result(s)
      type(mcmc_scale), intent(in), optional :: scales(:)
      integer, intent(in) :: i
      type(mcmc_scale) :: s
      s=scale_constant(1.0_dp)
      if (present(scales)) then
         if (size(scales)==1) then
            s=scales(1)
         else if (i<=size(scales)) then
            s=scales(i)
         else
            error stop "temper: scale list length must be 1 or ncomp"
         end if
      end if
   end function scale_for

   function temper_serial(lud,initial,neighbors,nbatch,blen,nspac,scales,out_dim,outfun,data) result(res)
      procedure(temper_density_callback) :: lud
      real(dp), intent(in) :: initial(:)
      logical, intent(in) :: neighbors(:,:)
      integer, intent(in) :: nbatch
      integer, intent(in), optional :: blen,nspac,out_dim
      type(mcmc_scale), intent(in), optional :: scales(:)
      procedure(output_callback), optional :: outfun
      class(*), intent(in), optional :: data
      type(temper_serial_result) :: res
      real(dp), allocatable :: state(:),trial(:),xprop(:),z(:),ld(:),bsum(:),obuf(:)
      real(dp), allocatable :: axn(:),axd(:),ain(:,:),aid(:,:),ibsum(:)
      real(dp) :: oldld,newld,logh
      integer :: k,p,bb,ns,ib,jb,is,i,j,comp,od
      logical :: accepted,within
      type(mcmc_scale) :: sc

      k=size(neighbors,1)
      if (.not. neighbors_valid(neighbors)) then
         res%status=1; return
      end if
      p=size(initial)-1
      if (p<=0 .or. nbatch<=0) then
         res%status=2; return
      end if
      comp=nint(initial(1))
      if (comp<1 .or. comp>k) then
         res%status=3; return
      end if
      bb=1; ns=1
      if (present(blen)) bb=blen
      if (present(nspac)) ns=nspac
      state=initial
      allocate(trial(p+1),xprop(p),z(p),ld(k),axn(k),axd(k),ain(k,k),aid(k,k),ibsum(k))
      ld=-huge(1.0_dp); axn=0.0_dp; axd=0.0_dp; ain=0.0_dp; aid=0.0_dp
      call call_lud(lud,state,ld(comp),data)
      if (.not. ieee_is_finite(ld(comp)) .or. ld(comp) <= -0.5_dp*huge(1.0_dp) .or. &
          ld(comp) >= 0.5_dp*huge(1.0_dp)) then
         res%status=4; return
      end if
      allocate(res%ibatch(nbatch,k),res%acceptx(k),res%accepti(k,k))
      if (present(outfun)) then
         if (.not. present(out_dim)) then
            res%status=5; return
         end if
         od=out_dim
         allocate(res%out_batch(nbatch,od),bsum(od),obuf(od))
      else
         allocate(res%batch(nbatch,p),bsum(p),obuf(p))
      end if
      allocate(res%initial(size(initial)))
      res%initial=initial

      do ib=1,nbatch
         bsum=0.0_dp; ibsum=0.0_dp
         do jb=1,bb
            do is=1,ns
               comp=nint(state(1))
               within=rand_uniform()<0.5_dp
               if (within) then
                  sc=scale_for(scales,comp)
                  call proposal_x(state(2:),xprop,z,sc)
                  trial(1)=real(comp,dp); trial(2:)=xprop
                  call call_lud(lud,trial,newld,data)
                  oldld=ld(comp)
                  accepted=.false.
                  if (newld > -huge(1.0_dp)/2.0_dp) then
                     logh=newld-oldld
                     if (logh>=0.0_dp) then
                        accepted=.true.
                     else
                        accepted=rand_uniform()<exp(logh)
                     end if
                  end if
                  if (accepted) then
                     state=trial; ld(comp)=newld; axn(comp)=axn(comp)+1.0_dp
                  end if
                  axd(comp)=axd(comp)+1.0_dp
               else
                  i=comp; j=choose_neighbor(neighbors,i)
                  trial=state; trial(1)=real(j,dp)
                  call call_lud(lud,trial,newld,data)
                  if (ld(i)<=-huge(1.0_dp)/2.0_dp) call call_lud(lud,state,ld(i),data)
                  oldld=ld(i)
                  logh=newld-oldld+log(real(count(neighbors(i,:)),dp)) - &
                       log(real(count(neighbors(j,:)),dp))
                  accepted=.false.
                  if (newld > -huge(1.0_dp)/2.0_dp) then
                     if (logh>=0.0_dp) then
                        accepted=.true.
                     else
                        accepted=rand_uniform()<exp(logh)
                     end if
                  end if
                  if (accepted) then
                     state=trial; ld(j)=newld; ain(i,j)=ain(i,j)+1.0_dp
                  end if
                  aid(i,j)=aid(i,j)+1.0_dp
               end if
            end do
            comp=nint(state(1))
            ibsum(comp)=ibsum(comp)+1.0_dp
            if (present(outfun)) then
               if (present(data)) then
                  call outfun(state,obuf,data)
               else
                  call outfun(state,obuf)
               end if
            else
               obuf=state(2:)
            end if
            bsum=bsum+obuf
         end do
         res%ibatch(ib,:)=ibsum/real(bb,dp)
         if (present(outfun)) then
            res%out_batch(ib,:)=bsum/real(bb,dp)
         else
            res%batch(ib,:)=bsum/real(bb,dp)
         end if
      end do
      do i=1,k
         if (axd(i)>0.0_dp) then
            res%acceptx(i)=axn(i)/axd(i)
         else
            res%acceptx(i)=-1.0_dp
         end if
         do j=1,k
            if (.not. neighbors(i,j) .or. aid(i,j)<=0.0_dp) then
               res%accepti(i,j)=-1.0_dp
            else
               res%accepti(i,j)=ain(i,j)/aid(i,j)
            end if
         end do
      end do
      res%final=state
   end function temper_serial

   function temper_parallel(lud,initial,neighbors,nbatch,blen,nspac,scales,out_dim,outfun,data) result(res)
      procedure(temper_density_callback) :: lud
      real(dp), intent(in) :: initial(:,:)
      logical, intent(in) :: neighbors(:,:)
      integer, intent(in) :: nbatch
      integer, intent(in), optional :: blen,nspac,out_dim
      type(mcmc_scale), intent(in), optional :: scales(:)
      procedure(parallel_output_callback), optional :: outfun
      class(*), intent(in), optional :: data
      type(temper_parallel_result) :: res
      real(dp), allocatable :: state(:,:),trial(:),cop(:),xprop(:),z(:),ld(:), &
         axn(:),axd(:),ain(:,:),aid(:,:),bsum(:,:),obuf(:),osum(:)
      real(dp) :: oldldi,oldldj,newldi,newldj,logh
      integer :: k,p,bb,ns,ib,jb,is,i,j,od
      logical :: accepted,within
      type(mcmc_scale) :: sc

      k=size(initial,1); p=size(initial,2)
      if (.not. neighbors_valid(neighbors) .or. size(neighbors,1)/=k) then
         res%status=1; return
      end if
      bb=1; ns=1
      if (present(blen)) bb=blen
      if (present(nspac)) ns=nspac
      state=initial
      allocate(trial(p+1),cop(p+1),xprop(p),z(p),ld(k),axn(k),axd(k),ain(k,k),aid(k,k))
      axn=0.0_dp; axd=0.0_dp; ain=0.0_dp; aid=0.0_dp
      do i=1,k
         trial(1)=real(i,dp); trial(2:)=state(i,:)
         call call_lud(lud,trial,ld(i),data)
         if (.not. ieee_is_finite(ld(i)) .or. ld(i) <= -0.5_dp*huge(1.0_dp) .or. &
             ld(i) >= 0.5_dp*huge(1.0_dp)) then
            res%status=2; return
         end if
      end do
      allocate(res%acceptx(k),res%accepti(k,k))
      od=1
      if (present(outfun)) then
         if (.not. present(out_dim)) then
            res%status=3; return
         end if
         od=out_dim
         allocate(res%out_batch(nbatch,od))
      else
         allocate(res%batch(nbatch,k,p))
      end if
      allocate(bsum(k,p),obuf(od),osum(od))
      allocate(res%initial(size(initial,1),size(initial,2)))
      res%initial=initial

      do ib=1,nbatch
         if (present(outfun)) then
            osum=0.0_dp
         else
            bsum=0.0_dp
         end if
         do jb=1,bb
            do is=1,ns
               within=rand_uniform()<0.5_dp
               i=1+int(rand_uniform()*real(k,dp)); i=min(k,max(1,i))
               if (within) then
                  sc=scale_for(scales,i)
                  call proposal_x(state(i,:),xprop,z,sc)
                  trial(1)=real(i,dp); trial(2:)=xprop
                  call call_lud(lud,trial,newldi,data)
                  oldldi=ld(i); logh=newldi-oldldi
                  accepted=.false.
                  if (newldi > -huge(1.0_dp)/2.0_dp) then
                     if (logh>=0.0_dp) then
                        accepted=.true.
                     else
                        accepted=rand_uniform()<exp(logh)
                     end if
                  end if
                  if (accepted) then
                     state(i,:)=xprop; ld(i)=newldi; axn(i)=axn(i)+1.0_dp
                  end if
                  axd(i)=axd(i)+1.0_dp
               else
                  j=choose_neighbor(neighbors,i)
                  trial(1)=real(i,dp); trial(2:)=state(j,:)
                  cop(1)=real(j,dp); cop(2:)=state(i,:)
                  call call_lud(lud,trial,newldi,data)
                  call call_lud(lud,cop,newldj,data)
                  oldldi=ld(i); oldldj=ld(j)
                  logh=newldi+newldj-oldldi-oldldj
                  accepted=.false.
                  if (newldi > -huge(1.0_dp)/2.0_dp .and. newldj > -huge(1.0_dp)/2.0_dp) then
                     if (logh>=0.0_dp) then
                        accepted=.true.
                     else
                        accepted=rand_uniform()<exp(logh)
                     end if
                  end if
                  if (accepted) then
                     xprop=state(i,:); state(i,:)=state(j,:); state(j,:)=xprop
                     ld(i)=newldi; ld(j)=newldj; ain(i,j)=ain(i,j)+1.0_dp
                  end if
                  aid(i,j)=aid(i,j)+1.0_dp
               end if
            end do
            if (present(outfun)) then
               if (present(data)) then
                  call outfun(state,obuf,data)
               else
                  call outfun(state,obuf)
               end if
               osum=osum+obuf
            else
               bsum=bsum+state
            end if
         end do
         if (present(outfun)) then
            res%out_batch(ib,:)=osum/real(bb,dp)
         else
            res%batch(ib,:,:)=bsum/real(bb,dp)
         end if
      end do
      do i=1,k
         if (axd(i)>0.0_dp) then
            res%acceptx(i)=axn(i)/axd(i)
         else
            res%acceptx(i)=-1.0_dp
         end if
         do j=1,k
            if (.not. neighbors(i,j) .or. aid(i,j)<=0.0_dp) then
               res%accepti(i,j)=-1.0_dp
            else
               res%accepti(i,j)=ain(i,j)/aid(i,j)
            end if
         end do
      end do
      res%final=state
   end function temper_parallel

end module mcmc_temper
