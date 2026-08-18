module mcmc_metrop
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use mcmc_kinds, only : dp
   use mcmc_numerics, only : rand_uniform, rand_normal_vec
   implicit none
   private

   integer, parameter, public :: SCALE_MODE_CONSTANT = 1
   integer, parameter, public :: SCALE_MODE_DIAGONAL = 2
   integer, parameter, public :: SCALE_MODE_FULL = 3

   abstract interface
      subroutine log_density_callback(state, value, data)
         import dp
         real(dp), intent(in) :: state(:)
         real(dp), intent(out) :: value
         class(*), intent(in), optional :: data
      end subroutine log_density_callback
      subroutine output_callback(state, value, data)
         import dp
         real(dp), intent(in) :: state(:)
         real(dp), intent(out) :: value(:)
         class(*), intent(in), optional :: data
      end subroutine output_callback
   end interface

   type, public :: mcmc_scale
      integer :: mode = SCALE_MODE_CONSTANT
      real(dp) :: scalar = 1.0_dp
      real(dp), allocatable :: diagonal(:)
      real(dp), allocatable :: full(:,:)
   end type mcmc_scale

   type, public :: metrop_result
      real(dp) :: accept = 0.0_dp
      real(dp), allocatable :: accept_batch(:)
      real(dp), allocatable :: batch(:,:)
      real(dp), allocatable :: initial(:)
      real(dp), allocatable :: final(:)
      integer :: nbatch = 0
      integer :: blen = 0
      integer :: nspac = 0
      type(mcmc_scale) :: scale
      logical :: debug = .false.
      real(dp), allocatable :: current(:,:)
      real(dp), allocatable :: proposal(:,:)
      real(dp), allocatable :: log_green(:)
      real(dp), allocatable :: u(:)
      real(dp), allocatable :: z(:,:)
      logical, allocatable :: debug_accept(:)
      integer :: status = 0
   end type metrop_result

   public :: log_density_callback, output_callback
   public :: scale_constant, scale_diagonal, scale_full
   public :: metrop

contains

   function scale_constant(value) result(s)
      real(dp), intent(in) :: value
      type(mcmc_scale) :: s
      s%mode = SCALE_MODE_CONSTANT
      s%scalar = value
   end function scale_constant

   function scale_diagonal(value) result(s)
      real(dp), intent(in) :: value(:)
      type(mcmc_scale) :: s
      s%mode = SCALE_MODE_DIAGONAL
      allocate(s%diagonal(size(value)))
      s%diagonal = value
   end function scale_diagonal

   function scale_full(value) result(s)
      real(dp), intent(in) :: value(:,:)
      type(mcmc_scale) :: s
      s%mode = SCALE_MODE_FULL
      allocate(s%full(size(value,1),size(value,2)))
      s%full = value
   end function scale_full

   logical function scale_valid(s,d) result(ok)
      type(mcmc_scale), intent(in) :: s
      integer, intent(in) :: d
      select case (s%mode)
      case (SCALE_MODE_CONSTANT)
         ok = ieee_is_finite(s%scalar)
      case (SCALE_MODE_DIAGONAL)
         ok = allocated(s%diagonal) .and. size(s%diagonal) == d
         if (ok) ok = all(ieee_is_finite(s%diagonal))
      case (SCALE_MODE_FULL)
         ok = allocated(s%full) .and. size(s%full,1) == d .and. size(s%full,2) == d
         if (ok) ok = all(ieee_is_finite(s%full))
      case default
         ok = .false.
      end select
   end function scale_valid

   subroutine propose(state,proposal,z,scale)
      real(dp), intent(in) :: state(:)
      real(dp), intent(out) :: proposal(size(state)),z(size(state))
      type(mcmc_scale), intent(in) :: scale
      integer :: i,j
      call rand_normal_vec(z)
      select case(scale%mode)
      case(SCALE_MODE_CONSTANT)
         proposal = state + scale%scalar*z
      case(SCALE_MODE_DIAGONAL)
         proposal = state + scale%diagonal*z
      case(SCALE_MODE_FULL)
         proposal = state
         ! C code reads the R column-major matrix sequentially for each z(i),
         ! which is mathematically scale %*% z.
         do i = 1, size(state)
            do j = 1, size(state)
               proposal(j) = proposal(j) + scale%full(j,i)*z(i)
            end do
         end do
      end select
   end subroutine propose

   subroutine call_lud(lud,state,value,data)
      procedure(log_density_callback) :: lud
      real(dp), intent(in) :: state(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: data
      if (present(data)) then
         call lud(state,value,data)
      else
         call lud(state,value)
      end if
   end subroutine call_lud

   subroutine call_out(outfun,state,value,data)
      procedure(output_callback), optional :: outfun
      real(dp), intent(in) :: state(:)
      real(dp), intent(out) :: value(:)
      class(*), intent(in), optional :: data
      if (present(outfun)) then
         if (present(data)) then
            call outfun(state,value,data)
         else
            call outfun(state,value)
         end if
      else
         value = state
      end if
   end subroutine call_out

   function metrop(lud,initial,nbatch,blen,nspac,scale,out_dim,outfun,data,debug) result(res)
      procedure(log_density_callback) :: lud
      real(dp), intent(in) :: initial(:)
      integer, intent(in) :: nbatch
      integer, intent(in), optional :: blen,nspac,out_dim
      type(mcmc_scale), intent(in), optional :: scale
      procedure(output_callback), optional :: outfun
      class(*), intent(in), optional :: data
      logical, intent(in), optional :: debug
      type(metrop_result) :: res
      type(mcmc_scale) :: sc
      real(dp), allocatable :: state(:),proposal(:),z(:),obuf(:),bsum(:)
      real(dp) :: current_ld,proposal_ld,green,u
      real(dp) :: accept_count,tries,abatch,tbatch
      integer :: ib,jb,is,d,od,bb,ns,iter,niter
      logical :: accepted,dbg

      res%status = 0
      d = size(initial)
      bb = 1
      ns = 1
      if (present(blen)) bb = blen
      if (present(nspac)) ns = nspac
      if (nbatch <= 0 .or. bb <= 0 .or. ns <= 0 .or. d <= 0 .or. &
          any(.not. ieee_is_finite(initial))) then
         res%status = 1
         return
      end if

      sc = scale_constant(1.0_dp)
      if (present(scale)) sc = scale
      if (.not. scale_valid(sc,d)) then
         res%status = 2
         return
      end if

      if (present(outfun)) then
         if (.not. present(out_dim)) then
            res%status = 3
            return
         end if
         od = out_dim
      else
         od = d
      end if
      if (od <= 0) then
         res%status = 4
         return
      end if
      dbg = .false.
      if (present(debug)) dbg = debug

      allocate(state(d),proposal(d),z(d),obuf(od),bsum(od))
      state = initial
      call call_lud(lud,state,current_ld,data)
      if (.not. ieee_is_finite(current_ld) .or. current_ld <= -0.5_dp*huge(1.0_dp) .or. &
          current_ld >= 0.5_dp*huge(1.0_dp)) then
         res%status = 5
         return
      end if

      res%initial = initial
      allocate(res%batch(nbatch,od),res%accept_batch(nbatch))
      res%nbatch = nbatch
      res%blen = bb
      res%nspac = ns
      res%scale = sc
      res%debug = dbg
      niter = nbatch*bb*ns
      if (dbg) then
         allocate(res%current(niter,d),res%proposal(niter,d),res%log_green(niter), &
                  res%u(niter),res%z(niter,d),res%debug_accept(niter))
      end if

      accept_count = 0.0_dp
      tries = 0.0_dp
      iter = 0
      do ib = 1, nbatch
         bsum = 0.0_dp
         abatch = 0.0_dp
         tbatch = 0.0_dp
         do jb = 1, bb
            do is = 1, ns
               iter = iter+1
               call propose(state,proposal,z,sc)
               call call_lud(lud,proposal,proposal_ld,data)
               accepted = .false.
               u = -1.0_dp
               green = proposal_ld-current_ld
               if (proposal_ld > -huge(1.0_dp)/2.0_dp) then
                  if (.not. ieee_is_finite(proposal_ld) .or. proposal_ld >= 0.5_dp*huge(1.0_dp)) then
                     res%status = 6
                     return
                  end if
                  if (green >= 0.0_dp) then
                     accepted = .true.
                  else
                     u = rand_uniform()
                     accepted = u < exp(green)
                  end if
               end if
               if (dbg) then
                  res%current(iter,:) = state
                  res%proposal(iter,:) = proposal
                  res%log_green(iter) = green
                  res%u(iter) = u
                  res%z(iter,:) = z
                  res%debug_accept(iter) = accepted
               end if
               if (accepted) then
                  state = proposal
                  current_ld = proposal_ld
                  accept_count = accept_count+1.0_dp
                  abatch = abatch+1.0_dp
               end if
               tries = tries+1.0_dp
               tbatch = tbatch+1.0_dp
            end do
            call call_out(outfun,state,obuf,data)
            bsum = bsum+obuf
         end do
         res%batch(ib,:) = bsum/real(bb,dp)
         res%accept_batch(ib) = abatch/tbatch
      end do
      res%accept = accept_count/tries
      res%final = state
   end function metrop

end module mcmc_metrop
