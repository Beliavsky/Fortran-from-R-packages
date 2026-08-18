module mcmc_morph
   use mcmc_kinds, only : dp
   use mcmc_numerics, only : euclid_norm, rand_uniform, rand_normal_vec
   use mcmc_metrop, only : mcmc_scale, metrop_result, log_density_callback, output_callback, &
      SCALE_MODE_CONSTANT, SCALE_MODE_DIAGONAL, SCALE_MODE_FULL, scale_constant
   implicit none
   private

   type, public :: morph_transform
      logical :: use_subexponential = .false.
      logical :: use_exponential = .false.
      real(dp) :: b = 1.0_dp
      real(dp) :: r = 0.0_dp
      real(dp) :: power = 3.0_dp
      real(dp), allocatable :: center(:)
   contains
      procedure :: transform => morph_forward
      procedure :: inverse => morph_inverse
      procedure :: log_jacobian => morph_log_jacobian
   end type morph_transform

   public :: morph_identity, morph_create, morph_metrop

contains

   function morph_identity(d,center) result(m)
      integer, intent(in) :: d
      real(dp), intent(in), optional :: center(:)
      type(morph_transform) :: m
      allocate(m%center(d))
      m%center = 0.0_dp
      if (present(center)) call set_center(m,center)
   end function morph_identity

   function morph_create(d,b,r,power,center) result(m)
      integer, intent(in) :: d
      real(dp), intent(in), optional :: b,r,power,center(:)
      type(morph_transform) :: m
      allocate(m%center(d))
      m%center = 0.0_dp
      if (present(center)) call set_center(m,center)
      if (present(b)) then
         if (b <= 0.0_dp) error stop "morph_create: b must be positive"
         m%use_subexponential = .true.
         m%b = b
      end if
      if (present(r) .or. present(power)) then
         m%use_exponential = .true.
         m%r = 0.0_dp
         m%power = 3.0_dp
         if (present(r)) m%r = r
         if (present(power)) m%power = power
         if (m%r < 0.0_dp .or. m%power <= 2.0_dp) &
            error stop "morph_create: r>=0 and power>2 required"
      end if
   end function morph_create

   subroutine set_center(m,c)
      type(morph_transform), intent(inout) :: m
      real(dp), intent(in) :: c(:)
      if (size(c) == 1) then
         m%center = c(1)
      else if (size(c) == size(m%center)) then
         m%center = c
      else
         error stop "morph: center must be scalar or state length"
      end if
   end subroutine set_center

   pure real(dp) function expo_expand(x,r,p) result(y)
      real(dp), intent(in) :: x,r,p
      if (x <= r) then
         y = x
      else
         y = x+(x-r)**p
      end if
   end function expo_expand

   pure real(dp) function expo_expand_deriv(x,r,p) result(y)
      real(dp), intent(in) :: x,r,p
      if (x <= r) then
         y = 1.0_dp
      else
         y = 1.0_dp+p*(x-r)**(p-1.0_dp)
      end if
   end function expo_expand_deriv

   real(dp) function expo_contract(y,r,p) result(x)
      real(dp), intent(in) :: y,r,p
      real(dp) :: cur,err,der,next
      integer :: i
      if (y <= r) then
         x = y
         return
      end if
      cur = r+(max(y-r,0.0_dp))**(1.0_dp/p)
      cur = max(r,cur)
      do i = 1, 100
         err = expo_expand(cur,r,p)-y
         der = expo_expand_deriv(cur,r,p)
         next = cur-err/der
         if (next < r) next = 0.5_dp*(cur+r)
         if (abs(next-cur) <= 4.0_dp*epsilon(1.0_dp)*(1.0_dp+abs(cur))) exit
         cur = next
      end do
      x = next
   end function expo_contract

   pure real(dp) function sub_expand(x,b) result(y)
      real(dp), intent(in) :: x,b
      if (x > 1.0_dp/b) then
         y = exp(b*x)-exp(1.0_dp)/3.0_dp
      else
         y = (x*b)**3*exp(1.0_dp)/6.0_dp+x*b*exp(1.0_dp)/2.0_dp
      end if
   end function sub_expand

   pure real(dp) function sub_expand_deriv(x,b) result(y)
      real(dp), intent(in) :: x,b
      if (x > 1.0_dp/b) then
         y = b*exp(b*x)
      else
         y = b*(x*b)**2*exp(1.0_dp)/2.0_dp+b*exp(1.0_dp)/2.0_dp
      end if
   end function sub_expand_deriv

   real(dp) function sub_contract(y,b) result(x)
      real(dp), intent(in) :: y,b
      real(dp) :: cur,err,der,next
      integer :: i
      if (y > 2.0_dp*exp(1.0_dp)/3.0_dp) then
         x = log(y+exp(1.0_dp)/3.0_dp)/b
         return
      end if
      cur = min(1.0_dp/b,max(0.0_dp,2.0_dp*y/(b*exp(1.0_dp))))
      do i = 1, 100
         err = sub_expand(cur,b)-y
         der = sub_expand_deriv(cur,b)
         next = max(0.0_dp,cur-err/der)
         if (abs(next-cur) <= 4.0_dp*epsilon(1.0_dp)*(1.0_dp+abs(cur))) exit
         cur = next
      end do
      x = next
   end function sub_contract

   real(dp) function radial_contract(m,radius) result(y)
      class(morph_transform), intent(in) :: m
      real(dp), intent(in) :: radius
      y = radius
      if (m%use_subexponential) y = sub_contract(y,m%b)
      if (m%use_exponential) y = expo_contract(y,m%r,m%power)
   end function radial_contract

   real(dp) function radial_expand(m,radius) result(y)
      class(morph_transform), intent(in) :: m
      real(dp), intent(in) :: radius
      y = radius
      if (m%use_exponential) y = expo_expand(y,m%r,m%power)
      if (m%use_subexponential) y = sub_expand(y,m%b)
   end function radial_expand

   real(dp) function radial_expand_deriv(m,radius) result(y)
      class(morph_transform), intent(in) :: m
      real(dp), intent(in) :: radius
      real(dp) :: temp
      y = 1.0_dp
      temp = radius
      if (m%use_exponential) then
         y = y*expo_expand_deriv(temp,m%r,m%power)
         temp = expo_expand(temp,m%r,m%power)
      end if
      if (m%use_subexponential) y = y*sub_expand_deriv(temp,m%b)
   end function radial_expand_deriv

   subroutine isotropic_apply(x,radial,y)
      real(dp), intent(in) :: x(:),radial
      real(dp), intent(out) :: y(size(x))
      real(dp) :: nrm
      nrm = euclid_norm(x)
      if (nrm <= tiny(1.0_dp)) then
         y = 0.0_dp
      else
         y = radial*x/nrm
      end if
   end subroutine isotropic_apply

   subroutine morph_forward(self,x,y)
      class(morph_transform), intent(in) :: self
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(size(x))
      real(dp), allocatable :: centered(:)
      real(dp) :: nrm,rad
      if (size(x) /= size(self%center)) error stop "morph_forward: dimension mismatch"
      centered = x-self%center
      nrm = euclid_norm(centered)
      rad = radial_contract(self,nrm)
      call isotropic_apply(centered,rad,y)
   end subroutine morph_forward

   subroutine morph_inverse(self,y,x)
      class(morph_transform), intent(in) :: self
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: x(size(y))
      real(dp), allocatable :: temp(:)
      real(dp) :: nrm,rad
      if (size(y) /= size(self%center)) error stop "morph_inverse: dimension mismatch"
      allocate(temp(size(y)))
      nrm = euclid_norm(y)
      rad = radial_expand(self,nrm)
      call isotropic_apply(y,rad,temp)
      x = temp+self%center
   end subroutine morph_inverse

   real(dp) function morph_log_jacobian(self,y) result(lj)
      class(morph_transform), intent(in) :: self
      real(dp), intent(in) :: y(:)
      real(dp) :: nrm,rad,der
      integer :: k
      k = size(y)
      nrm = euclid_norm(y)
      rad = radial_expand(self,nrm)
      der = radial_expand_deriv(self,nrm)
      if (nrm <= tiny(1.0_dp)) then
         lj = real(k,dp)*log(der)
      else
         lj = log(der)+real(k-1,dp)*(log(rad)-log(nrm))
      end if
   end function morph_log_jacobian

   function morph_metrop(lud,morph,initial,nbatch,blen,nspac,scale,out_dim,outfun,data,debug) result(res)
      procedure(log_density_callback) :: lud
      type(morph_transform), intent(in) :: morph
      real(dp), intent(in) :: initial(:)
      integer, intent(in) :: nbatch
      integer, intent(in), optional :: blen,nspac,out_dim
      type(mcmc_scale), intent(in), optional :: scale
      procedure(output_callback), optional :: outfun
      class(*), intent(in), optional :: data
      logical, intent(in), optional :: debug
      type(metrop_result) :: res
      type(mcmc_scale) :: sc
      real(dp), allocatable :: state(:),proposal(:),z(:),orig(:),p_orig(:),obuf(:),bsum(:)
      real(dp) :: cur_ld,prop_ld,green,u
      real(dp) :: accept_count,tries,abatch,tbatch
      integer :: d,od,bb,ns,ib,jb,is,iter,niter,i,j
      logical :: accepted,dbg

      d=size(initial)
      bb=1; ns=1
      if (present(blen)) bb=blen
      if (present(nspac)) ns=nspac
      if (size(morph%center) /= d .or. nbatch<=0 .or. bb<=0 .or. ns<=0) then
         res%status=1; return
      end if
      sc=scale_constant(1.0_dp)
      if (present(scale)) sc=scale
      if (present(outfun)) then
         if (.not. present(out_dim)) then
            res%status=2; return
         end if
         od=out_dim
      else
         od=d
      end if
      dbg=.false.; if (present(debug)) dbg=debug
      allocate(state(d),proposal(d),z(d),orig(d),p_orig(d),obuf(od),bsum(od))
      call morph%transform(initial,state)
      call morph%inverse(state,orig)
      if (present(data)) then
         call lud(orig,cur_ld,data)
      else
         call lud(orig,cur_ld)
      end if
      cur_ld=cur_ld+morph%log_jacobian(state)
      if (.not. ieee_is_finite_local(cur_ld)) then
         res%status=3; return
      end if
      res%initial=initial
      allocate(res%batch(nbatch,od),res%accept_batch(nbatch))
      res%nbatch=nbatch; res%blen=bb; res%nspac=ns; res%scale=sc; res%debug=dbg
      niter=nbatch*bb*ns
      if (dbg) allocate(res%current(niter,d),res%proposal(niter,d),res%log_green(niter), &
                        res%u(niter),res%z(niter,d),res%debug_accept(niter))
      accept_count=0.0_dp; tries=0.0_dp; iter=0
      do ib=1,nbatch
         bsum=0.0_dp; abatch=0.0_dp; tbatch=0.0_dp
         do jb=1,bb
            do is=1,ns
               iter=iter+1
               call rand_normal_vec(z)
               select case(sc%mode)
               case(SCALE_MODE_CONSTANT)
                  proposal=state+sc%scalar*z
               case(SCALE_MODE_DIAGONAL)
                  if (.not. allocated(sc%diagonal) .or. size(sc%diagonal)/=d) then
                     res%status=4; return
                  end if
                  proposal=state+sc%diagonal*z
               case(SCALE_MODE_FULL)
                  if (.not. allocated(sc%full) .or. size(sc%full,1)/=d .or. size(sc%full,2)/=d) then
                     res%status=4; return
                  end if
                  proposal=state
                  do i=1,d
                     do j=1,d
                        proposal(j)=proposal(j)+sc%full(j,i)*z(i)
                     end do
                  end do
               end select
               call morph%inverse(proposal,p_orig)
               if (present(data)) then
                  call lud(p_orig,prop_ld,data)
               else
                  call lud(p_orig,prop_ld)
               end if
               if (prop_ld > -huge(1.0_dp)/2.0_dp) prop_ld=prop_ld+morph%log_jacobian(proposal)
               green=prop_ld-cur_ld
               accepted=.false.; u=-1.0_dp
               if (prop_ld > -huge(1.0_dp)/2.0_dp) then
                  if (green>=0.0_dp) then
                     accepted=.true.
                  else
                     u=rand_uniform(); accepted=u<exp(green)
                  end if
               end if
               if (dbg) then
                  call morph%inverse(state,orig)
                  call morph%inverse(proposal,p_orig)
                  res%current(iter,:)=orig
                  res%proposal(iter,:)=p_orig
                  res%log_green(iter)=green; res%u(iter)=u; res%z(iter,:)=z
                  res%debug_accept(iter)=accepted
               end if
               if (accepted) then
                  state=proposal; cur_ld=prop_ld; accept_count=accept_count+1.0_dp; abatch=abatch+1.0_dp
               end if
               tries=tries+1.0_dp; tbatch=tbatch+1.0_dp
            end do
            call morph%inverse(state,orig)
            if (present(outfun)) then
               if (present(data)) then
                  call outfun(orig,obuf,data)
               else
                  call outfun(orig,obuf)
               end if
            else
               obuf=orig
            end if
            bsum=bsum+obuf
         end do
         res%batch(ib,:)=bsum/real(bb,dp)
         res%accept_batch(ib)=abatch/tbatch
      end do
      res%accept=accept_count/tries
      call morph%inverse(state,orig)
      res%final=orig
   end function morph_metrop

   pure logical function ieee_is_finite_local(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) < huge(1.0_dp)
   end function ieee_is_finite_local

end module mcmc_morph
