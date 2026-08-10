! SPDX-License-Identifier: GPL-3.0-only
! Cone projections translated from upstream SCS 3.x (MIT license).
module scs_cones
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
   use scs_kinds, only : dp, i4
   use scs_types, only : scs_cone, scs_scaling
   use scs_linalg, only : norm_2, norm_diff
   implicit none
   private
   public :: validate_cone, cone_dimension, set_r_y, project_dual_cone, project_primal_cone
   public :: enforce_cone_boundaries, scale_box_cone

   integer, parameter :: box_cone_max_iters = 25
   integer, parameter :: pow_cone_max_iters = 20
   real(dp), parameter :: pow_cone_tol = 1.0e-9_dp
   real(dp), parameter :: max_box_val = 1.0e15_dp
   real(dp), parameter :: exp_cone_infinity = 1.0e15_dp

contains

   pure integer(i4) function sd_cone_size(n) result(sz)
      integer(i4), intent(in) :: n
      sz = n * (n + 1_i4) / 2_i4
   end function sd_cone_size

   pure integer(i4) function cone_dimension(k) result(m)
      type(scs_cone), intent(in) :: k
      integer :: j
      m = k%z + k%l + k%bsize + 3_i4*(k%ep+k%ed)
      if (allocated(k%q)) m = m + sum(k%q)
      if (allocated(k%s)) then
         do j = 1, size(k%s)
            m = m + sd_cone_size(k%s(j))
         end do
      end if
      if (allocated(k%p)) m = m + 3_i4*int(size(k%p),i4)
   end function cone_dimension

   logical function validate_cone(k, m) result(ok)
      type(scs_cone), intent(in) :: k
      integer(i4), intent(in) :: m
      ok = .false.
      if (k%z < 0 .or. k%l < 0 .or. k%bsize < 0 .or. k%ep < 0 .or. k%ed < 0) return
      if (allocated(k%q)) then
         if (any(k%q < 0)) return
      end if
      if (allocated(k%s)) then
         if (any(k%s < 0)) return
      end if
      if (allocated(k%p)) then
         if (any(k%p < -1.0_dp) .or. any(k%p > 1.0_dp)) return
      end if
      if (k%bsize > 1) then
         if (.not. allocated(k%bl) .or. .not. allocated(k%bu)) return
         if (size(k%bl) /= k%bsize-1 .or. size(k%bu) /= k%bsize-1) return
         if (any(k%bl > k%bu)) return
      end if
      if (cone_dimension(k) /= m) return
      ok = .true.
   end function validate_cone

   subroutine set_r_y(k, scale, r_y)
      type(scs_cone), intent(in) :: k
      real(dp), intent(in) :: scale
      real(dp), intent(out) :: r_y(:)
      if (k%z > 0) r_y(1:k%z) = 1.0_dp/(1000.0_dp*scale)
      if (size(r_y) > k%z) r_y(k%z+1:) = 1.0_dp/scale
   end subroutine set_r_y

   subroutine enforce_cone_boundaries(k, vec, use_mean)
      type(scs_cone), intent(in) :: k
      real(dp), intent(inout) :: vec(:)
      logical, intent(in) :: use_mean
      integer :: count, j, nblk
      real(dp) :: w
      count = k%z + k%l + k%bsize
      if (allocated(k%q)) then
         do j=1,size(k%q)
            nblk=k%q(j); if (nblk>0) then
               if (use_mean) then
                  w = sum(vec(count+1:count+nblk)) / real(nblk,dp)
               else
                  w = maxval(abs(vec(count+1:count+nblk)))
               end if
               vec(count+1:count+nblk)=w
            end if
            count=count+nblk
         end do
      end if
      if (allocated(k%s)) then
         do j=1,size(k%s)
            nblk=sd_cone_size(k%s(j)); if (nblk>0) then
               if (use_mean) then
                  w = sum(vec(count+1:count+nblk)) / real(nblk,dp)
               else
                  w = maxval(abs(vec(count+1:count+nblk)))
               end if
               vec(count+1:count+nblk)=w
            end if
            count=count+nblk
         end do
      end if
      do j=1,k%ep+k%ed
         nblk=3; if (use_mean) then; w=sum(vec(count+1:count+3))/3.0_dp; else; w=maxval(abs(vec(count+1:count+3))); end if
         vec(count+1:count+3)=w; count=count+3
      end do
      if (allocated(k%p)) then
         do j=1,size(k%p)
            if (use_mean) then; w=sum(vec(count+1:count+3))/3.0_dp; else; w=maxval(abs(vec(count+1:count+3))); end if
            vec(count+1:count+3)=w; count=count+3
         end do
      end if
   end subroutine enforce_cone_boundaries

   subroutine scale_box_cone(k, scal)
      type(scs_cone), intent(inout) :: k
      type(scs_scaling), intent(in), optional :: scal
      integer :: j, off
      real(dp) :: d0, dj
      if (k%bsize <= 1) return
      off = k%z + k%l
      do j=1,k%bsize-1
         d0=1.0_dp; dj=1.0_dp
         if (present(scal)) then
            d0=scal%D(off+1); dj=scal%D(off+j+1)
         end if
         if (k%bu(j) >= max_box_val) then
            k%bu(j)=ieee_value(1.0_dp,ieee_positive_inf)
         else
            k%bu(j)=dj*k%bu(j)/d0
         end if
         if (k%bl(j) <= -max_box_val) then
            k%bl(j)=ieee_value(1.0_dp,ieee_negative_inf)
         else
            k%bl(j)=dj*k%bl(j)/d0
         end if
      end do
   end subroutine scale_box_cone

   subroutine project_dual_cone(x, k, r_y, box_t_warm)
      real(dp), intent(inout) :: x(:)
      type(scs_cone), intent(in) :: k
      real(dp), intent(in), optional :: r_y(:)
      real(dp), intent(inout), optional :: box_t_warm
      real(dp), allocatable :: orig(:), z(:)
      integer :: i
      allocate(orig(size(x)), z(size(x)))
      orig=x
      if (present(r_y)) then
         z=-r_y*x
         call project_primal_cone(z,k,r_y,box_t_warm)
         do i=1,size(x)
            x(i)=orig(i)+z(i)/r_y(i)
         end do
      else
         z=-x
         call project_primal_cone(z,k,box_t_warm=box_t_warm)
         x=orig+z
      end if
   end subroutine project_dual_cone

   subroutine project_primal_cone(x,k,r_y,box_t_warm)
      real(dp), intent(inout) :: x(:)
      type(scs_cone), intent(in) :: k
      real(dp), intent(in), optional :: r_y(:)
      real(dp), intent(inout), optional :: box_t_warm
      integer :: count,j,nblk
      real(dp) :: tw
      count=0
      if (k%z>0) then; x(1:k%z)=0.0_dp; count=count+k%z; end if
      if (k%l>0) then; x(count+1:count+k%l)=max(x(count+1:count+k%l),0.0_dp); count=count+k%l; end if
      if (k%bsize>0) then
         tw=1.0_dp; if (present(box_t_warm)) tw=box_t_warm
         if (present(r_y)) then
            call proj_box(x(count+1:count+k%bsize),k%bl,k%bu,tw,r_y(count+1:count+k%bsize))
         else
            call proj_box(x(count+1:count+k%bsize),k%bl,k%bu,tw)
         end if
         if (present(box_t_warm)) box_t_warm=tw
         count=count+k%bsize
      end if
      if (allocated(k%q)) then
         do j=1,size(k%q); nblk=k%q(j); call proj_soc(x(count+1:count+nblk)); count=count+nblk; end do
      end if
      if (allocated(k%s)) then
         do j=1,size(k%s); nblk=sd_cone_size(k%s(j)); call proj_psd(x(count+1:count+nblk),k%s(j)); count=count+nblk; end do
      end if
      do j=1,k%ep
         call proj_exp(x(count+1:count+3),.true.); count=count+3
      end do
      do j=1,k%ed
         call proj_exp(x(count+1:count+3),.false.); count=count+3
      end do
      if (allocated(k%p)) then
         do j=1,size(k%p)
            if (k%p(j)>=0.0_dp) then
               call proj_power(x(count+1:count+3),k%p(j))
            else
               block
                  real(dp) :: v(3)
                  v=-x(count+1:count+3); call proj_power(v,-k%p(j)); x(count+1:count+3)=x(count+1:count+3)+v
               end block
            end if
            count=count+3
         end do
      end if
   end subroutine project_primal_cone

   subroutine proj_soc(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: v1,s,alpha
      if (size(x)==0) return
      if (size(x)==1) then; x(1)=max(x(1),0.0_dp); return; end if
      v1=x(1); s=norm_2(x(2:)); alpha=(s+v1)/2.0_dp
      if (s<=v1) then
         return
      else if (s<=-v1) then
         x=0.0_dp
      else
         x(1)=alpha; x(2:)=x(2:)*(alpha/s)
      end if
   end subroutine proj_soc

   subroutine proj_box(tx,bl,bu,t,r_box)
      real(dp), intent(inout) :: tx(:)
      real(dp), intent(in) :: bl(:),bu(:)
      real(dp), intent(inout) :: t
      real(dp), intent(in), optional :: r_box(:)
      real(dp) :: gt,ht,tprev,rho_t,r
      integer :: iter,j
      if (size(tx)==1) then; tx(1)=max(tx(1),0.0_dp); t=tx(1); return; end if
      rho_t=1.0_dp; if (present(r_box)) rho_t=1.0_dp/r_box(1)
      do iter=1,box_cone_max_iters
         tprev=t; gt=rho_t*(t-tx(1)); ht=rho_t
         do j=1,size(tx)-1
            r=1.0_dp; if (present(r_box)) r=1.0_dp/r_box(j+1)
            if (tx(j+1)>t*bu(j)) then
               gt=gt+r*(t*bu(j)-tx(j+1))*bu(j); ht=ht+r*bu(j)*bu(j)
            else if (tx(j+1)<t*bl(j)) then
               gt=gt+r*(t*bl(j)-tx(j+1))*bl(j); ht=ht+r*bl(j)*bl(j)
            end if
         end do
         t=max(t-gt/max(ht,1.0e-8_dp),0.0_dp)
         if (abs(gt/max(ht,1.0e-6_dp))<1.0e-12_dp*max(t,1.0_dp) .or. abs(t-tprev)<1.0e-11_dp*max(t,1.0_dp)) exit
      end do
      do j=1,size(tx)-1
         tx(j+1)=min(max(tx(j+1),t*bl(j)),t*bu(j))
      end do
      tx(1)=t
   end subroutine proj_box

   subroutine proj_power(v,a)
      real(dp), intent(inout) :: v(3)
      real(dp), intent(in) :: a
      real(dp) :: xh,yh,rh,x,y,r,f,fp,dxdr,dydr
      integer :: it
      xh=v(1); yh=v(2); rh=abs(v(3))
      if (xh>=0 .and. yh>=0 .and. pow_cone_tol+xh**a*yh**(1.0_dp-a)>=rh) return
      if (xh<=0 .and. yh<=0 .and. pow_cone_tol+(-xh)**a*(-yh)**(1.0_dp-a)>=rh*a**a*(1.0_dp-a)**(1.0_dp-a)) then
         v=0.0_dp; return
      end if
      r=rh/2.0_dp; x=0.0_dp; y=0.0_dp
      do it=1,pow_cone_max_iters
         x=max(0.5_dp*(xh+sqrt(max(0.0_dp,xh*xh+4.0_dp*a*(rh-r)*r))),1.0e-12_dp)
         y=max(0.5_dp*(yh+sqrt(max(0.0_dp,yh*yh+4.0_dp*(1.0_dp-a)*(rh-r)*r))),1.0e-12_dp)
         f=x**a*y**(1.0_dp-a)-r
         if (abs(f)<pow_cone_tol) exit
         dxdr=a*(rh-2.0_dp*r)/(2.0_dp*x-xh)
         dydr=(1.0_dp-a)*(rh-2.0_dp*r)/(2.0_dp*y-yh)
         fp=x**a*y**(1.0_dp-a)*(a*dxdr/x+(1.0_dp-a)*dydr/y)-1.0_dp
         r=min(max(r-f/fp,0.0_dp),rh)
      end do
      v(1)=x; v(2)=y; v(3)=merge(-r,r,v(3)<0.0_dp)
   end subroutine proj_power

   subroutine proj_psd(x,n)
      real(dp), intent(inout) :: x(:)
      integer(i4), intent(in) :: n
      real(dp), allocatable :: a(:,:),eig(:),v(:,:)
      real(dp), parameter :: sqrt2=sqrt(2.0_dp)
      integer :: c,r,k,j
      if (n==0) return
      if (n==1) then; x(1)=max(x(1),0.0_dp); return; end if
      allocate(a(n,n),eig(n),v(n,n)); a=0.0_dp; k=0
      do c=1,n
         do r=c,n
            k=k+1
            if (r==c) then
               a(r,c)=x(k)
            else
               a(r,c)=x(k)/sqrt2; a(c,r)=a(r,c)
            end if
         end do
      end do
      call jacobi_eigen(a,eig,v)
      a=0.0_dp
      do j=1,n
         if (eig(j)>0.0_dp) then
            do c=1,n
               do r=1,n
                  a(r,c)=a(r,c)+eig(j)*v(r,j)*v(c,j)
               end do
            end do
         end if
      end do
      k=0
      do c=1,n
         do r=c,n
            k=k+1
            if (r==c) then; x(k)=a(r,c); else; x(k)=sqrt2*a(r,c); end if
         end do
      end do
   end subroutine proj_psd

   subroutine jacobi_eigen(a,eig,v)
      real(dp), intent(inout) :: a(:,:)
      real(dp), intent(out) :: eig(:),v(:,:)
      integer :: n,p,q,i,it,maxit
      real(dp) :: apq,app,aqq,tau,t,c,s,tmp,off
      n=size(a,1); v=0.0_dp
      do i=1,n; v(i,i)=1.0_dp; end do
      maxit=max(50,50*n*n)
      do it=1,maxit
         off=0.0_dp; p=1; q=min(2,n)
         do i=1,n-1
            block
               integer :: j
               do j=i+1,n
                  if (abs(a(i,j))>off) then; off=abs(a(i,j)); p=i; q=j; end if
               end do
            end block
         end do
         if (off<=1.0e-13_dp*max(1.0_dp,maxval(abs(a)))) exit
         apq=a(p,q); app=a(p,p); aqq=a(q,q); tau=(aqq-app)/(2.0_dp*apq)
         if (tau>=0.0_dp) then; t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau)); else; t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau)); end if
         c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
         do i=1,n
            if (i/=p .and. i/=q) then
               tmp=a(i,p); a(i,p)=c*tmp-s*a(i,q); a(p,i)=a(i,p); a(i,q)=s*tmp+c*a(i,q); a(q,i)=a(i,q)
            end if
         end do
         a(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
         a(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
         a(p,q)=0.0_dp; a(q,p)=0.0_dp
         do i=1,n
            tmp=v(i,p); v(i,p)=c*tmp-s*v(i,q); v(i,q)=s*tmp+c*v(i,q)
         end do
      end do
      do i=1,n; eig(i)=a(i,i); end do
   end subroutine jacobi_eigen

   pure logical function exp_finite(x)
      real(dp),intent(in)::x
      exp_finite=abs(x)<exp_cone_infinity
   end function exp_finite

   subroutine hfun(v0,rho,f,df)
      real(dp),intent(in)::v0(3),rho
      real(dp),intent(out)::f,df
      real(dp)::t0,s0,r0,ep,en
      r0=v0(1); s0=v0(2); t0=v0(3); ep=exp(rho); en=exp(-rho)
      f=((rho-1.0_dp)*r0+s0)*ep-(r0-rho*s0)*en-(rho*(rho-1.0_dp)+1.0_dp)*t0
      df=(rho*r0+s0)*ep+(r0-(rho-1.0_dp)*s0)*en-(2.0_dp*rho-1.0_dp)*t0
   end subroutine hfun

   real(dp) function root_binary(v0,xl0,xu0,x0) result(xp)
      real(dp),intent(in)::v0(3),xl0,xu0,x0
      real(dp)::xl,xu,x,f,df
      integer::i
      xl=xl0; xu=xu0; x=x0; xp=x
      do i=1,40
         call hfun(v0,x,f,df)
         if(f<0) then; xl=x; else; xu=x; end if
         xp=0.5_dp*(xl+xu)
         if(abs(xp-x)<=1.0e-12_dp*max(1.0_dp,abs(xp)) .or. &
            abs(xp-xl)<=epsilon(1.0_dp)*max(1.0_dp,abs(xp)) .or. &
            abs(xp-xu)<=epsilon(1.0_dp)*max(1.0_dp,abs(xp))) exit
         x=xp
      end do
   end function root_binary

   real(dp) function root_newton(v0,xl0,xu0,x0) result(root)
      real(dp),intent(in)::v0(3),xl0,xu0,x0
      real(dp)::xl,xu,x,xp,f,df
      integer::i
      xl=xl0; xu=xu0; x=x0
      do i=1,20
         call hfun(v0,x,f,df)
         if(abs(f)<=1.0e-15_dp) then; root=max(xl,min(xu,x)); return; end if
         if(f<0) then; xl=x; else; xu=x; end if
         if(xu<=xl) exit
         if(.not.exp_finite(f) .or. df<1.0e-13_dp) exit
         xp=x-f/df
         if(abs(xp-x)<=1.0e-15_dp*max(1.0_dp,abs(xp))) then; root=max(xl,min(xu,x)); return; end if
         if(xp>=xu) then; x=min(0.05_dp*x+0.95_dp*xu,xu); else if(xp<=xl) then; x=max(0.05_dp*x+0.95_dp*xl,xl); else; x=xp; end if
      end do
      root=root_binary(v0,xl,xu,x)
   end function root_newton

   real(dp) function ppsi(v0) result(v)
      real(dp),intent(in)::v0(3); real(dp)::psi,r0,s0,d
      r0=v0(1); s0=v0(2); d=sqrt(r0*r0+s0*s0-r0*s0)
      if(r0>s0) then; psi=(r0-s0+d)/r0; else; psi=-s0/(r0-s0-d); end if
      v=((psi-1.0_dp)*r0+s0)/(psi*(psi-1.0_dp)+1.0_dp)
   end function ppsi
   real(dp) function pomega(rho) result(v)
      real(dp),intent(in)::rho; v=exp(rho)/(rho*(rho-1.0_dp)+1.0_dp); if(rho<2.0_dp)v=min(v,exp(2.0_dp)/3.0_dp)
   end function pomega
   real(dp) function dpsi(v0) result(v)
      real(dp),intent(in)::v0(3); real(dp)::psi,r0,s0,d
      r0=v0(1); s0=v0(2); d=sqrt(r0*r0+s0*s0-r0*s0)
      if(s0>r0) then; psi=(r0-d)/s0; else; psi=(r0-s0)/(r0+d); end if
      v=(r0-psi*s0)/(psi*(psi-1.0_dp)+1.0_dp)
   end function dpsi
   real(dp) function domega(rho) result(v)
      real(dp),intent(in)::rho; v=-exp(-rho)/(rho*(rho-1.0_dp)+1.0_dp); if(rho>-1.0_dp)v=max(v,-exp(1.0_dp)/3.0_dp)
   end function domega

   subroutine exp_heuristics(v0,vp,vd,pdist,ddist)
      real(dp),intent(in)::v0(3); real(dp),intent(out)::vp(3),vd(3),pdist,ddist
      real(dp)::r0,s0,t0,tp,td,nd
      r0=v0(1);s0=v0(2);t0=v0(3)
      vp=[min(r0,0.0_dp),0.0_dp,max(t0,0.0_dp)]; pdist=norm_diff(v0,vp)
      if(s0>0.0_dp) then; tp=max(t0,s0*exp(r0/s0)); nd=tp-t0; if(nd<pdist) then; vp=[r0,s0,tp]; pdist=nd; end if; end if
      vd=[0.0_dp,min(s0,0.0_dp),min(t0,0.0_dp)]; ddist=norm_diff(v0,vd)
      if(r0>0.0_dp) then; td=min(t0,-r0*exp(s0/r0-1.0_dp)); nd=t0-td; if(nd<ddist) then; vd=[r0,s0,td]; ddist=nd; end if; end if
   end subroutine exp_heuristics

   subroutine exp_bracket(v0,pdist,ddist,low,upr)
      real(dp),intent(in)::v0(3),pdist,ddist; real(dp),intent(out)::low,upr
      real(dp)::t0,s0,r0,baselow,baseupr,dpv,ddv,cur,tpu,tdl,fl,fu,df
      r0=v0(1);s0=v0(2);t0=v0(3); baselow=-exp_cone_infinity; baseupr=exp_cone_infinity; low=baselow;upr=baseupr
      dpv=sqrt(max(0.0_dp,pdist*pdist-min(s0,0.0_dp)**2)); ddv=sqrt(max(0.0_dp,ddist*ddist-min(r0,0.0_dp)**2))
      if(t0>0) then; cur=log(t0/ppsi(v0)); low=max(low,cur); else if(t0<0) then; cur=-log(-t0/dpsi(v0)); upr=min(upr,cur); end if
      if (r0 > 0.0_dp) then
         baselow = 1.0_dp - s0/r0
         low = max(low,baselow)
         tpu = max(1.0e-12_dp,min(ddv,dpv+t0))
         cur = max(low,baselow+tpu/r0/pomega(low))
         upr = min(upr,cur)
      end if
      if (s0 > 0.0_dp) then
         baseupr = r0/s0
         upr = min(upr,baseupr)
         tdl = -max(1.0e-12_dp,min(dpv,ddv-t0))
         cur = min(upr,baseupr-tdl/s0/domega(upr))
         low = max(low,cur)
      end if
      low=max(baselow,min(baseupr,min(low,upr))); upr=max(baselow,min(baseupr,max(low,upr)))
      if (abs(low-upr) > epsilon(1.0_dp)*max(1.0_dp,abs(low),abs(upr))) then
         call hfun(v0,low,fl,df)
         call hfun(v0,upr,fu,df)
         if (fl*fu > 0.0_dp) then
            if (abs(fl) < abs(fu)) then
               upr = low
            else
               low = upr
            end if
         end if
      end if
   end subroutine exp_bracket

   real(dp) function exp_sol_primal(v0,rho,vp) result(dist)
      real(dp),intent(in)::v0(3),rho;real(dp),intent(out)::vp(3);real(dp)::lin,e,q
      lin=(rho-1.0_dp)*v0(1)+v0(2);e=exp(rho)
      if (lin > 0.0_dp .and. exp_finite(e)) then
         q = rho*(rho-1.0_dp)+1.0_dp
         vp = [rho*lin/q,lin/q,e*lin/q]
         dist = norm_diff(vp,v0)
      else
         vp = [0.0_dp,0.0_dp,exp_cone_infinity]
         dist = exp_cone_infinity
      end if
   end function exp_sol_primal
   real(dp) function exp_sol_polar(v0,rho,vd) result(dist)
      real(dp),intent(in)::v0(3),rho;real(dp),intent(out)::vd(3);real(dp)::lin,e,q
      lin=v0(1)-rho*v0(2);e=exp(-rho)
      if (lin > 0.0_dp .and. exp_finite(e)) then
         q = rho*(rho-1.0_dp)+1.0_dp
         vd = [lin/q,(1.0_dp-rho)*lin/q,-e*lin/q]
         dist = norm_diff(v0,vd)
      else
         vd = [0.0_dp,0.0_dp,-exp_cone_infinity]
         dist = exp_cone_infinity
      end if
   end function exp_sol_polar

   subroutine proj_exp(v,primal)
      real(dp),intent(inout)::v(3);logical,intent(in)::primal
      real(dp)::v0(3),vp(3),vd(3),vh(3),pd,dd,err,xl,xh,rho,dh
      logical::opt
      v0=v; if(.not.primal)v0=-v0
      call exp_heuristics(v0,vp,vd,pd,dd)
      err = maxval(abs(vp+vd-v0))
      opt = (v0(2)<=0.0_dp .and. v0(1)<=0.0_dp) .or. min(pd,dd)<=1.0e-8_dp .or. &
         (err<=1.0e-8_dp .and. dot_product(vp,vd)<=1.0e-8_dp)
      if(.not.opt) then
         call exp_bracket(v0,pd,dd,xl,xh);rho=root_newton(v0,xl,xh,0.5_dp*(xl+xh))
         if(primal) then; dh=exp_sol_primal(v0,rho,vh);if(dh<=pd)vp=vh;else; dh=exp_sol_polar(v0,rho,vh);if(dh<=dd)vd=vh;end if
      end if
      if(primal) then;v=vp;else;v=-vd;end if
   end subroutine proj_exp

end module scs_cones
