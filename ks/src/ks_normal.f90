! SPDX-License-Identifier: GPL-2.0-only
module ks_normal
   use ks_kinds, only: dp, pi, sqrt2, log2pi
   use ks_linalg, only: spd_inverse, spd_logdet, spd_cholesky
   use ks_rng, only: rng_state, rng_uniform, rng_normal, rng_chisq
   implicit none
   private
   public :: normal_pdf, normal_cdf, normal_quantile, normal_derivative
   public :: mvn_logpdf, mvn_pdf, mvn_derivative_tensor, mvn_sample
   public :: student_t_pdf, mvt_logpdf, mvt_pdf, mvt_sample
   public :: dnorm_mixture, dmvnorm_mixture, dmvt_mixture
   public :: rnorm_mixture, rmvnorm_mixture, rmvt_mixture
   public :: psins_1d
contains
   elemental function normal_pdf(x,mu,sigma) result(v)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: mu,sigma
      real(dp) :: v,m,s
      m=0.0_dp; s=1.0_dp
      if (present(mu)) m=mu
      if (present(sigma)) s=sigma
      if (s<=0.0_dp) then
         v=0.0_dp
      else
         v=exp(-0.5_dp*((x-m)/s)**2)/(sqrt(2.0_dp*pi)*s)
      end if
   end function normal_pdf

   elemental function normal_cdf(x,mu,sigma) result(v)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: mu,sigma
      real(dp) :: v,m,s
      m=0.0_dp; s=1.0_dp
      if (present(mu)) m=mu
      if (present(sigma)) s=sigma
      if (s<=0.0_dp) then
         if (x<m) then; v=0.0_dp; else; v=1.0_dp; end if
      else
         v=0.5_dp*erfc(-(x-m)/(s*sqrt2))
      end if
   end function normal_cdf

   elemental function normal_quantile(p,mu,sigma) result(x)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: mu,sigma
      real(dp) :: x,m,s,q,r,z
      real(dp), parameter :: a1=-3.969683028665376e1_dp,a2=2.209460984245205e2_dp
      real(dp), parameter :: a3=-2.759285104469687e2_dp,a4=1.383577518672690e2_dp
      real(dp), parameter :: a5=-3.066479806614716e1_dp,a6=2.506628277459239_dp
      real(dp), parameter :: b1=-5.447609879822406e1_dp,b2=1.615858368580409e2_dp
      real(dp), parameter :: b3=-1.556989798598866e2_dp,b4=6.680131188771972e1_dp,b5=-1.328068155288572e1_dp
      real(dp), parameter :: c1=-7.784894002430293e-3_dp,c2=-3.223964580411365e-1_dp
      real(dp), parameter :: c3=-2.400758277161838_dp,c4=-2.549732539343734_dp
      real(dp), parameter :: c5=4.374664141464968_dp,c6=2.938163982698783_dp
      real(dp), parameter :: d1=7.784695709041462e-3_dp,d2=3.224671290700398e-1_dp
      real(dp), parameter :: d3=2.445134137142996_dp,d4=3.754408661907416_dp
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
      m=0.0_dp; s=1.0_dp
      if (present(mu)) m=mu
      if (present(sigma)) s=sigma
      if (p<=0.0_dp) then; x=-huge(1.0_dp); return; end if
      if (p>=1.0_dp) then; x=huge(1.0_dp); return; end if
      if (p<plow) then
         q=sqrt(-2.0_dp*log(p))
         z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p<=phigh) then
         q=p-0.5_dp; r=q*q
         z=((((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q)/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      else
         q=sqrt(-2.0_dp*log(1.0_dp-p))
         z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      end if
      ! one Halley correction
      q=normal_cdf(z)-p
      z=z-q/(normal_pdf(z)+0.5_dp*z*q)
      x=m+s*z
   end function normal_quantile

   elemental function normal_derivative(x,mu,sigma,order) result(v)
      real(dp), intent(in) :: x,mu,sigma
      integer, intent(in) :: order
      real(dp) :: v,arg,h0,h1,hnew
      integer :: i
      if (sigma<=0.0_dp .or. order<0) then; v=0.0_dp; return; end if
      arg=(x-mu)/sigma
      h0=1.0_dp
      if (order==0) then
         hnew=h0
      else
         h1=arg
         hnew=h1
         do i=2,order
            hnew=arg*h1-real(i-1,dp)*h0
            h0=h1; h1=hnew
         end do
      end if
      v=(-1.0_dp)**order*normal_pdf(x,mu,sigma)*hnew/sigma**order
   end function normal_derivative

   function mvn_logpdf(x,mu,sigma) result(v)
      real(dp), intent(in) :: x(:),mu(:),sigma(:,:)
      real(dp) :: v,ainv(size(x),size(x)),dx(size(x)),ld
      integer :: info
      call spd_inverse(sigma,ainv,info)
      if (info/=0) then; v=-huge(1.0_dp); return; end if
      ld=spd_logdet(sigma,info)
      if (info/=0) then; v=-huge(1.0_dp); return; end if
      dx=x-mu
      v=-0.5_dp*(real(size(x),dp)*log2pi+ld+dot_product(dx,matmul(ainv,dx)))
   end function mvn_logpdf

   function mvn_pdf(x,mu,sigma) result(v)
      real(dp), intent(in) :: x(:),mu(:),sigma(:,:)
      real(dp) :: v,lv
      lv=mvn_logpdf(x,mu,sigma)
      if (lv<-log(huge(1.0_dp))) then; v=0.0_dp; else; v=exp(lv); end if
   end function mvn_pdf

   recursive function hermite_poly(indices,u,a) result(v)
      integer, intent(in) :: indices(:)
      real(dp), intent(in) :: u(:),a(:,:)
      real(dp) :: v
      integer :: n,k,i
      integer, allocatable :: rest(:),sub(:)
      n=size(indices)
      if (n==0) then; v=1.0_dp; return; end if
      i=indices(1)
      if (n==1) then; v=u(i); return; end if
      allocate(rest(n-1)); rest=indices(2:n)
      v=u(i)*hermite_poly(rest,u,a)
      do k=1,n-1
         allocate(sub(n-2))
         if (k>1) sub(1:k-1)=rest(1:k-1)
         if (k<n-1) sub(k:n-2)=rest(k+1:n-1)
         v=v-a(i,rest(k))*hermite_poly(sub,u,a)
         deallocate(sub)
      end do
   end function hermite_poly

   subroutine decode_tuple(code,d,r,idx)
      integer, intent(in) :: code,d,r
      integer, intent(out) :: idx(r)
      integer :: c,k
      c=code
      do k=1,r
         idx(k)=mod(c,d)+1
         c=c/d
      end do
   end subroutine decode_tuple

   subroutine mvn_derivative_tensor(x,mu,sigma,order,deriv,info)
      real(dp), intent(in) :: x(:),mu(:),sigma(:,:)
      integer, intent(in) :: order
      real(dp), allocatable, intent(out) :: deriv(:)
      integer, intent(out), optional :: info
      integer :: d,nout,j,ierr
      real(dp) :: ainv(size(x),size(x)),u(size(x)),phi
      integer, allocatable :: idx(:)
      d=size(x)
      if (order<0) then
         allocate(deriv(0)); if(present(info)) info=-1; return
      end if
      nout=d**order
      allocate(deriv(nout),idx(order))
      call spd_inverse(sigma,ainv,ierr)
      if (ierr/=0) then; deriv=0.0_dp; if(present(info)) info=ierr; return; end if
      u=matmul(ainv,x-mu)
      phi=mvn_pdf(x,mu,sigma)
      if (order==0) then
         deriv(1)=phi
      else
         do j=0,nout-1
            call decode_tuple(j,d,order,idx)
            deriv(j+1)=(-1.0_dp)**order*phi*hermite_poly(idx,u,ainv)
         end do
      end if
      if(present(info)) info=0
   end subroutine mvn_derivative_tensor

   subroutine mvn_sample(rng,mu,sigma,x,info)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: mu(:),sigma(:,:)
      real(dp), intent(out) :: x(size(mu))
      integer, intent(out), optional :: info
      real(dp) :: l(size(mu),size(mu)),z(size(mu))
      integer :: ierr,i
      call spd_cholesky(sigma,l,ierr)
      if (ierr/=0) then; x=mu; if(present(info)) info=ierr; return; end if
      do i=1,size(mu); z(i)=rng_normal(rng); end do
      x=mu+matmul(l,z)
      if(present(info)) info=0
   end subroutine mvn_sample

   elemental function student_t_pdf(x,mu,sigma,df) result(v)
      real(dp), intent(in) :: x,mu,sigma,df
      real(dp) :: v,z
      if (sigma<=0.0_dp .or. df<=0.0_dp) then; v=0.0_dp; return; end if
      z=(x-mu)/sigma
      v=exp(log_gamma(0.5_dp*(df+1.0_dp))-log_gamma(0.5_dp*df))/ &
        (sqrt(df*pi)*sigma)*(1.0_dp+z*z/df)**(-0.5_dp*(df+1.0_dp))
   end function student_t_pdf

   function mvt_logpdf(x,mu,sigma,df) result(v)
      real(dp), intent(in) :: x(:),mu(:),sigma(:,:),df
      real(dp) :: v,ainv(size(x),size(x)),dx(size(x)),ld,q,dreal
      integer :: info
      if (df<=0.0_dp) then; v=-huge(1.0_dp); return; end if
      call spd_inverse(sigma,ainv,info); if(info/=0) then; v=-huge(1.0_dp); return; end if
      ld=spd_logdet(sigma,info); if(info/=0) then; v=-huge(1.0_dp); return; end if
      dx=x-mu; q=dot_product(dx,matmul(ainv,dx)); dreal=real(size(x),dp)
      v=log_gamma(0.5_dp*(df+dreal))-log_gamma(0.5_dp*df)-0.5_dp*(dreal*log(df*pi)+ld) &
        -0.5_dp*(df+dreal)*log(1.0_dp+q/df)
   end function mvt_logpdf

   function mvt_pdf(x,mu,sigma,df) result(v)
      real(dp), intent(in) :: x(:),mu(:),sigma(:,:),df
      real(dp) :: v
      v=exp(mvt_logpdf(x,mu,sigma,df))
   end function mvt_pdf

   subroutine mvt_sample(rng,mu,sigma,df,x,info)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: mu(:),sigma(:,:),df
      real(dp), intent(out) :: x(size(mu))
      integer, intent(out), optional :: info
      real(dp) :: z(size(mu)),q
      integer :: ierr
      call mvn_sample(rng,0.0_dp*mu,sigma,z,ierr)
      if (ierr/=0 .or. df<=0.0_dp) then; x=mu; if(present(info)) info=-1; return; end if
      q=rng_chisq(rng,df)
      x=mu+z*sqrt(df/q)
      if(present(info)) info=0
   end subroutine mvt_sample

   function dnorm_mixture(x,mus,sigmas,props) result(v)
      real(dp), intent(in) :: x,mus(:),sigmas(:),props(:)
      real(dp) :: v
      integer :: k
      v=0.0_dp
      do k=1,size(props); v=v+props(k)*normal_pdf(x,mus(k),sigmas(k)); end do
   end function dnorm_mixture

   function dmvnorm_mixture(x,mus,sigmas,props) result(v)
      real(dp), intent(in) :: x(:),mus(:,:),sigmas(:,:,:),props(:)
      real(dp) :: v
      integer :: k
      v=0.0_dp
      do k=1,size(props); v=v+props(k)*mvn_pdf(x,mus(k,:),sigmas(:,:,k)); end do
   end function dmvnorm_mixture

   function dmvt_mixture(x,mus,sigmas,dfs,props) result(v)
      real(dp), intent(in) :: x(:),mus(:,:),sigmas(:,:,:),dfs(:),props(:)
      real(dp) :: v
      integer :: k
      v=0.0_dp
      do k=1,size(props); v=v+props(k)*mvt_pdf(x,mus(k,:),sigmas(:,:,k),dfs(k)); end do
   end function dmvt_mixture

   integer function draw_component(rng,props) result(k)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: props(:)
      real(dp) :: u,c,s
      integer :: j
      s=sum(props); u=rng_uniform(rng)*s; c=0.0_dp
      k=size(props)
      do j=1,size(props)
         c=c+props(j)
         if (u<=c) then; k=j; exit; end if
      end do
   end function draw_component

   subroutine rnorm_mixture(rng,mus,sigmas,props,x,labels)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: mus(:),sigmas(:),props(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out), optional :: labels(size(x))
      integer :: i,k
      do i=1,size(x)
         k=draw_component(rng,props)
         x(i)=mus(k)+sigmas(k)*rng_normal(rng)
         if(present(labels)) labels(i)=k
      end do
   end subroutine rnorm_mixture

   subroutine rmvnorm_mixture(rng,mus,sigmas,props,x,labels)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: mus(:,:),sigmas(:,:,:),props(:)
      real(dp), intent(out) :: x(:,:)
      integer, intent(out), optional :: labels(size(x,1))
      integer :: i,k
      do i=1,size(x,1)
         k=draw_component(rng,props)
         call mvn_sample(rng,mus(k,:),sigmas(:,:,k),x(i,:))
         if(present(labels)) labels(i)=k
      end do
   end subroutine rmvnorm_mixture

   subroutine rmvt_mixture(rng,mus,sigmas,dfs,props,x,labels)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: mus(:,:),sigmas(:,:,:),dfs(:),props(:)
      real(dp), intent(out) :: x(:,:)
      integer, intent(out), optional :: labels(size(x,1))
      integer :: i,k
      do i=1,size(x,1)
         k=draw_component(rng,props)
         call mvt_sample(rng,mus(k,:),sigmas(:,:,k),dfs(k),x(i,:))
         if(present(labels)) labels(i)=k
      end do
   end subroutine rmvt_mixture

   elemental function psins_1d(r,sigma) result(v)
      integer, intent(in) :: r
      real(dp), intent(in) :: sigma
      real(dp) :: v
      if (mod(r,2)==1 .or. sigma<=0.0_dp) then
         v=0.0_dp
      else
         v=(-1.0_dp)**(r/2)*exp(log_gamma(real(r+1,dp))-log_gamma(real(r/2+1,dp))) &
           /((2.0_dp*sigma)**(r+1)*sqrt(pi))
      end if
   end function psins_1d
end module ks_normal
