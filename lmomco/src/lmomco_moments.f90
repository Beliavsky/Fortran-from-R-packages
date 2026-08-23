module lmomco_moments
   use lmomco_kinds, only : dp, pi
   use lmomco_types, only : lmomco_params, make_params
   use lmomco_distributions, only : lmomco_quantile
   implicit none
   private
   public :: sample_lmoments, pwm_unbiased, lmom_ratios, theoretical_lmoments
   public :: fit_lmoments

contains

   pure real(dp) function choose_real(n,k) result(c)
      integer, intent(in) :: n,k
      integer :: i, kk
      if(k<0 .or. k>n) then; c=0.0_dp; return; end if
      kk=min(k,n-k); c=1.0_dp
      do i=1,kk
         c=c*real(n-kk+i,dp)/real(i,dp)
      end do
   end function choose_real

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i,j
      real(dp) :: key
      do i=2,size(x)
         key=x(i); j=i-1
         do while(j>=1)
            if(x(j)<=key) exit
            x(j+1)=x(j); j=j-1
         end do
         x(j+1)=key
      end do
   end subroutine sort_real

   subroutine pwm_unbiased(x, rmax, b)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: rmax
      real(dp), intent(out) :: b(0:rmax)
      real(dp), allocatable :: s(:)
      integer :: n,r,k
      n=size(x); allocate(s(n)); s=x; call sort_real(s)
      b=0.0_dp
      do r=0,rmax
         if(n-1<r) cycle
         do k=r+1,n
            b(r)=b(r)+choose_real(k-1,r)/choose_real(n-1,r)*s(k)
         end do
         b(r)=b(r)/real(n,dp)
      end do
   end subroutine pwm_unbiased

   subroutine sample_lmoments(x, nmom, lmom)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: nmom
      real(dp), intent(out) :: lmom(nmom)
      real(dp), allocatable :: b(:)
      integer :: r,k
      if(nmom<1) return
      allocate(b(0:nmom-1)); call pwm_unbiased(x,nmom-1,b)
      do r=0,nmom-1
         lmom(r+1)=0.0_dp
         do k=0,r
            lmom(r+1)=lmom(r+1)+(-1.0_dp)**real(r-k,dp)*choose_real(r,k)*choose_real(r+k,k)*b(k)
         end do
      end do
   end subroutine sample_lmoments

   subroutine lmom_ratios(lmom, ratios)
      real(dp), intent(in) :: lmom(:)
      real(dp), intent(out) :: ratios(size(lmom))
      integer :: r
      ratios=lmom
      if(size(lmom)>=2 .and. abs(lmom(1))>tiny(1.0_dp)) ratios(2)=lmom(2)/lmom(1)
      if(size(lmom)>=3 .and. abs(lmom(2))>tiny(1.0_dp)) then
         do r=3,size(lmom)
            ratios(r)=lmom(r)/lmom(2)
         end do
      end if
   end subroutine lmom_ratios

   pure real(dp) function shifted_legendre(r, p) result(v)
      integer, intent(in) :: r
      real(dp), intent(in) :: p
      integer :: k
      if(r==0) then; v=1.0_dp; return; end if
      v=0.0_dp
      do k=0,r
         v=v+(-1.0_dp)**real(r-k,dp)*choose_real(r,k)*choose_real(r+k,k)*p**k
      end do
   end function shifted_legendre

   subroutine theoretical_lmoments(par, nmom, lmom, ngrid)
      type(lmomco_params), intent(in) :: par
      integer, intent(in) :: nmom
      real(dp), intent(out) :: lmom(nmom)
      integer, intent(in), optional :: ngrid
      integer :: n,i,r
      real(dp) :: h,p,w,q
      n=4000; if(present(ngrid)) n=max(200,ngrid); if(mod(n,2)/=0)n=n+1
      h=(1.0_dp-2.0e-8_dp)/real(n,dp); lmom=0.0_dp
      do i=0,n
         p=1.0e-8_dp+h*real(i,dp); q=lmomco_quantile(p,par)
         if(i==0 .or. i==n) then; w=1.0_dp
         else if(mod(i,2)==0) then; w=2.0_dp
         else; w=4.0_dp; end if
         do r=0,nmom-1
            lmom(r+1)=lmom(r+1)+w*q*shifted_legendre(r,p)
         end do
      end do
      lmom=lmom*h/3.0_dp
   end subroutine theoretical_lmoments

   function fit_lmoments(family, lmom) result(par)
      character(len=*), intent(in) :: family
      real(dp), intent(in) :: lmom(:)
      type(lmomco_params) :: par
      character(len=8) :: f
      real(dp) :: l1,l2,t3,c,k,a,u,gam1,euler
      f=adjustl(family); l1=lmom(1); l2=lmom(2); euler=0.5772156649015328606_dp
      select case(trim(f))
      case('nor')
         par=make_params('nor',[l1,l2*sqrt(pi)])
      case('exp')
         a=2.0_dp*l2; par=make_params('exp',[l1-a,a])
      case('gum')
         a=l2/log(2.0_dp); u=l1-euler*a; par=make_params('gum',[u,a])
      case('gev')
         if(size(lmom)<3) return
         t3=lmom(3)/l2
         c=2.0_dp/(3.0_dp+t3)-log(2.0_dp)/log(3.0_dp)
         k=7.8590_dp*c+2.9554_dp*c*c
         if(abs(k)<1.0e-7_dp) then
            a=l2/log(2.0_dp); u=l1-euler*a; k=0.0_dp
         else
            gam1=gamma(1.0_dp+k)
            a=l2*k/((1.0_dp-2.0_dp**(-k))*gam1)
            u=l1-a*(1.0_dp-gam1)/k
         end if
         par=make_params('gev',[u,a,k])
      case default
         par%family=''; par%npar=0
      end select
   end function fit_lmoments

end module lmomco_moments
