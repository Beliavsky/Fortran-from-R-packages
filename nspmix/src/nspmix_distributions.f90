module nspmix_distributions
   use nspmix_kinds, only : dp
   use nspmix_types, only : disc_dist
   implicit none
   private
   public :: dnpnorm, pnpnorm, dnppois, pnppois, dnpgeom, pnpgeom, dnpnbinom, pnpnbinom
   public :: rnpnorm, rnppois, rnpgeom, rnpnbinom, dmix_disc
contains
   pure real(dp) function normal_pdf(x,mu,sd)
      real(dp),intent(in)::x,mu,sd
      normal_pdf=exp(-0.5_dp*((x-mu)/sd)**2)/(sd*sqrt(2.0_dp*acos(-1.0_dp)))
   end function
   pure real(dp) function normal_cdf(x,mu,sd)
      real(dp),intent(in)::x,mu,sd
      normal_cdf=0.5_dp*erfc(-(x-mu)/(sd*sqrt(2.0_dp)))
   end function

   pure real(dp) function dnpnorm(x,mix,sd)
      real(dp),intent(in)::x,sd; type(disc_dist),intent(in)::mix; integer::j
      dnpnorm=0.0_dp; do j=1,size(mix%pt); dnpnorm=dnpnorm+mix%pr(j)*normal_pdf(x,mix%pt(j),sd); end do
   end function
   pure real(dp) function pnpnorm(x,mix,sd)
      real(dp),intent(in)::x,sd; type(disc_dist),intent(in)::mix; integer::j
      pnpnorm=0.0_dp; do j=1,size(mix%pt); pnpnorm=pnpnorm+mix%pr(j)*normal_cdf(x,mix%pt(j),sd); end do
   end function
   pure real(dp) function dnppois(x,mix)
      integer,intent(in)::x; type(disc_dist),intent(in)::mix; integer::j; real(dp)::la
      dnppois=0.0_dp; if(x<0)return
      do j=1,size(mix%pt)
         la=max(mix%pt(j),1.0e-300_dp)
         dnppois=dnppois+mix%pr(j)*exp(real(x,dp)*log(la)-la-log_gamma(real(x+1,dp)))
      end do
   end function
   pure real(dp) function pnppois(x,mix)
      integer,intent(in)::x; type(disc_dist),intent(in)::mix; integer::j,k; real(dp)::la,t,s
      pnppois=0.0_dp; if(x<0)return
      do j=1,size(mix%pt)
         la=max(mix%pt(j),0.0_dp); t=exp(-la); s=t
         do k=1,x; t=t*la/real(k,dp); s=s+t; end do
         pnppois=pnppois+mix%pr(j)*min(s,1.0_dp)
      end do
   end function
   pure real(dp) function dnpgeom(x,mix)
      integer,intent(in)::x; type(disc_dist),intent(in)::mix; integer::j; real(dp)::p
      dnpgeom=0.0_dp; if(x<0)return
      do j=1,size(mix%pt); p=min(max(mix%pt(j),0.0_dp),1.0_dp); dnpgeom=dnpgeom+mix%pr(j)*p*(1.0_dp-p)**x; end do
   end function
   pure real(dp) function pnpgeom(x,mix)
      integer,intent(in)::x; type(disc_dist),intent(in)::mix; integer::j; real(dp)::p
      pnpgeom=0.0_dp; if(x<0)return
      do j=1,size(mix%pt); p=min(max(mix%pt(j),0.0_dp),1.0_dp); pnpgeom=pnpgeom+mix%pr(j)*(1.0_dp-(1.0_dp-p)**(x+1)); end do
   end function
   pure real(dp) function dnpnbinom(x,size_par,mix)
      integer,intent(in)::x; real(dp),intent(in)::size_par; type(disc_dist),intent(in)::mix; integer::j; real(dp)::p
      dnpnbinom=0.0_dp; if(x<0)return
      do j=1,size(mix%pt)
         p=min(max(mix%pt(j),1.0e-300_dp),1.0_dp-1.0e-15_dp)
         dnpnbinom=dnpnbinom+mix%pr(j)*exp(log_gamma(size_par+real(x,dp))-log_gamma(size_par)-log_gamma(real(x+1,dp)) &
              +size_par*log(p)+real(x,dp)*log(1.0_dp-p))
      end do
   end function
   pure real(dp) function pnpnbinom(x,size_par,mix)
      integer,intent(in)::x; real(dp),intent(in)::size_par; type(disc_dist),intent(in)::mix; integer::j,k; real(dp)::p,t,s
      pnpnbinom=0.0_dp; if(x<0)return
      do j=1,size(mix%pt)
         p=min(max(mix%pt(j),1.0e-300_dp),1.0_dp); t=p**size_par; s=t
         do k=1,x; t=t*(size_par+real(k-1,dp))/real(k,dp)*(1.0_dp-p); s=s+t; end do
         pnpnbinom=pnpnbinom+mix%pr(j)*min(s,1.0_dp)
      end do
   end function

   subroutine rnpnorm(n,mix,sd,x)
      integer,intent(in)::n; type(disc_dist),intent(in)::mix; real(dp),intent(in)::sd; real(dp),intent(out)::x(n)
      integer::i,j; real(dp)::u1,u2,r,z
      do i=1,n
         j=draw_component(mix); call random_number(u1); call random_number(u2); u1=max(u1,1.0e-300_dp)
         r=sqrt(-2.0_dp*log(u1)); z=r*cos(2.0_dp*acos(-1.0_dp)*u2); x(i)=mix%pt(j)+sd*z
      end do
   end subroutine
   subroutine rnppois(n,mix,x)
      integer,intent(in)::n; type(disc_dist),intent(in)::mix; integer,intent(out)::x(n)
      integer::i,j,k; real(dp)::u,t,s,la
      do i=1,n
         j=draw_component(mix); la=max(mix%pt(j),0.0_dp); call random_number(u); t=exp(-la); s=t; k=0
         do while(u>s .and. k<100000); k=k+1; t=t*la/real(k,dp); s=s+t; end do; x(i)=k
      end do
   end subroutine
   subroutine rnpgeom(n,mix,x)
      integer,intent(in)::n; type(disc_dist),intent(in)::mix; integer,intent(out)::x(n)
      integer::i,j; real(dp)::u,p
      do i=1,n; j=draw_component(mix); p=min(max(mix%pt(j),1.0e-15_dp),1.0_dp); call random_number(u)
         if(p>=1.0_dp) then; x(i)=0; else; x(i)=int(log(max(1.0_dp-u,1.0e-300_dp))/log(1.0_dp-p)); end if
      end do
   end subroutine
   subroutine rnpnbinom(n,size_par,mix,x)
      integer,intent(in)::n; real(dp),intent(in)::size_par; type(disc_dist),intent(in)::mix; integer,intent(out)::x(n)
      integer::i,j,k; real(dp)::u,p,t,s
      do i=1,n
         j=draw_component(mix); p=min(max(mix%pt(j),1.0e-12_dp),1.0_dp); call random_number(u); t=p**size_par; s=t; k=0
         do while(u>s .and. k<100000); k=k+1; t=t*(size_par+real(k-1,dp))/real(k,dp)*(1.0_dp-p); s=s+t; end do; x(i)=k
      end do
   end subroutine
   pure real(dp) function dmix_disc(x,mix)
      real(dp),intent(in)::x; type(disc_dist),intent(in)::mix; integer::j
      dmix_disc=0.0_dp
      do j=1,size(mix%pt)
         if(abs(x-mix%pt(j)) <= epsilon(1.0_dp)*max(1.0_dp,abs(x),abs(mix%pt(j)))) then
            dmix_disc=dmix_disc+mix%pr(j)
         end if
      end do
   end function

   integer function draw_component(mix)
      type(disc_dist),intent(in)::mix; real(dp)::u,s; integer::j
      call random_number(u); s=0.0_dp; draw_component=size(mix%pt)
      do j=1,size(mix%pt); s=s+mix%pr(j); if(u<=s) then; draw_component=j; return; end if; end do
   end function
end module nspmix_distributions
