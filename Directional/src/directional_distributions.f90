module directional_distributions
   use directional_kinds, only : dp, pi
   use directional_special, only : log_i0, log_bessel_i, normal_pdf, normal_cdf, log_beta
   implicit none
   private
   public :: dvm, dmmvm, dvmf, dcardio, dwrapcauchy, dwrapnormal, dcircbeta, dcircexp
   public :: dcircpurka, dspcauchy, dpkbd, dpurka, iagd, desag3, pwrapcauchy, pvm_numeric
contains
   elemental real(dp) function dvm(x,m,k,rads,logden) result(v)
      real(dp),intent(in)::x,m,k; logical,intent(in),optional::rads,logden
      real(dp)::z,l; logical::rr,ll
      rr=.false.;if(present(rads))rr=rads;ll=.false.;if(present(logden))ll=logden
      z=x;if(.not.rr)z=z*pi/180.0_dp;l=k*cos(z-m)-log(2*pi)-log_i0(k);v=merge(l,exp(l),ll)
   end function dvm
   elemental real(dp) function dmmvm(x,m,k,nmode,rads,logden) result(v)
      real(dp),intent(in)::x,m,k;integer,intent(in)::nmode;logical,intent(in),optional::rads,logden
      real(dp)::z,l;logical::rr,ll;rr=.false.;if(present(rads))rr=rads;ll=.false.;if(present(logden))ll=logden
      z=x;if(.not.rr)z=z*pi/180;l=k*cos(real(nmode,dp)*(z-m))-log(2*pi)-log_i0(k);v=merge(l,exp(l),ll)
   end function dmmvm
   function dvmf(y,mu,k,logden) result(v)
      real(dp),intent(in)::y(:,:),mu(:),k;logical,intent(in),optional::logden;real(dp)::v(size(y,1)),l,nu
      integer::i,p;logical::ll;p=size(y,2);nu=0.5_dp*p-1;l=nu*log(k)-0.5_dp*p*log(2*pi)-log_bessel_i(nu,k)
      ll=.false.
      if(present(logden))ll=logden
      do i=1,size(y,1)
      v(i)=l+k*dot_product(mu,y(i,:))
      if(.not.ll)v(i)=exp(v(i))
      end do
   end function dvmf
   elemental real(dp) function dcardio(x,m,rho,rads,logden) result(v)
      real(dp),intent(in)::x,m,rho;logical,intent(in),optional::rads,logden;real(dp)::z,l;logical::rr,ll
      rr=.false.;if(present(rads))rr=rads;ll=.false.;if(present(logden))ll=logden;z=x;if(.not.rr)z=z*pi/180
      l=-log(2*pi)+log(1+2*rho*cos(z-m));v=merge(l,exp(l),ll)
   end function dcardio
   elemental real(dp) function dwrapcauchy(x,m,rho,rads,logden) result(v)
      real(dp),intent(in)::x,m,rho;logical,intent(in),optional::rads,logden;real(dp)::z,l;logical::rr,ll
      rr=.false.;if(present(rads))rr=rads;ll=.false.;if(present(logden))ll=logden;z=x;if(.not.rr)z=z*pi/180
      l=-log(2*pi)+log(1-rho*rho)-log(1+rho*rho-2*rho*cos(z-m));v=merge(l,exp(l),ll)
   end function dwrapcauchy
   elemental real(dp) function dwrapnormal(x,m,rho,rads,logden) result(v)
      real(dp),intent(in)::x,m,rho;logical,intent(in),optional::rads,logden;real(dp)::z,s,l;logical::rr,ll;integer::j
      rr=.false.;if(present(rads))rr=rads;ll=.false.;if(present(logden))ll=logden;z=x;if(.not.rr)z=z*pi/180;s=0
      do j=1,100;s=s+rho**(j*j)*cos(real(j,dp)*(z-m));end do;l=-log(2*pi)+log(1+2*s);v=merge(l,exp(l),ll)
   end function dwrapnormal
   elemental real(dp) function dcircbeta(x,m,a,b,rads,logden) result(v)
      real(dp),intent(in)::x,m,a,b;logical,intent(in),optional::rads,logden;real(dp)::z,c,l;logical::rr,ll
      rr=.false.;if(present(rads))rr=rads;ll=.false.;if(present(logden))ll=logden;z=x;if(.not.rr)z=z*pi/180;c=cos(z-m)
      l=-(a+b)*log(2.0_dp)-log_beta(a,b)+(a-.5_dp)*log(1+c)+(b-.5_dp)*log(1-c);v=merge(l,exp(l),ll)
   end function dcircbeta
   elemental real(dp) function dcircexp(x,lambda,rads,logden) result(v)
      real(dp),intent(in)::x,lambda;logical,intent(in),optional::rads,logden;real(dp)::z,l;logical::rr,ll
      rr=.false.;if(present(rads))rr=rads;ll=.false.;if(present(logden))ll=logden;z=x;if(.not.rr)z=z*pi/180
      l=log(lambda)-lambda*z-log(1-exp(-2*pi*lambda));v=merge(l,exp(l),ll)
   end function dcircexp
   elemental real(dp) function dcircpurka(x,m,a,rads,logden) result(v)
      real(dp),intent(in)::x,m,a;logical,intent(in),optional::rads,logden;real(dp)::z,l,cs;logical::rr,ll
      rr=.false.;if(present(rads))rr=rads;ll=.false.;if(present(logden))ll=logden;z=x;if(.not.rr)z=z*pi/180
      cs=max(-1.0_dp,min(1.0_dp,cos(z-m)));l=log(a)-log(2.0_dp)-log(1-exp(-a*pi))-a*acos(cs);v=merge(l,exp(l),ll)
   end function dcircpurka
   function dspcauchy(y,mu,rho,logden) result(v)
      real(dp),intent(in)::y(:,:),mu(:),rho
      logical,intent(in),optional::logden
      real(dp)::v(size(y,1)),a,l
      integer::i,d
      logical::ll
      d=size(y,2)-1;ll=.false.;if(present(logden))ll=logden;do i=1,size(y,1);a=dot_product(y(i,:),mu)
      l=d*log(1-rho*rho)-d*log(1+rho*rho-2*rho*a)+log_gamma(.5_dp*(d+1))-.5_dp*(d+1)*log(pi)-log(2.0_dp)
      v(i)=merge(l,exp(l),ll);end do
   end function dspcauchy
   function dpkbd(y,mu,rho,logden) result(v)
      real(dp),intent(in)::y(:,:),mu(:),rho
      logical,intent(in),optional::logden
      real(dp)::v(size(y,1)),a,l
      integer::i,d
      logical::ll
      d=size(y,2);ll=.false.;if(present(logden))ll=logden;do i=1,size(y,1);a=dot_product(y(i,:),mu)
      l=log(1-rho*rho)-.5_dp*d*log(1+rho*rho-2*rho*a)+log_gamma(.5_dp*d)-.5_dp*d*log(pi)-log(2.0_dp)
      v(i)=merge(l,exp(l),ll);end do
   end function dpkbd
   function dpurka(y,theta,a,logden) result(v)
      real(dp),intent(in)::y(:,:),theta(:),a
      logical,intent(in),optional::logden
      real(dp)::v(size(y,1)),ang,l
      integer::i,p
      logical::ll
      p=size(y,2)
      ll=.false.
      if(present(logden))ll=logden
      do i=1,size(y,1)
      ang=acos(max(-1.0_dp,min(1.0_dp,dot_product(y(i,:),theta))))
      if(p==3)then;l=log(a*a+1)-log(2*pi)-log(1+exp(-a*pi))-a*ang
      else;l=log_gamma(.5_dp*p)-.5_dp*p*log(pi)+log_bessel_i(real(p-1,dp),a)-a*ang;end if
      v(i)=merge(l,exp(l),ll);end do
   end function dpurka
   function iagd(y,mu,logden) result(v)
      real(dp),intent(in)::y(:,:),mu(:)
      logical,intent(in),optional::logden
      real(dp)::v(size(y,1)),aa,mp0,mp1,mp,l
      integer::i,j,p
      logical::ll
      p=size(y,2)-1;ll=.false.;if(present(logden))ll=logden
      do i=1,size(y,1);aa=dot_product(y(i,:),mu);mp0=1;mp1=aa*normal_cdf(aa)+normal_pdf(aa)
         if(p==1)then;mp=mp1;else;mp=(1+aa*aa)*normal_cdf(aa)+aa*normal_pdf(aa)
            do j=3,p; l=mp;mp=aa*mp+real(j-1,dp)*mp1;mp1=l;end do
         end if
         l=-.5_dp*p*log(2*pi)+.5_dp*(aa*aa-sum(mu*mu))+log(mp);v(i)=merge(l,exp(l),ll)
      end do
   end function iagd
   function desag3(y,mu,gam,logden) result(v)
      real(dp),intent(in)::y(:,:),mu(3),gam(2)
      logical,intent(in),optional::logden
      real(dp)::v(size(y,1)),m0,rl,x1(3),x2(3),va(3,3),g1,g2,a,m2,l,h
      integer::i,j,k
      logical::ll
      m0=sqrt(mu(2)**2+mu(3)**2)
      rl=sqrt(sum(mu*mu))
      x1=[-m0*m0,mu(1)*mu(2),mu(1)*mu(3)]/(m0*rl)
      x2=[0.0_dp,-mu(3),mu(2)]/m0
      h=sqrt(sum(gam*gam)+1)-1
      va=0
      do i=1,3
      va(i,i)=1
      do j=1,3
      va(i,j)=va(i,j)+gam(1)*(x1(i)*x1(j)-x2(i)*x2(j))+gam(2)*(x1(i)*x2(j)+x2(i)*x1(j))+h*(x1(i)*x1(j)+x2(i)*x2(j))
      end do
      end do
      ll=.false.
      if(present(logden))ll=logden
      do i=1,size(y,1)
      g2=dot_product(y(i,:),mu)
      g1=0
      do j=1,3
      do k=1,3
      g1=g1+y(i,j)*va(j,k)*y(i,k)
      end do
      end do
      a=g2/sqrt(g1)
      m2=(1+a*a)*normal_cdf(a)+a*normal_pdf(a)
      l=-log(2*pi)+.5_dp*a*a-.5_dp*sum(mu*mu)-1.5_dp*log(g1)+log(m2)
      v(i)=merge(l,exp(l),ll)
      end do
   end function desag3
   elemental real(dp) function pwrapcauchy(x,m,rho,rads) result(v)
      real(dp),intent(in)::x,m,rho;logical,intent(in),optional::rads;real(dp)::z,t;logical::rr
      rr=.false.
      if(present(rads))rr=rads
      z=x
      if(.not.rr)z=z*pi/180
      t=atan2((1+rho)*sin(.5_dp*(z-m)),(1-rho)*cos(.5_dp*(z-m)))
      v=modulo(t/pi+1.0_dp,2.0_dp)/2.0_dp
   end function pwrapcauchy
   real(dp) function pvm_numeric(x,m,k,rads) result(v)
      real(dp),intent(in)::x,m,k;logical,intent(in),optional::rads;real(dp)::z,h,s,t;logical::rr;integer::i,n
      rr=.false.
      if(present(rads))rr=rads
      z=x
      if(.not.rr)z=z*pi/180
      z=modulo(z,2*pi)
      n=800
      h=z/n
      s=.5_dp*(dvm(0.0_dp,m,k,.true.)+dvm(z,m,k,.true.))
      do i=1,n-1
      t=i*h
      s=s+dvm(t,m,k,.true.)
      end do
      v=max(0.0_dp,min(1.0_dp,h*s))
   end function pvm_numeric
end module directional_distributions
