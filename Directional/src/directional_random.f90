module directional_random
   use directional_kinds, only : dp, pi
   use directional_geometry, only : rotation_matrix, normalize_rows
   implicit none
   private
   public :: rvmf, rvonmises, rspcauchy, rpkbd
contains
   function randn() result(z)
      real(dp)::z,u1,u2
      call random_number(u1)
      call random_number(u2)
      u1=max(u1,tiny(1.0_dp))
      z=sqrt(-2*log(u1))*cos(2*pi*u2)
   end function
   function rvmf(n,mu,kappa) result(x)
      integer,intent(in)::n
      real(dp),intent(in)::mu(:),kappa
      real(dp)::x(n,size(mu)),v(n,size(mu)-1),w(n),b,x0,c,z,u,ta,nrm
      real(dp)::ini(size(mu)),rot(size(mu),size(mu))
      integer::i,j,p
      p=size(mu)
      if(kappa<=0)then
      do i=1,n
      do j=1,p
      x(i,j)=randn()
      end do
      nrm=sqrt(sum(x(i,:)**2))
      x(i,:)=x(i,:)/nrm
      end do
      return
      end if
      do i=1,n;do j=1,p-1;v(i,j)=randn();end do;nrm=sqrt(sum(v(i,:)**2));v(i,:)=v(i,:)/nrm;end do
      b=(-2*kappa+sqrt(4*kappa*kappa+real(p-1,dp)**2))/real(p-1,dp);x0=(1-b)/(1+b);c=kappa*x0+real(p-1,dp)*log(1-x0*x0)
      do i=1,n
      do
      z=beta_sym(.5_dp*real(p-1,dp))
      call random_number(u)
      w(i)=(1-(1+b)*z)/(1-(1-b)*z)
      ta=kappa*w(i)+real(p-1,dp)*log(1-x0*w(i))
      if(ta-c>=log(max(u,tiny(1.0_dp))))exit
      end do
      end do
      ini=0
      ini(p)=1
      rot=rotation_matrix(ini,mu/sqrt(sum(mu*mu)))
      do i=1,n
      x(i,1:p-1)=v(i,:)*sqrt(max(0.0_dp,1-w(i)*w(i)))
      x(i,p)=w(i)
      x(i,:)=matmul(rot,x(i,:))
      end do
   end function
   function beta_sym(a) result(x)
      real(dp),intent(in)::a;real(dp)::x,g1,g2;g1=gamma_rand(a);g2=gamma_rand(a);x=g1/(g1+g2)
   end function
   recursive function gamma_rand(a) result(g)
      real(dp),intent(in)::a;real(dp)::g,d,c,x,v,u
      if(a<1)then
      call random_number(u)
      g=gamma_rand(a+1)*u**(1/a)
      return
      end if
      d=a-1.0_dp/3
      c=1/sqrt(9*d)
      do
      x=randn()
      v=(1+c*x)**3
      if(v<=0)cycle
      call random_number(u)
      if(u<1-.0331_dp*x**4 .or. log(u)<.5_dp*x*x+d*(1-v+log(v)))exit
      end do
      g=d*v
   end function
   function rvonmises(n,m,kappa,rads) result(u)
      integer,intent(in)::n
      real(dp),intent(in)::m,kappa
      logical,intent(in),optional::rads
      real(dp)::u(n),mu(2),x(n,2)
      logical::rr
      integer::i
      rr=.true.
      if(present(rads))rr=rads
      mu=[cos(m),sin(m)]
      x=rvmf(n,mu,kappa)
      do i=1,n
      u(i)=modulo(atan2(x(i,2),x(i,1)),2*pi)
      end do
      if(.not.rr)u=u*180/pi
   end function
   function rspcauchy(n,mu,rho) result(x)
      integer,intent(in)::n
      real(dp),intent(in)::mu(:),rho
      real(dp)::x(n,size(mu)),z(size(mu)),v(size(mu)),nr
      integer::i,j
      do i=1,n
      do j=1,size(mu)
      z(j)=randn()
      end do
      z=z/sqrt(sum(z*z))
      v=(1-rho*rho)*z+2*rho*mu
      nr=sqrt(sum(v*v))
      x(i,:)=v/nr
      end do
   end function
   function rpkbd(n,mu,rho) result(x)
      integer,intent(in)::n;real(dp),intent(in)::mu(:),rho;real(dp)::x(n,size(mu))
      real(dp)::muv(size(mu)),z(size(mu)),mz,zz,com,qa,u,lam,bstar,b1,b2,lo,hi,c,d,fc,fd,gr,lr
      integer::i,j,p
      p=size(mu);muv=mu/max(sqrt(sum(mu*mu)),tiny(1.0_dp));lam=2.0_dp*rho/(1.0_dp+rho*rho)
      if(rho<=1e-12_dp)then
         do i=1,n;do j=1,p;x(i,j)=randn();end do;x(i,:)=x(i,:)/sqrt(sum(x(i,:)**2));end do;return
      end if
      lo=max(lam*(2.0_dp-lam)+1e-10_dp,1e-10_dp);hi=1.0_dp-1e-10_dp;gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
      c=hi-gr*(hi-lo);d=lo+gr*(hi-lo);fc=pkbd_env_obj(c,lam,p);fd=pkbd_env_obj(d,lam,p)
      do i=1,120
         if(abs(hi-lo)<1e-10_dp)exit
         if(fc<fd)then;hi=d;d=c;fd=fc;c=hi-gr*(hi-lo);fc=pkbd_env_obj(c,lam,p)
         else;lo=c;c=d;fc=fd;d=lo+gr*(hi-lo);fd=pkbd_env_obj(d,lam,p);end if
      end do
      bstar=0.5_dp*(lo+hi);b1=bstar/(1.0_dp-bstar);b2=-1.0_dp+1.0_dp/sqrt(1.0_dp-bstar)
      i=0
      do while(i<n)
         do j=1,p;z(j)=randn();end do;mz=dot_product(z,muv);zz=sum(z*z);com=sqrt(zz+b1*mz*mz);qa=mz*(1.0_dp+b2)/com
         lr=0.5_dp*real(p,dp)*(-log(max(tiny(1.0_dp),1.0_dp-lam*qa))+log(max(tiny(1.0_dp),1.0_dp-bstar*qa*qa)) &
            -log(2.0_dp/(1.0_dp+sqrt(max(0.0_dp,1.0_dp-lam*lam/bstar)))))
         call random_number(u)
         if(log(max(u,tiny(1.0_dp)))<=lr)then;i=i+1;x(i,:)=(z+b2*mz*muv)/com;end if
      end do
   end function
   pure real(dp) function pkbd_env_obj(be,lam,p) result(v)
      real(dp),intent(in)::be,lam;integer,intent(in)::p
      v=real(p,dp)*log((1.0_dp+sqrt(max(0.0_dp,1.0_dp-lam*lam)))/(1.0_dp+sqrt(max(0.0_dp,1.0_dp-lam*lam/be))))-log(max(tiny(1.0_dp),1.0_dp-be))
   end function

end module directional_random
