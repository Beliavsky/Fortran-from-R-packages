module directional_advanced
   use directional_kinds, only : dp, pi
   use directional_inference, only : vmf_mle, vmf_mle_result
   use directional_linalg, only : solve_linear
   implicit none
   private
   type, public :: vector_mle_result
      real(dp), allocatable :: mu(:)
      real(dp) :: loglik=-huge(1.0_dp)
      real(dp) :: gamma=0.0_dp
      real(dp) :: angle=0.0_dp
   end type
   public :: iag_mle, sipc_mle, cipc_mle, gcpc_mle
contains
   function iag_mle(x,tol,maxit) result(res)
      real(dp),intent(in)::x(:,:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxit
      type(vector_mle_result)::res
      real(dp)::m(size(x,2)),cand(size(x,2)),best,val,step,eps,nr
      integer::j,it,nit,d,n
      n=size(x,1);d=size(x,2);eps=1e-7_dp;if(present(tol))eps=tol;nit=1200;if(present(maxit))nit=maxit
      allocate(res%mu(d));m=sum(x,dim=1)/real(n,dp);nr=sqrt(sum(m*m));if(nr<1e-8_dp)m=1.0_dp/sqrt(real(d,dp))
      m=m/max(sqrt(sum(m*m)),tiny(1.0_dp));m=2.0_dp*m
      best=iag_loglik(x,m);step=0.5_dp
      do it=1,nit
         val=best
         do j=1,d
            cand=m;cand(j)=cand(j)+step
            if(iag_loglik(x,cand)>best)then;m=cand;best=iag_loglik(x,cand);cycle;end if
            cand=m;cand(j)=cand(j)-step
            if(iag_loglik(x,cand)>best)then;m=cand;best=iag_loglik(x,cand);end if
         end do
         if(best<=val+eps*max(1.0_dp,abs(best)))then
            step=0.5_dp*step
            if(step<eps)exit
         end if
      end do
      res%mu=m;res%gamma=sqrt(sum(m*m));res%loglik=best
   end function

   pure real(dp) function iag_loglik(x,mu) result(ll)
      real(dp),intent(in)::x(:,:),mu(:)
      integer::n,d,p,i,j
      real(dp)::a,g2,phid,ncdf,m0,m1,m2,mp
      n=size(x,1);d=size(x,2);p=d-1;g2=sum(mu*mu);ll=-0.5_dp*real(n*p,dp)*log(2.0_dp*pi)
      do i=1,n
         a=dot_product(x(i,:),mu);phid=exp(-0.5_dp*a*a)/sqrt(2.0_dp*pi);ncdf=0.5_dp*erfc(-a/sqrt(2.0_dp))
         m0=ncdf
         if(p==1)then;mp=a*m0+phid
         else
            m1=a*m0+phid
            if(p==2)then;mp=(1.0_dp+a*a)*ncdf+a*phid
            else
               m2=(1.0_dp+a*a)*ncdf+a*phid
               do j=3,p
                  mp=a*m2+real(j-1,dp)*m1;m1=m2;m2=mp
               end do
               mp=m2
            end if
         end if
         ll=ll+0.5_dp*(a*a-g2)+log(max(mp,tiny(1.0_dp)))
      end do
   end function

   function sipc_mle(x,tol,maxit) result(res)
      real(dp),intent(in)::x(:,:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxit
      type(vector_mle_result)::res
      type(vector_mle_result)::ini
      real(dp)::m(size(x,2)),cand(size(x,2)),best,val,step,eps
      integer::j,it,nit
      eps=1e-6_dp;if(present(tol))eps=tol;nit=1500;if(present(maxit))nit=maxit
      ini=iag_mle(x);m=ini%mu;best=sipc_loglik(x,m);step=0.5_dp
      do it=1,nit
         val=best
         do j=1,size(m)
            cand=m;cand(j)=cand(j)+step
            if(sipc_loglik(x,cand)>best)then;m=cand;best=sipc_loglik(x,cand);cycle;end if
            cand=m;cand(j)=cand(j)-step
            if(sipc_loglik(x,cand)>best)then;m=cand;best=sipc_loglik(x,cand);end if
         end do
         if(best<=val+eps*max(1.0_dp,abs(best)))then;step=0.5_dp*step;if(step<eps)exit;end if
      end do
      allocate(res%mu(size(m)));res%mu=m;res%gamma=sqrt(sum(m*m));res%loglik=best-real(size(x,1),dp)*log(4.0_dp*pi*pi)
   end function

   pure real(dp) function sipc_loglik(x,m) result(ll)
      real(dp),intent(in)::x(:,:),m(:)
      real(dp)::a,d,sqd,up;integer::i
      ll=0.0_dp
      do i=1,size(x,1)
         a=dot_product(x(i,:),m);d=sum(m*m)+1.0_dp-a*a;sqd=sqrt(max(d,tiny(1.0_dp)))
         up=(sum(m*m)+1.0_dp)*sqd*(atan2(sqd,-a)-atan2(sqd,a)+pi)+2.0_dp*a*d
         ll=ll+log(max(up,tiny(1.0_dp)))-2.0_dp*log(max(d,tiny(1.0_dp)))
      end do
   end function

   function cipc_mle(u,rads,tol,maxit) result(res)
      real(dp),intent(in)::u(:);logical,intent(in),optional::rads;real(dp),intent(in),optional::tol;integer,intent(in),optional::maxit
      type(vector_mle_result)::res
      real(dp)::z(size(u)),x(size(u),2),m(2),a(size(u)),c,ci(size(u)),up(size(u),2),g(2),h(2,2),delta(2),old,ll,eps
      integer::i,j,k,it,nit,info,n
      logical::rr
      rr=.false.;if(present(rads))rr=rads;z=u;if(.not.rr)z=z*pi/180.0_dp
      x(:,1)=cos(z);x(:,2)=sin(z);n=size(u);m=sum(x,dim=1)/real(n,dp);eps=1e-6_dp;if(present(tol))eps=tol;nit=300;if(present(maxit))nit=maxit
      old=-huge(1.0_dp)
      do it=1,nit
         c=sqrt(sum(m*m)+1.0_dp);a=matmul(x,m);ci=1.0_dp/max(c-a,tiny(1.0_dp))
         do i=1,n;up(i,:)=x(i,:)-m/c;end do
         g=sum(up*spread(ci,2,2),dim=1);h=0.0_dp
         do j=1,2;do k=1,2
            h(j,k)=sum(up(:,j)*up(:,k)*ci*ci)-sum(ci)*(merge(c,0.0_dp,j==k)-m(j)*m(k)/c)/(c*c)
         end do;end do
         call solve_linear(h,g,delta,info);if(info/=0)exit;m=m-delta
         c=sqrt(sum(m*m)+1.0_dp);ll=-sum(log(max(c-matmul(x,m),tiny(1.0_dp))))
         if(abs(ll-old)<eps)exit;old=ll
      end do
      allocate(res%mu(2));res%mu=m;res%gamma=sqrt(sum(m*m));res%angle=modulo(atan2(m(2),m(1)),2*pi);if(.not.rr)res%angle=res%angle*180.0_dp/pi
      res%loglik=ll-real(n,dp)*log(2.0_dp*pi)
   end function

   function gcpc_mle(u,rads,tol,maxit) result(res)
      real(dp),intent(in)::u(:);logical,intent(in),optional::rads
      real(dp),intent(in),optional::tol;integer,intent(in),optional::maxit
      type(vector_mle_result)::res
      type(vector_mle_result)::ini
      real(dp)::z(size(u)),x(size(u),2),par(3),cand(3),step(3),best,val,eps,rho
      integer::i,j,it,nit
      logical::rr
      rr=.false.;if(present(rads))rr=rads;z=u;if(.not.rr)z=z*pi/180.0_dp
      x(:,1)=cos(z);x(:,2)=sin(z);ini=cipc_mle(u,rads=rr)
      par(1:2)=ini%mu;par(3)=0.0_dp;step=[0.4_dp,0.4_dp,0.4_dp]
      eps=1e-6_dp;if(present(tol))eps=tol;nit=1000;if(present(maxit))nit=maxit
      best=gcpc_loglik(x,par)
      do it=1,nit
         val=best
         do j=1,3
            cand=par;cand(j)=cand(j)+step(j)
            if(gcpc_loglik(x,cand)>best)then;par=cand;best=gcpc_loglik(x,cand);cycle;end if
            cand=par;cand(j)=cand(j)-step(j)
            if(gcpc_loglik(x,cand)>best)then;par=cand;best=gcpc_loglik(x,cand);end if
         end do
         if(best<=val+eps*max(1.0_dp,abs(best)))then
            step=0.6_dp*step
            if(maxval(step)<eps)exit
         end if
      end do
      allocate(res%mu(3));rho=exp(par(3));res%mu=[par(1),par(2),rho]
      res%gamma=sqrt(sum(par(1:2)**2));res%angle=modulo(atan2(par(2),par(1)),2*pi)
      if(.not.rr)res%angle=res%angle*180.0_dp/pi
      res%loglik=best-real(size(u),dp)*log(2.0_dp*pi)
   end function

   pure real(dp) function gcpc_loglik(x,par) result(ll)
      real(dp),intent(in)::x(:,:),par(3)
      real(dp)::mu(2),rho,g2,nr,ksi(2),sinv(2,2),a,b,den
      integer::i,n
      mu=par(1:2);rho=exp(max(-12.0_dp,min(12.0_dp,par(3))));g2=sum(mu*mu);nr=sqrt(max(g2,tiny(1.0_dp)))
      ksi=mu/nr
      sinv(1,1)=ksi(1)**2+ksi(2)**2/rho
      sinv(2,2)=ksi(2)**2+ksi(1)**2/rho
      sinv(1,2)=ksi(1)*ksi(2)*(1.0_dp-1.0_dp/rho);sinv(2,1)=sinv(1,2)
      n=size(x,1);ll=-0.5_dp*real(n,dp)*log(rho)
      do i=1,n
         a=dot_product(x(i,:),mu);b=dot_product(x(i,:),matmul(sinv,x(i,:)))
         den=b*sqrt(g2+1.0_dp)-a*sqrt(max(b,tiny(1.0_dp)))
         ll=ll-log(max(den,tiny(1.0_dp)))
      end do
   end function
end module directional_advanced
