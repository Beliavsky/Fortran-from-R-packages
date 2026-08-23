module directional_fisher_bingham
   use directional_kinds, only : dp, pi
   use directional_linalg, only : symmetric_eigen, sort_eigen_desc, svd3, det3
   use directional_geometry, only : rotation_matrix
   implicit none
   private
   public :: fb_saddle, dkent, rbingham, rkent, kent_mle, matrixfisher_mle, rmatrixfisher
contains
   function fb_saddle(gam,lam_in) result(logcon)
      real(dp),intent(in)::gam(:),lam_in(:);real(dp)::logcon(3),lam(size(lam_in)),mina,shift,low,up,tau,mid,r3,r4,ta,c1
      integer::p,it
      p=size(gam);lam=lam_in;call sort_ascending(lam);mina=minval(lam);shift=0;if(mina<=0)then;shift=abs(mina)+1;lam=lam+shift;end if
      low=lam(1)-0.25_dp*p-0.5_dp*sqrt(0.25_dp*p*p+p*maxval(gam)**2);up=lam(1)-0.25_dp-0.5_dp*sqrt(0.25_dp+minval(gam)**2)
      do it=1,200;mid=0.5_dp*(low+up);if(saddle_eq(mid,gam,lam)>0)then;up=mid;else;low=mid;end if;if(abs(up-low)<1e-10_dp)exit;end do;tau=0.5_dp*(low+up)
      r3=kfb(3,gam,lam,tau)/kfb(2,gam,lam,tau)**1.5_dp;r4=kfb(4,gam,lam,tau)/kfb(2,gam,lam,tau)**2;ta=r4/8-5*r3*r3/24
      c1=0.5_dp*log(2.0_dp)+0.5_dp*(p-1)*log(pi)-0.5_dp*log(kfb(2,gam,lam,tau))-0.5_dp*sum(log(lam-tau))-tau+0.25_dp*sum(gam*gam/(lam-tau))+shift
      logcon=[c1,c1+log(max(tiny(1.0_dp),1+ta)),c1+ta]
   end function
   pure real(dp) function saddle_eq(t,g,l) result(v)
      real(dp),intent(in)::t,g(:),l(:);v=sum(0.5_dp/(l-t)+0.25_dp*g*g/(l-t)**2)-1
   end function
   pure real(dp) function kfb(j,g,l,t) result(v)
      integer,intent(in)::j;real(dp),intent(in)::g(:),l(:),t;integer::f1,f2,k
      if(j==1)then;v=sum(0.5_dp/(l-t)+0.25_dp*g*g/(l-t)**2);return;end if
      f1=1;do k=2,j-1;f1=f1*k;end do;f2=f1*j;v=sum(0.5_dp*f1/(l-t)**j+0.25_dp*f2*g*g/(l-t)**(j+1))
   end function
   function dkent(y,g,kappa,beta,logden) result(v)
      real(dp),intent(in)::y(:,:),g(3,3),kappa,beta;logical,intent(in),optional::logden;real(dp)::v(size(y,1)),gam(3),lam(3),lc(3),l;integer::i;logical::ll
      gam=[0.0_dp,kappa,0.0_dp];lam=[0.0_dp,-beta,beta];lc=fb_saddle(gam,lam);ll=.false.;if(present(logden))ll=logden
      do i=1,size(y,1);l=-lc(3)+kappa*dot_product(y(i,:),g(:,1))+beta*dot_product(y(i,:),g(:,2))**2-beta*dot_product(y(i,:),g(:,3))**2;v(i)=merge(l,exp(l),ll);end do
   end function

   function rbingham(n,a) result(x)
      integer,intent(in)::n;real(dp),intent(in)::a(:,:);real(dp)::x(n,size(a,1)),eval(size(a,1)),vec(size(a,1),size(a,1)),lam(size(a,1)),sig(size(a,1)),z(size(a,1)),y(size(a,1)),u,lratio,nr;integer::i,j,p
      p=size(a,1);call symmetric_eigen(a,eval,vec);call sort_eigen_desc(eval,vec);lam=eval-eval(p);sig=sqrt(1/(1+2*max(lam,0.0_dp)))
      do i=1,n
         do
            do j=1,p;z(j)=randn()*sig(j);end do;nr=sqrt(sum(z*z));y=z/nr;lratio=-sum(y*y*lam)-0.5_dp*p*log(real(p,dp))+0.5_dp*(p-1)+0.5_dp*p*log(sum(y*y*(1+2*lam)));call random_number(u);if(log(max(u,tiny(1.0_dp)))<lratio)exit
         end do
         x(i,:)=matmul(vec,y)
      end do
   end function

   function rkent(n,kappa,m,beta) result(x)
      integer,intent(in)::n;real(dp),intent(in)::kappa,m(3),beta;real(dp)::x(n,3),g(3,3),a(3,3),cand(1,3),u,lcur,lmax,mm(3),e1(3),e2(3);integer::i
      mm=m/sqrt(sum(m*m));e1=[1.0_dp,0.0_dp,0.0_dp];if(abs(dot_product(e1,mm))>0.9_dp)e1=[0.0_dp,1.0_dp,0.0_dp];e1=e1-dot_product(e1,mm)*mm;e1=e1/sqrt(sum(e1*e1));e2=cross3(mm,e1);g(:,1)=mm;g(:,2)=e1;g(:,3)=e2
      a=0;a=reshape([ -beta,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,beta],[3,3]);lmax=kappa+abs(beta)
      do i=1,n
         do;cand=uniform_sphere(1,3);lcur=kappa*dot_product(cand(1,:),mm)+beta*dot_product(cand(1,:),e1)**2-beta*dot_product(cand(1,:),e2)**2;call random_number(u);if(log(max(u,tiny(1.0_dp)))<=lcur-lmax)exit;end do;x(i,:)=cand(1,:);end do
   end function

   subroutine kent_mle(x,g,kappa,beta,psi,logcon,loglik)
      real(dp),intent(in)::x(:,:);real(dp),intent(out)::g(3,3),kappa,beta,psi,logcon,loglik
      real(dp)::xb(3),s(3,3),h(3,3),bmat(3,3),km(3,3),xg(size(x,1),3),rbar
      real(dp)::theta,phi,ct,st,cp,sp,xg1,xg2,xg3,best,val,sk,sb,nk,nb
      integer::i,j,it,n
      n=size(x,1);xb=sum(x,dim=1)/real(n,dp);rbar=sqrt(sum(xb*xb));xb=xb/max(rbar,tiny(1.0_dp))
      s=0.0_dp
      do i=1,n
         do j=1,3;s(:,j)=s(:,j)+x(i,:)*x(i,j);end do
      end do
      s=s/real(n,dp)
      theta=acos(max(-1.0_dp,min(1.0_dp,xb(1))));phi=modulo(atan2(xb(3),xb(2)),2.0_dp*pi)
      ct=cos(theta);st=sin(theta);cp=cos(phi);sp=sin(phi)
      h=reshape([ct,st*cp,st*sp,-st,ct*cp,ct*sp,0.0_dp,-sp,cp],[3,3])
      bmat=matmul(transpose(h),matmul(s,h))
      psi=0.5_dp*atan2(2.0_dp*bmat(2,3),bmat(2,2)-bmat(3,3))
      km=reshape([1.0_dp,0.0_dp,0.0_dp,0.0_dp,cos(psi),sin(psi),0.0_dp,-sin(psi),cos(psi)],[3,3])
      g=matmul(h,km);xg=matmul(x,g);xg1=sum(xg(:,1));xg2=sum(xg(:,2)**2);xg3=sum(xg(:,3)**2)
      kappa=max(1e-4_dp,rbar*(3.0_dp-rbar*rbar)/max(1e-6_dp,1.0_dp-rbar*rbar));beta=min(0.45_dp*kappa,max(0.0_dp,0.25_dp*kappa))
      best=kent_nll(n,xg1,xg2,xg3,kappa,beta);sk=max(0.25_dp*kappa,0.5_dp);sb=max(0.15_dp*kappa,0.1_dp)
      do it=1,160
         nk=max(1e-6_dp,kappa+sk);nb=min(0.499_dp*nk,max(0.0_dp,beta));val=kent_nll(n,xg1,xg2,xg3,nk,nb);if(val<best)then;kappa=nk;beta=nb;best=val;cycle;end if
         nk=max(1e-6_dp,kappa-sk);nb=min(0.499_dp*nk,max(0.0_dp,beta));val=kent_nll(n,xg1,xg2,xg3,nk,nb);if(val<best)then;kappa=nk;beta=nb;best=val;cycle;end if
         nb=min(0.499_dp*kappa,beta+sb);val=kent_nll(n,xg1,xg2,xg3,kappa,nb);if(val<best)then;beta=nb;best=val;cycle;end if
         nb=max(0.0_dp,beta-sb);val=kent_nll(n,xg1,xg2,xg3,kappa,nb);if(val<best)then;beta=nb;best=val;cycle;end if
         sk=0.7_dp*sk;sb=0.7_dp*sb;if(max(sk,sb)<1e-7_dp)exit
      end do
      block
         real(dp) :: lc3(3)
         lc3=fb_saddle([0.0_dp,kappa,0.0_dp],[0.0_dp,-beta,beta]);logcon=lc3(3)
      end block
      loglik=-best
   end subroutine

   real(dp) function kent_nll(n,xg1,xg2,xg3,k,b) result(v)
      integer,intent(in)::n;real(dp),intent(in)::xg1,xg2,xg3,k,b;real(dp)::lc(3)
      if(k<=0.0_dp .or. b<0.0_dp .or. 2.0_dp*b>=k)then;v=huge(1.0_dp);return;end if
      lc=fb_saddle([0.0_dp,k,0.0_dp],[0.0_dp,-b,b]);v=real(n,dp)*lc(3)-k*xg1-b*(xg2-xg3)
   end function

   subroutine matrixfisher_mle(x,u,s,v)
      real(dp),intent(in)::x(:,:,:);real(dp),intent(out)::u(3,3),s(3),v(3,3);real(dp)::xb(3,3);integer::i
      xb=0;do i=1,size(x,3);xb=xb+x(:,:,i);end do;xb=xb/size(x,3);call svd3(xb,u,s,v)
   end subroutine

   function rmatrixfisher(n,f) result(x)
      integer,intent(in)::n;real(dp),intent(in)::f(3,3);real(dp)::x(3,3,n),u(3,3),v(3,3),s(3),a4(4,4),q(n,4),qq(4),r(3,3);integer::i
      call svd3(f,u,s,v);if(det3(matmul(u,transpose(v)))<0)then;u(:,3)=-u(:,3);s(3)=-s(3);end if
      a4=0;a4(2,2)=2*(s(2)+s(3));a4(3,3)=2*(s(1)+s(3));a4(4,4)=2*(s(1)+s(2));q=rbingham(n,a4)
      do i=1,n;qq=q(i,:);r=quat_to_rot(qq);x(:,:,i)=matmul(u,matmul(r,transpose(v)));end do
   end function
   pure function quat_to_rot(q) result(r)
      real(dp),intent(in)::q(4);real(dp)::r(3,3);real(dp)::x1,x2,x3,x4;x1=q(1);x2=q(2);x3=q(3);x4=q(4)
      r(1,1)=x1*x1+x2*x2-x3*x3-x4*x4;r(2,1)=2*(x1*x4+x3*x2);r(3,1)=-2*(x1*x3-x2*x4)
      r(1,2)=-2*(x1*x4-x3*x2);r(2,2)=x1*x1+x3*x3-x2*x2-x4*x4;r(3,2)=2*(x1*x2+x3*x4)
      r(1,3)=2*(x1*x3+x2*x4);r(2,3)=-2*(x1*x2-x3*x4);r(3,3)=x1*x1+x4*x4-x2*x2-x3*x3
   end function
   function uniform_sphere(n,p) result(x)
      integer,intent(in)::n,p;real(dp)::x(n,p),nr;integer::i,j;do i=1,n;do j=1,p;x(i,j)=randn();end do;nr=sqrt(sum(x(i,:)**2));x(i,:)=x(i,:)/nr;end do
   end function
   real(dp) function randn() result(z)
      real(dp)::u1,u2;call random_number(u1);call random_number(u2);z=sqrt(-2*log(max(u1,tiny(1.0_dp))))*cos(2*pi*u2)
   end function
   pure function cross3(a,b) result(c)
      real(dp),intent(in)::a(3),b(3);real(dp)::c(3);c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
   end function
   subroutine sort_ascending(x)
      real(dp),intent(inout)::x(:);integer::i,j;real(dp)::t;do i=2,size(x);t=x(i);j=i-1;do while(j>=1)
         if(x(j)<=t) exit
         x(j+1)=x(j);j=j-1
      end do
      x(j+1)=t
   end do
   end subroutine
end module directional_fisher_bingham
