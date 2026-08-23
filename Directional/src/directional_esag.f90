module directional_esag
   use directional_kinds, only : dp, pi
   use directional_special, only : normal_pdf, normal_cdf
   implicit none
   private
   public :: desag, resag, esag_parameters
contains
   subroutine esag_parameters(gam,lambda,rot,info)
      real(dp),intent(in)::gam(:)
      real(dp),allocatable,intent(out)::lambda(:),rot(:,:)
      integer,intent(out),optional::info
      integer::g,d,i,j,m,n,c1,c2,idx
      real(dp),allocatable::kk(:),theta(:),phi(:),grp(:),lo(:,:),la(:,:),rr(:,:),tmp(:,:)
      real(dp)::prodv,normg
      g=size(gam)
      d=nint(0.5_dp*(1.0_dp+sqrt(1.0_dp+8.0_dp*real(1+g,dp))))
      if(d<3 .or. d*(d-1)/2-1/=g)then
         allocate(lambda(0),rot(0,0));if(present(info))info=1;return
      end if
      allocate(lambda(d-1),rot(d-1,d-1),kk(d-2),theta(d-2))
      if(d>3)allocate(phi((d-2)*(d-3)/2));if(d==3)allocate(phi(0))
      kk=1.0_dp;theta=0.0_dp;if(size(phi)>0)phi=0.0_dp
      kk(1)=sqrt(gam(1)*gam(1)+gam(2)*gam(2))+1.0_dp
      theta(1)=atan2(gam(2),gam(1))
      idx=1
      if(d>=4)then
         do i=2,d-2
            c1=i*(i+1)/2
            c2=(i+1)*(i+2)/2-1
            n=c2-c1+1
            allocate(grp(n));grp=gam(c1:c2)
            normg=sqrt(sum(grp*grp));kk(i)=normg+1.0_dp
            do j=1,n-2
               if(sum(grp(j:n)*grp(j:n))<=tiny(1.0_dp))then
                  phi(idx)=0.0_dp
               else
                  phi(idx)=acos(max(-1.0_dp,min(1.0_dp,grp(j)/sqrt(sum(grp(j:n)*grp(j:n))))))
               end if
               idx=idx+1
            end do
            theta(i)=atan2(grp(n),grp(n-1));deallocate(grp)
         end do
      end if
      prodv=1.0_dp
      do i=1,d-2;prodv=prodv*kk(i)**real(d-1-i,dp);end do
      lambda(1)=(1.0_dp/prodv)**(1.0_dp/real(d-1,dp))
      do j=2,d-1;lambda(j)=kk(j-1)*lambda(j-1);end do
      allocate(lo(d-1,d-1),la(d-1,d-1),rr(d-1,d-1),tmp(d-1,d-1))
      call eye(rr)
      if(d>=4)then
         do m=1,d-3
            call eye(lo)
            i=d-m-1
            lo(1,1)=cos(theta(i));lo(2,2)=lo(1,1);lo(1,2)=-sin(theta(i));lo(2,1)=-lo(1,2)
            call eye(la)
            do j=1,d-m-2
               call eye(tmp)
               idx=1-j+(d-m-1)*(d-m-2)/2
               tmp(j+1,j+1)=cos(phi(idx));tmp(j+2,j+2)=tmp(j+1,j+1)
               tmp(j+1,j+2)=-sin(phi(idx));tmp(j+2,j+1)=-tmp(j+1,j+2)
               la=matmul(la,tmp)
            end do
            rr=matmul(matmul(rr,lo),la)
         end do
      end if
      call eye(lo);lo(1,1)=cos(theta(1));lo(2,2)=lo(1,1);lo(1,2)=-sin(theta(1));lo(2,1)=-lo(1,2)
      rot=matmul(rr,lo)
      if(present(info))info=0
   end subroutine

   function desag(y,mu,gam,logden) result(v)
      real(dp),intent(in)::y(:,:),mu(:),gam(:)
      logical,intent(in),optional::logden
      real(dp)::v(size(y,1)),a,yp,mp0,mp1,mp,l
      real(dp),allocatable::lambda(:),rot(:,:),basis(:,:),pmap(:,:),zinv(:),q(:)
      integer::d,p,i,j,info
      logical::ll
      d=size(mu);p=d-1;ll=.false.;if(present(logden))ll=logden
      if(size(y,2)/=d)then;v=merge(-huge(1.0_dp),0.0_dp,ll);return;end if
      call esag_parameters(gam,lambda,rot,info)
      if(info/=0 .or. size(lambda)/=p)then;v=merge(-huge(1.0_dp),0.0_dp,ll);return;end if
      allocate(basis(d,d),pmap(d,d),zinv(d),q(d));call esag_onb(mu,basis)
      pmap(:,1:p)=matmul(basis(:,1:p),rot);pmap(:,d)=basis(:,d)
      zinv(1:p)=1.0_dp/lambda;zinv(d)=1.0_dp
      do i=1,size(y,1)
         q=matmul(transpose(pmap),y(i,:));yp=sum(zinv*q*q)
         a=dot_product(y(i,:),mu)/sqrt(max(yp,tiny(1.0_dp)))
         mp0=1.0_dp;mp1=a*normal_cdf(a)+normal_pdf(a)
         if(p==1)then
            mp=mp1
         else
            mp=(1.0_dp+a*a)*normal_cdf(a)+a*normal_pdf(a)
            do j=3,p;l=mp;mp=a*mp+real(j-1,dp)*mp1;mp1=l;end do
         end if
         l=-0.5_dp*real(p,dp)*log(2.0_dp*pi)-0.5_dp*real(d,dp)*log(max(yp,tiny(1.0_dp))) &
           +0.5_dp*(a*a-sum(mu*mu))+log(max(mp,tiny(1.0_dp)))
         v(i)=merge(l,exp(l),ll)
      end do
   end function

   function resag(n,mu,gam) result(y)
      integer,intent(in)::n
      real(dp),intent(in)::mu(:),gam(:)
      real(dp)::y(n,size(mu)),z(size(mu)),q(size(mu)),nr
      real(dp),allocatable::lambda(:),rot(:,:),basis(:,:),pmap(:,:)
      integer::d,p,i,j,info
      d=size(mu);p=d-1;call esag_parameters(gam,lambda,rot,info)
      allocate(basis(d,d),pmap(d,d));call esag_onb(mu,basis)
      if(info==0)then;pmap(:,1:p)=matmul(basis(:,1:p),rot);pmap(:,d)=basis(:,d);end if
      do i=1,n
         do j=1,d;q(j)=randn();end do
         if(info==0)then
            q(1:p)=q(1:p)*sqrt(lambda);z=mu+matmul(pmap,q)
         else
            z=mu+q
         end if
         nr=sqrt(sum(z*z));y(i,:)=z/max(nr,tiny(1.0_dp))
      end do
   end function

   subroutine esag_onb(mu,basis)
      real(dp),intent(in)::mu(:);real(dp),intent(out)::basis(size(mu),size(mu))
      real(dp)::u(size(mu),size(mu)),v(size(mu),size(mu)),nr
      integer::d,i,j
      d=size(mu);u=0.0_dp;u(:,1)=mu
      if(d>=2)u(:,2)=[-mu(2),mu(1),(0.0_dp,j=3,d)]
      do i=3,d
         u(1:i-1,i)=mu(1:i-1)*mu(i);u(i,i)=-sum(mu(1:i-1)**2)
      end do
      do j=1,d
         nr=sqrt(sum(u(:,j)*u(:,j)))
         if(nr>sqrt(tiny(1.0_dp)))then
            v(:,j)=u(:,j)/nr
         else
            v(:,j)=0.0_dp;v(j,j)=1.0_dp
         end if
      end do
      basis(:,1:d-1)=v(:,2:d);basis(:,d)=v(:,1)
   end subroutine

   subroutine eye(a)
      real(dp),intent(out)::a(:,:);integer::i;a=0.0_dp;do i=1,min(size(a,1),size(a,2));a(i,i)=1.0_dp;end do
   end subroutine

   real(dp) function randn() result(z)
      real(dp)::u1,u2;call random_number(u1);call random_number(u2);z=sqrt(-2.0_dp*log(max(u1,tiny(1.0_dp))))*cos(2.0_dp*pi*u2)
   end function
end module directional_esag
