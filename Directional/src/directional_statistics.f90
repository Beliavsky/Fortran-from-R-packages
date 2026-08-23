module directional_statistics
   use directional_kinds, only : dp, pi
   implicit none
   private
   type, public :: circ_summary_result
      real(dp) :: mean_direction=0, mrl=0, circvariance=0, circstd=0, kappa=0, loglik=0
   end type
   public :: circ_summary, circ_cor2, spher_cor_rsq, mean_direction, median_direction
   public :: distance_cor, circ_dcor, spher_dcor
contains
   pure real(dp) function mean_direction(u) result(mu)
      real(dp),intent(in)::u(:);mu=atan2(sum(sin(u)),sum(cos(u)));if(mu<0)mu=mu+2*pi
   end function
   function circ_summary(u,rads) result(res)
      real(dp),intent(in)::u(:)
      logical,intent(in),optional::rads
      type(circ_summary_result)::res
      real(dp)::x(size(u)),c,s,r,k,a1,der
      integer::it
      logical::rr
      rr=.false.
      if(present(rads))rr=rads
      x=u
      if(.not.rr)x=x*pi/180
      c=sum(cos(x))/size(x)
      s=sum(sin(x))/size(x)
      r=sqrt(c*c+s*s)
      res%mean_direction=atan2(s,c)
      if(res%mean_direction<0)res%mean_direction=res%mean_direction+2*pi
      res%mrl=r
      res%circvariance=1-r
      res%circstd=sqrt(max(0.0_dp,-2*log(max(r,tiny(1.0_dp)))))
      k=max(1e-8_dp,r*(2-r*r)/max(1e-8_dp,1-r*r))
      do it=1,100
      a1=i1_over_i0(k)
      der=1-a1*a1-a1/max(k,1e-12_dp)
      if(abs(der)<1e-12_dp)exit
      k=max(0.0_dp,k-(a1-r)/der)
      end do
      res%kappa=k
      res%loglik=k*sum(cos(x-res%mean_direction))-size(x)*(log(2*pi)+log_i0_local(k))
      if(.not.rr)res%mean_direction=res%mean_direction*180/pi
   end function
   pure real(dp) function i1_over_i0(x) result(r)
      real(dp),intent(in)::x;real(dp)::h
      h=1e-5_dp*max(1.0_dp,abs(x));r=(log_i0_local(x+h)-log_i0_local(max(0.0_dp,x-h)))/(merge(2*h,h,x>h))
   end function
   pure real(dp) function log_i0_local(x) result(v)
      real(dp),intent(in)::x
      real(dp)::term,s
      integer::j
      term=1
      s=1
      do j=1,200
      term=term*(x*x/4)/(real(j,dp)**2)
      s=s+term
      if(abs(term)<1e-15_dp*s)exit
      end do
      v=log(s)
   end function
   pure real(dp) function cor(a,b) result(r)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::aa(size(a)),bb(size(b)),sa,sb
      aa=a-sum(a)/size(a)
      bb=b-sum(b)/size(b)
      sa=sqrt(sum(aa*aa))
      sb=sqrt(sum(bb*bb))
      if(sa*sb>0)then
      r=sum(aa*bb)/(sa*sb)
      else
      r=0
      end if
   end function
   function circ_cor2(theta,phi,rads) result(out)
      real(dp),intent(in)::theta(:),phi(:)
      logical,intent(in),optional::rads
      real(dp)::out(2),t(size(theta)),p(size(phi)),rcc,rcs,rss,rsc,r1,r2,up,down
      logical::rr
      rr=.false.
      if(present(rads))rr=rads
      t=theta
      p=phi
      if(.not.rr)then
      t=t*pi/180
      p=p*pi/180
      end if
      rcc=cor(cos(t),cos(p))
      rcs=cor(cos(t),sin(p))
      rss=cor(sin(t),sin(p))
      rsc=cor(sin(t),cos(p))
      r1=cor(cos(t),sin(t))
      r2=cor(cos(p),sin(p))
      up=rcc*rcc+rcs*rcs+rsc*rsc+rss*rss+2*(rcc*rss+rcs*rsc)*r1*r2-2*(rcc*rcs+rsc*rss)*r2-2*(rcc*rsc+rcs*rss)*r1
      down=(1-r1*r1)*(1-r2*r2)
      out=[up/down,-1.0_dp]
   end function
   function spher_cor_rsq(x,y) result(rsq)
      real(dp),intent(in)::x(:,:),y(:,:)
      real(dp)::rsq,xc(size(x,1),size(x,2)),yc(size(y,1),size(y,2))
      integer::i,j
      rsq=0
      xc=x
      yc=y
      do j=1,size(x,2)
      xc(:,j)=xc(:,j)-sum(xc(:,j))/size(x,1)
      end do
      do j=1,size(y,2)
      yc(:,j)=yc(:,j)-sum(yc(:,j))/size(y,1)
      end do
      do i=1,size(x,2);do j=1,size(y,2);rsq=rsq+cor(xc(:,i),yc(:,j))**2;end do;end do;rsq=rsq/min(size(x,2),size(y,2))
   end function
   function median_direction(x) result(m)
      real(dp),intent(in)::x(:,:)
      real(dp)::m(size(x,2)),score,best
      integer::i,j,ib
      best=huge(1.0_dp)
      ib=1
      do i=1,size(x,1)
      score=0
      do j=1,size(x,1)
      score=score+acos(max(-1.0_dp,min(1.0_dp,dot_product(x(i,:),x(j,:)))))
      end do
      if(score<best)then
      best=score
      ib=i
      end if
      end do
      m=x(ib,:)
   end function
   function distance_cor(x,y) result(r)
      real(dp),intent(in)::x(:,:),y(:,:)
      real(dp)::r,a(size(x,1),size(x,1)),b(size(y,1),size(y,1))
      real(dp)::ar(size(x,1)),br(size(y,1)),am,bm,dc2,dvx,dvy
      integer::i,j,n
      n=size(x,1)
      if(size(y,1)/=n .or. n<2)then;r=0.0_dp;return;end if
      do i=1,n
         do j=1,n
            a(i,j)=sqrt(sum((x(i,:)-x(j,:))**2))
            b(i,j)=sqrt(sum((y(i,:)-y(j,:))**2))
         end do
      end do
      ar=sum(a,dim=2)/real(n,dp);br=sum(b,dim=2)/real(n,dp)
      am=sum(ar)/real(n,dp);bm=sum(br)/real(n,dp)
      do i=1,n
         do j=1,n
            a(i,j)=a(i,j)-ar(i)-ar(j)+am
            b(i,j)=b(i,j)-br(i)-br(j)+bm
         end do
      end do
      dc2=sum(a*b)/real(n*n,dp);dvx=sum(a*a)/real(n*n,dp);dvy=sum(b*b)/real(n*n,dp)
      if(dvx<=tiny(1.0_dp) .or. dvy<=tiny(1.0_dp))then
         r=0.0_dp
      else
         r=sqrt(max(0.0_dp,dc2)/sqrt(dvx*dvy))
      end if
   end function

   function circ_dcor(theta,phi,rads) result(r)
      real(dp),intent(in)::theta(:),phi(:);logical,intent(in),optional::rads
      real(dp)::r,t(size(theta)),p(size(phi)),xt(size(theta),2),xp(size(phi),2)
      logical::rr
      rr=.false.;if(present(rads))rr=rads;t=theta;p=phi
      if(.not.rr)then;t=t*pi/180.0_dp;p=p*pi/180.0_dp;end if
      xt(:,1)=cos(t);xt(:,2)=sin(t);xp(:,1)=cos(p);xp(:,2)=sin(p)
      r=distance_cor(xt,xp)
   end function

   function spher_dcor(x,y) result(r)
      real(dp),intent(in)::x(:,:),y(:,:);real(dp)::r
      r=distance_cor(x,y)
   end function

end module directional_statistics
