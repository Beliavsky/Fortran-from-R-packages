module urca_breaks
   use urca_kinds, only : dp
   use urca_types, only : johansen_result
   use urca_regression, only : lm_fit_multi
   use urca_linalg, only : determinant_spd, invert_spd, invert_matrix, chol_lower, symmetric_eigen
   use urca_cointegration, only : JO_NONE, JO_TRANSITORY
   implicit none
   private
   public :: johansen_level_shift
contains
   subroutine seasonal_dummies_local(n,season,d)
      integer,intent(in)::n,season
      real(dp),allocatable,intent(out)::d(:,:)
      integer::i,j,k
      if(season<=1)then
      allocate(d(n,0))
      return
      end if
      allocate(d(n,season-1))
      d=0.0_dp
      do i=1,n
      k=mod(i-1,season)+1
      do j=1,season-1
      d(i,j)=-1.0_dp/real(season,dp)
      if(k==j)d(i,j)=d(i,j)+1.0_dp
      end do
      end do
   end subroutine seasonal_dummies_local

   function johansen_level_shift(x,trend,k,season) result(out)
      real(dp),intent(in)::x(:,:)
      logical,intent(in)::trend
      integer,intent(in)::k
      integer,intent(in),optional::season
      type(johansen_result)::out
      integer::nt,p,n,ss,ns,nrhs,i,j,l,c,best,bp,naux,info,ne,rank
      real(dp)::detv,bestdet
      real(dp),allocatable::sd(:,:),lhs(:,:),rhs(:,:),design(:,:),dt(:),b(:,:),res(:,:),sigma(:,:),rhsaux(:,:), &
         & coef(:,:),yfit(:,:),dx(:,:),z0(:,:),z1(:,:),zk(:,:)
      real(dp),allocatable::m00(:,:),m11(:,:),mkk(:,:),m01(:,:),m0k(:,:),mk0(:,:),m10(:,:),m1k(:,:),mk1(:,:),m11i(:, &
         & :),r0(:,:),rk(:,:),s00(:,:),s0k(:,:),sk0(:,:),skk(:,:),chol(:,:),ci(:,:),s00i(:,:),a(:,:),ev(:),e(:,:), &
         & vorg(:,:),v(:,:),vv(:,:),vvi(:,:),skki(:,:)
      nt=size(x,1)
      p=size(x,2)
      ss=0
      if(present(season))ss=season
      ns=max(0,ss-1)
      if(k<2.or.nt-k<max(8,p+1))then
      out%info=-1
      return
      end if
      call seasonal_dummies_local(nt,ss,sd)
      n=nt-k
      nrhs=k*p+ns+merge(1,0,trend)
      allocate(lhs(n,p),rhs(n,nrhs))
      c=0
      if(ns>0)then
      rhs(:,1:ns)=sd(k+1:nt,:)
      c=ns
      end if
      if(trend)then
      c=c+1
      do i=1,n
      rhs(i,c)=real(k+i,dp)
      end do
      end if
      do l=1,k
         rhs(:,c+(l-1)*p+1:c+l*p)=x(k+1-l:nt-l,:)
      end do
      lhs=x(k+1:nt,:)
      bestdet=huge(1.0_dp)
      best=1
      do j=1,n-1
         allocate(dt(n),design(n,2+nrhs))
         dt=0.0_dp
         dt(j+1:n)=1.0_dp
         design(:,1)=1.0_dp
         design(:,2)=dt
         design(:,3:)=rhs
         call lm_fit_multi(design,lhs,b,res,sigma,info)
         if(info==0)then
         detv=determinant_spd(matmul(transpose(res),res),info)
         if(info==0.and.detv<bestdet)then
         bestdet=detv
         best=j
         end if
         end if
         deallocate(dt,design)
      end do
      bp=best+k+1
! Shift-correction regression: intercept + break dummy + deterministic terms only.
      naux=2+ns+merge(1,0,trend)
      allocate(design(n,naux),dt(n))
      dt=0.0_dp
      dt(best+k+1:n)=1.0_dp
! R's tau.opt is best+K, and its sample dt is c(0*tau.opt,1*(N-tau.opt));
! because N here equals nt-K, cap the zero segment at N.
      dt=0.0_dp
      if(best+k<n)dt(best+k+1:n)=1.0_dp
      design(:,1)=1.0_dp
      design(:,2)=dt
      c=2
      if(ns>0)then
      design(:,c+1:c+ns)=sd(k+1:nt,:)
      c=c+ns
      end if
      if(trend)then
      c=c+1
      do i=1,n
      design(i,c)=real(k+i,dp)
      end do
      end if
      call lm_fit_multi(design,lhs,coef,res,sigma,info)
      if(info/=0)then
      out%info=100+info
      return
      end if
      allocate(yfit(nt,p))
      yfit=x
      do i=1,nt
         yfit(i,:)=yfit(i,:)-coef(1,:)
         if(i>best+k)yfit(i,:)=yfit(i,:)-coef(2,:)
         c=2
         if(ns>0)then
         yfit(i,:)=yfit(i,:)-matmul(sd(i,:),coef(c+1:c+ns,:))
         c=c+ns
         end if
         if(trend)then
         c=c+1
         yfit(i,:)=yfit(i,:)-real(i,dp)*coef(c,:)
         end if
      end do
      allocate(dx(nt-1,p))
      dx=yfit(2:nt,:)-yfit(1:nt-1,:)
      ne=nt-k
      allocate(z0(ne,p),z1(ne,(k-1)*p),zk(ne,p))
      z0=dx(k:nt-1,:)
      zk=yfit(k:nt-1,:)
      do l=1,k-1
      z1(:,(l-1)*p+1:l*p)=dx(k-l:nt-1-l,:)
      end do
      m00=matmul(transpose(z0),z0)/real(ne,dp)
      m11=matmul(transpose(z1),z1)/real(ne,dp)
      mkk=matmul(transpose(zk),zk)/real(ne,dp)
      m01=matmul(transpose(z0),z1)/real(ne,dp)
      m10=transpose(m01)
      m0k=matmul(transpose(z0),zk)/real(ne,dp)
      mk0=transpose(m0k)
      m1k=matmul(transpose(z1),zk)/real(ne,dp)
      mk1=transpose(m1k)
      call invert_spd(m11,m11i,info)
      if(info/=0)call invert_matrix(m11,m11i,info)
      if(info/=0)then
      out%info=200+info
      return
      end if
      r0=z0-matmul(z1,matmul(m11i,m10))
      rk=zk-matmul(z1,matmul(m11i,m1k))
      s00=m00-matmul(m01,matmul(m11i,m10))
      s0k=m0k-matmul(m01,matmul(m11i,m1k))
      sk0=transpose(s0k)
      skk=mkk-matmul(mk1,matmul(m11i,m1k))
      call chol_lower(skk,chol,info)
      if(info/=0)then
      out%info=300+info
      return
      end if
      call invert_matrix(chol,ci,info)
      if(info/=0)then
      out%info=301+info
      return
      end if
      call invert_spd(s00,s00i,info)
      if(info/=0)call invert_matrix(s00,s00i,info)
      if(info/=0)then
      out%info=400+info
      return
      end if
      a=matmul(ci,matmul(sk0,matmul(s00i,matmul(s0k,transpose(ci)))))
      call symmetric_eigen(a,ev,e,info,.true.)
      if(info/=0)then
      out%info=500+info
      return
      end if
      vorg=matmul(transpose(ci),e)
      v=vorg
      do j=1,p
      if(abs(v(1,j))>1e-14_dp)v(:,j)=v(:,j)/v(1,j)
      end do
      vv=matmul(transpose(v),matmul(skk,v))
      call invert_matrix(vv,vvi,info)
      if(info/=0)then
      out%info=600+info
      return
      end if
      out%w=matmul(s0k,matmul(v,vvi))
      call invert_spd(skk,skki,info)
      if(info/=0)call invert_matrix(skk,skki,info)
      if(info/=0)then
      out%info=700+info
      return
      end if
      out%pi=matmul(s0k,skki)
      out%delta=s00-matmul(s0k,matmul(v,matmul(vvi,matmul(transpose(v),sk0))))
      out%gamma=matmul(m01,m11i)-matmul(out%pi,matmul(mk1,m11i))
      allocate(out%teststat(p))
      do i=1,p
      rank=p-i
      out%teststat(i)=real(ne,dp)*sum(log(1.0_dp+max(0.0_dp,ev(rank+1:p))))
      end do
      allocate(out%critical_values(p,3))
      call shift_critical_values(p,trend,out%critical_values)
      out%x=yfit
      out%z0=z0
      out%z1=z1
      out%zk=zk
      out%lambda=ev
      out%vorg=vorg
      out%v=v
      out%r0=r0
      out%rk=rk
      out%p=p
      out%lag=k
      out%ecdet=JO_NONE
      out%spec=JO_TRANSITORY
      out%break_point=bp
      out%info=0
   end function johansen_level_shift

   subroutine shift_critical_values(p,trend,cv)
      integer,intent(in)::p
      logical,intent(in)::trend
      real(dp),intent(out)::cv(p,3)
      real(dp),parameter::yes(5,3)=reshape([5.423_dp,13.784_dp,25.931_dp,42.083_dp,61.918_dp,6.785_dp,15.826_dp, &
         & 28.455_dp,45.204_dp,65.662_dp,10.042_dp,19.854_dp,33.757_dp,51.601_dp,73.116_dp],[5,3])
      real(dp),parameter::no(5,3)=reshape([2.996_dp,10.446_dp,21.801_dp,36.903_dp,55.952_dp,4.118_dp,12.276_dp, &
         & 24.282_dp,40.067_dp,59.749_dp,6.888_dp,16.420_dp,29.467_dp,46.305_dp,67.170_dp],[5,3])
      integer::i
      cv=0.0_dp
      if(p>5)return
      do i=1,p
      if(trend)then
      cv(i,:)=yes(p-i+1,:)
      else
      cv(i,:)=no(p-i+1,:)
      end if
      end do
   end subroutine shift_critical_values
end module urca_breaks
