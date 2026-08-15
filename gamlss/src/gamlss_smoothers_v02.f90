! Additional matrix-first smoothers corresponding to fp(), lo(), pvc(), and pbm().
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_smoothers_v02
   use gamlss_kinds, only : dp
   use gamlss_linalg, only : weighted_least_squares, penalized_weighted_least_squares, invert_matrix
   use gamlss_smoothers, only : p_spline_spec_t, fit_p_spline_basis, predict_p_spline_basis, difference_matrix
   implicit none
   private
   public :: fp_spec_t, fractional_polynomial_basis, select_fractional_polynomial
   public :: predict_fractional_polynomial
   public :: loess_spec_t, fit_loess, predict_loess
   public :: varying_coefficient_p_spline
   public :: monotone_spline_result_t, fit_monotone_p_spline, predict_monotone_p_spline

   type,public :: fp_spec_t
      integer :: npoly=0
      real(dp) :: powers(3)=0.0_dp
      real(dp) :: shift=0.0_dp
      real(dp) :: scale=1.0_dp
      real(dp),allocatable :: coefficients(:)
      real(dp) :: deviance=huge(1.0_dp)
      real(dp) :: edf=0.0_dp
   end type fp_spec_t

   type,public :: loess_spec_t
      real(dp),allocatable :: x(:,:),y(:),weights(:)
      real(dp),allocatable :: center(:),scale(:)
      real(dp) :: span=0.75_dp
      integer :: degree=2
      logical :: normalize=.true.
      real(dp) :: edf=0.0_dp
   end type loess_spec_t

   type,public :: monotone_spline_result_t
      type(p_spline_spec_t) :: spec
      real(dp),allocatable :: coefficients(:),fitted(:)
      real(dp) :: lambda=1.0_dp
      real(dp) :: edf=0.0_dp
      integer :: iterations=0
      integer :: status=0
      logical :: converged=.false.
   end type monotone_spline_result_t
contains

   subroutine fractional_polynomial_basis(x,powers,basis,shift,scale,status)
      real(dp),intent(in)::x(:),powers(:)
      real(dp),allocatable,intent(out)::basis(:,:)
      real(dp),intent(in),optional::shift,scale
      integer,intent(out),optional::status
      real(dp)::sh,sc,xx,prev
      integer::i,j,istat
      sh=0.0_dp;sc=1.0_dp
      if(present(shift))sh=shift
      if(present(scale))sc=scale
      istat=0
      if(sc<=0.0_dp .or. size(powers)<1 .or. size(powers)>3)istat=1
      if(istat==0 .and. any(x+sh<=0.0_dp))istat=2
      if(istat/=0)then;allocate(basis(0,0));if(present(status))status=istat;return;end if
      allocate(basis(size(x),size(powers)))
      do i=1,size(x)
         xx=(x(i)+sh)/sc
         prev=1.0_dp
         do j=1,size(powers)
            if(j>1)then
               if(abs(powers(j)-powers(j-1))<1.0e-14_dp)then
                  basis(i,j)=log(xx)*prev
               else if(abs(powers(j))<1.0e-14_dp)then
                  basis(i,j)=log(xx)
               else
                  basis(i,j)=xx**powers(j)
               end if
            else if(abs(powers(j))<1.0e-14_dp)then
               basis(i,j)=log(xx)
            else
               basis(i,j)=xx**powers(j)
            end if
            prev=basis(i,j)
         end do
      end do
      if(present(status))status=0
   end subroutine fractional_polynomial_basis

   subroutine select_fractional_polynomial(x,y,w,npoly,spec,fitted,status,shift,scale)
      real(dp),intent(in)::x(:),y(:),w(:)
      integer,intent(in)::npoly
      type(fp_spec_t),intent(out)::spec
      real(dp),allocatable,intent(out)::fitted(:)
      integer,intent(out),optional::status
      real(dp),intent(in),optional::shift,scale
      real(dp),parameter::pgrid(8)=[-2.0_dp,-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp,2.0_dp,3.0_dp]
      real(dp),allocatable::b(:,:),design(:,:),beta(:),cov(:,:)
      real(dp)::sh,sc,best,dev
      real(dp)::pp(3)
      integer::i,j,k,istat
      if(size(y)/=size(x).or.size(w)/=size(x).or.npoly<1.or.npoly>3)then
         allocate(fitted(0));if(present(status))status=1;return
      end if
      call fp_scale(x,sh,sc,istat)
      if(present(shift))sh=shift;if(present(scale))sc=scale
      if(istat/=0 .or. sc<=0.0_dp .or. any(x+sh<=0.0_dp))then
         allocate(fitted(0));if(present(status))status=2;return
      end if
      best=huge(1.0_dp);pp=0.0_dp
      if(npoly==1)then
         do i=1,8
            call try_fit([pgrid(i)])
         end do
      else if(npoly==2)then
         do i=1,8;do j=i,8
            call try_fit([pgrid(i),pgrid(j)])
         end do;end do
      else
         do i=1,8;do j=i,8;do k=j,8
            call try_fit([pgrid(i),pgrid(j),pgrid(k)])
         end do;end do;end do
      end if
      if(best>=huge(1.0_dp)/2.0_dp)then
         allocate(fitted(0));if(present(status))status=3;return
      end if
      spec%npoly=npoly;spec%powers=0.0_dp;spec%powers(1:npoly)=pp(1:npoly)
      spec%shift=sh;spec%scale=sc;spec%deviance=best;spec%edf=real(npoly+1,dp)
      call fractional_polynomial_basis(x,pp(1:npoly),b,sh,sc,istat)
      allocate(design(size(x),npoly+1));design(:,1)=1.0_dp;design(:,2:)=b
      call weighted_least_squares(design,y,w,beta,cov,istat)
      spec%coefficients=beta;allocate(fitted(size(x)));fitted=matmul(design,beta)
      if(present(status))status=istat
   contains
      subroutine try_fit(powers)
         real(dp),intent(in)::powers(:)
         real(dp),allocatable::bb(:,:),xx(:,:),bt(:),cc(:,:)
         integer::st
         call fractional_polynomial_basis(x,powers,bb,sh,sc,st);if(st/=0)return
         allocate(xx(size(x),size(powers)+1));xx(:,1)=1.0_dp;xx(:,2:)=bb
         call weighted_least_squares(xx,y,w,bt,cc,st);if(st/=0)return
         dev=sum(w*(y-matmul(xx,bt))**2)
         if(dev<best)then;best=dev;pp=0.0_dp;pp(1:size(powers))=powers;end if
      end subroutine try_fit
   end subroutine select_fractional_polynomial

   subroutine predict_fractional_polynomial(x,spec,pred,status)
      real(dp),intent(in)::x(:)
      type(fp_spec_t),intent(in)::spec
      real(dp),allocatable,intent(out)::pred(:)
      integer,intent(out),optional::status
      real(dp),allocatable::b(:,:),xx(:,:)
      integer::istat
      if(spec%npoly<1.or..not.allocated(spec%coefficients))then
         allocate(pred(0));if(present(status))status=1;return
      end if
      call fractional_polynomial_basis(x,spec%powers(1:spec%npoly),b,spec%shift,spec%scale,istat)
      if(istat/=0)then;allocate(pred(0));if(present(status))status=istat;return;end if
      allocate(xx(size(x),spec%npoly+1));xx(:,1)=1.0_dp;xx(:,2:)=b
      allocate(pred(size(x)));pred=matmul(xx,spec%coefficients)
      if(present(status))status=0
   end subroutine predict_fractional_polynomial

   subroutine fp_scale(x,shift,scale,status)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::shift,scale
      integer,intent(out)::status
      real(dp)::mindiff,rng,d
      integer::i,j
      status=0;rng=maxval(x)-minval(x)
      if(rng<=0.0_dp)then;shift=0.0_dp;scale=1.0_dp;status=1;return;end if
      shift=0.0_dp
      if(minval(x)<=0.0_dp)then
         mindiff=huge(1.0_dp)
         do i=1,size(x)-1;do j=i+1,size(x)
            d=abs(x(i)-x(j));if(d>0.0_dp)mindiff=min(mindiff,d)
         end do;end do
         if(mindiff>=huge(1.0_dp)/2.0_dp)then;status=1;return;end if
         shift=mindiff-minval(x)
      end if
      scale=10.0_dp**(sign(1.0_dp,log10(rng))*real(int(abs(log10(rng))),dp))
      if(.not.(scale>0.0_dp))scale=1.0_dp
   end subroutine fp_scale

   subroutine fit_loess(x,y,w,spec,fitted,status,span,degree,normalize,target_edf)
      real(dp),intent(in)::x(:,:),y(:),w(:)
      type(loess_spec_t),intent(out)::spec
      real(dp),allocatable,intent(out)::fitted(:)
      integer,intent(out),optional::status
      real(dp),intent(in),optional::span,target_edf
      integer,intent(in),optional::degree
      logical,intent(in),optional::normalize
      real(dp),allocatable::xs(:,:),trial(:)
      real(dp)::sp,best_sp,best_delta,edf,delta,cand
      integer::deg,istat,j
      logical::norm
      if(size(x,1)/=size(y).or.size(w)/=size(y).or.size(x,2)<1.or.size(x,2)>2)then
         allocate(fitted(0));if(present(status))status=1;return
      end if
      sp=0.75_dp;if(present(span))sp=span;sp=min(1.0_dp,max(0.05_dp,sp))
      deg=2;if(present(degree))deg=degree;deg=max(0,min(2,deg))
      norm=.true.;if(present(normalize))norm=normalize
      spec%x=x;spec%y=y;spec%weights=w;spec%degree=deg;spec%normalize=norm
      allocate(spec%center(size(x,2)),spec%scale(size(x,2)));spec%center=0.0_dp;spec%scale=1.0_dp
      xs=x
      if(norm)then
         do j=1,size(x,2)
            spec%center(j)=sum(x(:,j))/real(size(x,1),dp)
            spec%scale(j)=sqrt(max(1.0e-12_dp,sum((x(:,j)-spec%center(j))**2)/real(max(1,size(x,1)-1),dp)))
            xs(:,j)=(x(:,j)-spec%center(j))/spec%scale(j)
         end do
      end if
      if(present(target_edf))then
         best_delta=huge(1.0_dp);best_sp=sp
         do j=0,24
            cand=0.08_dp+0.92_dp*real(j,dp)/24.0_dp
            call loess_compute(xs,y,w,xs,cand,deg,trial,edf,istat,.true.)
            if(istat/=0)cycle
            delta=abs(edf-target_edf)
            if(delta<best_delta)then;best_delta=delta;best_sp=cand;end if
         end do
         sp=best_sp
      end if
      spec%span=sp
      call loess_compute(xs,y,w,xs,sp,deg,fitted,edf,istat,.true.);spec%edf=edf
      if(present(status))status=istat
   end subroutine fit_loess

   subroutine predict_loess(xnew,spec,pred,status)
      real(dp),intent(in)::xnew(:,:)
      type(loess_spec_t),intent(in)::spec
      real(dp),allocatable,intent(out)::pred(:)
      integer,intent(out),optional::status
      real(dp),allocatable::xt(:,:),xn(:,:)
      real(dp)::edf
      integer::j,istat
      if(size(xnew,2)/=size(spec%x,2))then;allocate(pred(0));if(present(status))status=1;return;end if
      xt=spec%x;xn=xnew
      if(spec%normalize)then
         do j=1,size(xt,2)
            xt(:,j)=(xt(:,j)-spec%center(j))/spec%scale(j)
            xn(:,j)=(xn(:,j)-spec%center(j))/spec%scale(j)
         end do
      end if
      call loess_compute(xt,spec%y,spec%weights,xn,spec%span,spec%degree,pred,edf,istat,.false.)
      if(present(status))status=istat
   end subroutine predict_loess

   subroutine loess_compute(x,y,w,target,span,degree,pred,edf,status,training)
      real(dp),intent(in)::x(:,:),y(:),w(:),target(:,:),span
      integer,intent(in)::degree
      real(dp),allocatable,intent(out)::pred(:)
      real(dp),intent(out)::edf
      integer,intent(out)::status
      logical,intent(in)::training
      real(dp),allocatable::dist(:),wl(:),phi(:,:),beta(:),cov(:,:)
      real(dp)::h,u,lev
      integer::n,d,k,i,j,istat,nearest
      n=size(x,1);d=size(x,2);k=max(local_terms(d,degree)+1,int(ceiling(span*real(n,dp))))
      k=min(n,k);allocate(pred(size(target,1)),dist(n),wl(n));pred=0.0_dp;edf=0.0_dp;status=0
      do i=1,size(target,1)
         do j=1,n;dist(j)=sqrt(sum((x(j,:)-target(i,:))**2));end do
         h=kth_value(dist,k)
         if(h<=1.0e-14_dp)then
            wl=merge(w,0.0_dp,dist<=1.0e-14_dp)
         else
            do j=1,n
               u=dist(j)/h
               if(u<1.0_dp)then;wl(j)=w(j)*(1.0_dp-u**3)**3;else;wl(j)=0.0_dp;end if
            end do
         end if
         call local_design(x,target(i,:),degree,phi)
         call weighted_least_squares(phi,y,wl,beta,cov,istat)
         if(istat/=0)then
            ! Fall back to the nearest observation for degenerate local neighborhoods.
            nearest=minloc(dist,dim=1);pred(i)=y(nearest);cycle
         end if
         pred(i)=beta(1)
         if(training .and. i<=n)then
            lev=wl(i)*cov(1,1);edf=edf+max(0.0_dp,min(1.0_dp,lev))
         end if
      end do
   end subroutine loess_compute

   integer pure function local_terms(d,degree) result(p)
      integer,intent(in)::d,degree
      p=1
      if(degree>=1)p=p+d
      if(degree>=2)p=p+d*(d+1)/2
   end function local_terms

   subroutine local_design(x,target,degree,phi)
      real(dp),intent(in)::x(:,:),target(:)
      integer,intent(in)::degree
      real(dp),allocatable,intent(out)::phi(:,:)
      real(dp),allocatable::dx(:,:)
      integer::i,j,k,col,d
      d=size(x,2);allocate(phi(size(x,1),local_terms(d,degree)),dx(size(x,1),d))
      do j=1,d;dx(:,j)=x(:,j)-target(j);end do
      phi(:,1)=1.0_dp;col=1
      if(degree>=1)then
         do j=1,d;col=col+1;phi(:,col)=dx(:,j);end do
      end if
      if(degree>=2)then
         do j=1,d;do k=j,d;col=col+1;phi(:,col)=dx(:,j)*dx(:,k);end do;end do
      end if
   end subroutine local_design

   real(dp) function kth_value(a,k) result(v)
      real(dp),intent(in)::a(:)
      integer,intent(in)::k
      real(dp),allocatable::b(:)
      real(dp)::tmp
      integer::i,j
      b=a
      do i=1,k
         j=i-1+minloc(b(i:),dim=1);tmp=b(i);b(i)=b(j);b(j)=tmp
      end do
      v=b(k)
   end function kth_value

   subroutine varying_coefficient_p_spline(x,z,basis,penalty,spec,df,degree,order,status)
      real(dp),intent(in)::x(:),z(:)
      real(dp),allocatable,intent(out)::basis(:,:),penalty(:,:)
      type(p_spline_spec_t),intent(out)::spec
      integer,intent(in),optional::df,degree,order
      integer,intent(out),optional::status
      real(dp),allocatable::b(:,:)
      integer::i,istat
      if(size(z)/=size(x))then;allocate(basis(0,0),penalty(0,0));if(present(status))status=1;return;end if
      call fit_p_spline_basis(x,spec,b,df,degree,order,istat)
      if(istat/=0)then;allocate(basis(0,0),penalty(0,0));if(present(status))status=istat;return;end if
      allocate(basis(size(b,1),size(b,2)));basis=b
      do i=1,size(x);basis(i,:)=z(i)*basis(i,:);end do
      penalty=spec%penalty
      if(present(status))status=0
   end subroutine varying_coefficient_p_spline

   subroutine fit_monotone_p_spline(x,y,w,result,increasing,lambda,df,degree,order,kappa,max_iter,tol)
      real(dp),intent(in)::x(:),y(:),w(:)
      type(monotone_spline_result_t),intent(out)::result
      logical,intent(in),optional::increasing
      real(dp),intent(in),optional::lambda,kappa,tol
      integer,intent(in),optional::df,degree,order,max_iter
      real(dp),allocatable::b(:,:),d1(:,:),pen(:,:),active_pen(:,:),beta(:),cov(:,:),diff(:),old(:)
      real(dp)::kap,crit,lam
      integer::istat,it,nit,i
      logical::up,changed
      if(size(y)/=size(x).or.size(w)/=size(x))then;result%status=1;return;end if
      up=.true.;if(present(increasing))up=increasing
      lam=1.0_dp;if(present(lambda))lam=max(0.0_dp,lambda)
      kap=1.0e8_dp;if(present(kappa))kap=max(1.0_dp,kappa)
      nit=50;if(present(max_iter))nit=max_iter;crit=1.0e-7_dp;if(present(tol))crit=tol
      call fit_p_spline_basis(x,result%spec,b,df,degree,order,istat)
      if(istat/=0)then;result%status=istat;return;end if
      call difference_matrix(size(b,2),1,d1);allocate(active_pen(size(b,2),size(b,2)));active_pen=0.0_dp
      allocate(beta(size(b,2)));beta=0.0_dp
      do it=1,nit
         pen=lam*result%spec%penalty+kap*active_pen
         old=beta
         call penalized_weighted_least_squares(b,y,w,pen,beta,cov,istat)
         if(istat/=0)then;result%status=2;return;end if
         diff=matmul(d1,beta);active_pen=0.0_dp;changed=.false.
         do i=1,size(diff)
            if((up.and.diff(i)<0.0_dp).or.((.not.up).and.diff(i)>0.0_dp))then
               active_pen=active_pen+outer_row(d1(i,:));changed=.true.
            end if
         end do
         result%iterations=it
         if(maxval(abs(beta-old))<crit .and. .not.changed)then;result%converged=.true.;exit;end if
         if(maxval(abs(beta-old))<crit .and. changed)then
            ! Active set is stable if the violating set no longer changes materially.
            result%converged=.true.;exit
         end if
      end do
      result%coefficients=beta;allocate(result%fitted(size(x)));result%fitted=matmul(b,beta)
      result%lambda=lam;result%edf=penalized_edf(b,w,lam*result%spec%penalty+kap*active_pen)
      result%status=0
   end subroutine fit_monotone_p_spline

   subroutine predict_monotone_p_spline(x,result,pred,status)
      real(dp),intent(in)::x(:)
      type(monotone_spline_result_t),intent(in)::result
      real(dp),allocatable,intent(out)::pred(:)
      integer,intent(out),optional::status
      real(dp),allocatable::b(:,:)
      integer::istat
      call predict_p_spline_basis(x,result%spec,b,istat)
      if(istat/=0)then;allocate(pred(0));if(present(status))status=istat;return;end if
      allocate(pred(size(x)));pred=matmul(b,result%coefficients);if(present(status))status=0
   end subroutine predict_monotone_p_spline

   function outer_row(v) result(a)
      real(dp),intent(in)::v(:)
      real(dp)::a(size(v),size(v))
      integer::i
      do i=1,size(v);a(i,:)=v(i)*v;end do
   end function outer_row

   real(dp) function penalized_edf(x,w,pen) result(edf)
      real(dp),intent(in)::x(:,:),w(:),pen(:,:)
      real(dp),allocatable::xtwx(:,:),a(:,:),ainv(:,:),h(:,:)
      integer::i,j,k,p,istat
      p=size(x,2);allocate(xtwx(p,p));xtwx=0.0_dp
      do i=1,size(x,1);do j=1,p;do k=1,p
         xtwx(j,k)=xtwx(j,k)+w(i)*x(i,j)*x(i,k)
      end do;end do;end do
      a=xtwx+pen;call invert_matrix(a,ainv,istat)
      if(istat/=0)then;edf=real(p,dp);return;end if
      h=matmul(ainv,xtwx);edf=0.0_dp;do i=1,p;edf=edf+h(i,i);end do
   end function penalized_edf

end module gamlss_smoothers_v02
