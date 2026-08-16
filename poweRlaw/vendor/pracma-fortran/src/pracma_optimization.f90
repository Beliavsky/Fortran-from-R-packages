! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_optimization
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use pracma_kinds, only : dp, eps_dp
   use pracma_status
   use pracma_types, only : optimization_result, qpspecial_result, linprog_result, regression_result
   use pracma_callbacks
   use pracma_differentiation, only : grad, jacobian
   use pracma_linalg, only : mldivide, identity_matrix, pinv, nearest_spd, outer_product
   use quadprog, only : solve_qp, quadprog_result => qp_result, qp_success
   implicit none
   private

   public :: fminbnd, fibsearch, nelder_mead, fminsearch, anms
   public :: hooke_jeeves, steep_descent, fminunc, fletcher_powell
   public :: gaussNewton, curvefit, lsqnonlin, lsqcurvefit, lsqnonneg
   public :: L1linreg, geo_median, qpspecial, qpsolve, quadprog
   public :: linprog, fmincon, qpsolve_projection

contains

   function fminbnd(f,a,b,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(optimization_result)::res
      real(dp)::left,right,x,w,v,u,fx,fw,fv,fu,d,e,p,q,r,tol,tol1,tol2,m
      real(dp),parameter::cgold=0.3819660112501051_dp
      integer::iter,niter
      tol=1.0e-10_dp; if(present(tolerance))tol=tolerance
      niter=500; if(present(max_iter))niter=max_iter
      allocate(res%x(1),res%history(niter)); res%history=0.0_dp
      if(b<=a)then
         res%x=ieee_value(0.0_dp,ieee_quiet_nan); res%status=pracma_invalid_argument; return
      end if
      left=a; right=b; x=left+cgold*(right-left); w=x; v=x
      fx=f(x); fw=fx; fv=fx; res%evaluations=1; d=0.0_dp; e=0.0_dp
      do iter=1,niter
         m=0.5_dp*(left+right); tol1=tol*abs(x)+sqrt(eps_dp); tol2=2.0_dp*tol1
         if(abs(x-m)<=tol2-0.5_dp*(right-left))exit
         p=0.0_dp; q=0.0_dp; r=0.0_dp
         if(abs(e)>tol1)then
            r=(x-w)*(fx-fv); q=(x-v)*(fx-fw); p=(x-v)*q-(x-w)*r
            q=2.0_dp*(q-r); if(q>0.0_dp)p=-p; q=abs(q)
            r=e; e=d
            if(abs(p)>=abs(0.5_dp*q*r) .or. p<=q*(left-x) .or. p>=q*(right-x))then
               if(x>=m)then; e=left-x; else; e=right-x; end if
               d=cgold*e
            else
               d=p/q; u=x+d
               if(u-left<tol2 .or. right-u<tol2)d=sign(tol1,m-x)
            end if
         else
            if(x>=m)then; e=left-x; else; e=right-x; end if
            d=cgold*e
         end if
         if(abs(d)>=tol1)then; u=x+d; else; u=x+sign(tol1,d); end if
         fu=f(u); res%evaluations=res%evaluations+1
         if(fu<=fx)then
            if(u>=x)then; left=x; else; right=x; end if
            v=w; fv=fw; w=x; fw=fx; x=u; fx=fu
         else
            if(u<x)then; left=u; else; right=u; end if
            if(fu<=fw .or. abs(w-x)<=eps_dp)then
               v=w; fv=fw; w=u; fw=fu
            else if(fu<=fv .or. abs(v-x)<=eps_dp .or. abs(v-w)<=eps_dp)then
               v=u; fv=fu
            end if
         end if
         res%history(iter)=fx
      end do
      res%x(1)=x; res%value=fx; res%iterations=min(iter,niter)
      res%converged=iter<=niter; res%status=merge(pracma_ok,pracma_not_converged,res%converged)
      call shrink_history(res)
   end function fminbnd

   function fibsearch(f,a,b,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(optimization_result)::res
      res=fminbnd(f,a,b,tolerance,max_iter)
   end function fibsearch

   function nelder_mead(f,x0,step,tolerance,max_iter) result(res)
      procedure(objective_function)::f
      real(dp),intent(in)::x0(:)
      real(dp),intent(in),optional::step,tolerance
      integer,intent(in),optional::max_iter
      type(optimization_result)::res
      real(dp),allocatable::simplex(:,:),values(:),centroid(:),xr(:),xe(:),xc(:)
      real(dp)::s,tol,fr,fe,fc,spreadx
      integer::n,i,iter,niter
      n=size(x0); s=0.05_dp; if(present(step))s=step
      tol=1.0e-9_dp; if(present(tolerance))tol=tolerance
      niter=max(500,200*n); if(present(max_iter))niter=max_iter
      allocate(simplex(n,n+1),values(n+1),centroid(n),xr(n),xe(n),xc(n))
      allocate(res%x(n),res%history(niter)); res%history=0.0_dp
      simplex(:,1)=x0
      do i=1,n
         simplex(:,i+1)=x0
         simplex(i,i+1)=x0(i)+s*max(1.0_dp,abs(x0(i)))
      end do
      do i=1,n+1
         values(i)=f(simplex(:,i)); res%evaluations=res%evaluations+1
      end do
      do iter=1,niter
         call sort_simplex(simplex,values)
         res%history(iter)=values(1)
         spreadx=maxval(abs(simplex-spread(simplex(:,1),2,n+1)))
         if(maxval(abs(values-values(1)))<=tol*(1.0_dp+abs(values(1))) .and. &
            spreadx<=sqrt(tol)*(1.0_dp+maxval(abs(simplex(:,1)))))exit
         centroid=sum(simplex(:,1:n),dim=2)/real(n,dp)
         xr=centroid+(centroid-simplex(:,n+1)); fr=f(xr); res%evaluations=res%evaluations+1
         if(fr<values(1))then
            xe=centroid+2.0_dp*(xr-centroid); fe=f(xe); res%evaluations=res%evaluations+1
            if(fe<fr)then; simplex(:,n+1)=xe; values(n+1)=fe
            else; simplex(:,n+1)=xr; values(n+1)=fr; end if
         else if(fr<values(n))then
            simplex(:,n+1)=xr; values(n+1)=fr
         else
            if(fr<values(n+1))then
               xc=centroid+0.5_dp*(xr-centroid)
            else
               xc=centroid+0.5_dp*(simplex(:,n+1)-centroid)
            end if
            fc=f(xc); res%evaluations=res%evaluations+1
            if(fc<min(fr,values(n+1)))then
               simplex(:,n+1)=xc; values(n+1)=fc
            else
               do i=2,n+1
                  simplex(:,i)=simplex(:,1)+0.5_dp*(simplex(:,i)-simplex(:,1))
                  values(i)=f(simplex(:,i)); res%evaluations=res%evaluations+1
               end do
            end if
         end if
      end do
      call sort_simplex(simplex,values)
      res%x=simplex(:,1); res%value=values(1); res%iterations=min(iter,niter)
      res%converged=iter<=niter; res%status=merge(pracma_ok,pracma_not_converged,res%converged)
      call shrink_history(res)
   end function nelder_mead

   function fminsearch(f,x0,step,tolerance,max_iter) result(res)
      procedure(objective_function)::f
      real(dp),intent(in)::x0(:)
      real(dp),intent(in),optional::step,tolerance
      integer,intent(in),optional::max_iter
      type(optimization_result)::res
      res=nelder_mead(f,x0,step,tolerance,max_iter)
   end function fminsearch

   function anms(f,x0,step,tolerance,max_iter) result(res)
      procedure(objective_function)::f
      real(dp),intent(in)::x0(:)
      real(dp),intent(in),optional::step,tolerance
      integer,intent(in),optional::max_iter
      type(optimization_result)::res
      res=nelder_mead(f,x0,step,tolerance,max_iter)
   end function anms

   subroutine sort_simplex(simplex,values)
      real(dp),intent(inout)::simplex(:,:),values(:)
      real(dp),allocatable::col(:)
      real(dp)::v
      integer::i,j,k
      allocate(col(size(simplex,1)))
      do i=1,size(values)-1
         k=i
         do j=i+1,size(values)
            if(values(j)<values(k))k=j
         end do
         if(k/=i)then
            v=values(i); values(i)=values(k); values(k)=v
            col=simplex(:,i); simplex(:,i)=simplex(:,k); simplex(:,k)=col
         end if
      end do
   end subroutine sort_simplex

   function hooke_jeeves(f,x0,step,tolerance,max_iter) result(res)
      procedure(objective_function)::f
      real(dp),intent(in)::x0(:)
      real(dp),intent(in),optional::step,tolerance
      integer,intent(in),optional::max_iter
      type(optimization_result)::res
      real(dp),allocatable::x(:),xb(:),xn(:),delta(:)
      real(dp)::fx,fb,fn,tol,s
      integer::n,i,iter,niter
      n=size(x0); s=0.5_dp; if(present(step))s=step
      tol=1.0e-8_dp; if(present(tolerance))tol=tolerance
      niter=5000; if(present(max_iter))niter=max_iter
      allocate(x(n),xb(n),xn(n),delta(n),res%x(n),res%history(niter))
      x=x0; xb=x0; delta=s*max(1.0_dp,abs(x0)); fb=f(xb); res%evaluations=1
      do iter=1,niter
         xn=xb; fn=fb
         do i=1,n
            x=xn; x(i)=x(i)+delta(i); fx=f(x); res%evaluations=res%evaluations+1
            if(fx<fn)then
               xn=x; fn=fx
            else
               x=xn; x(i)=x(i)-delta(i); fx=f(x); res%evaluations=res%evaluations+1
               if(fx<fn)then; xn=x; fn=fx; end if
            end if
         end do
         if(fn<fb)then
            x=xn+(xn-xb); xb=xn; fb=fn
            fx=f(x); res%evaluations=res%evaluations+1
            if(fx<fb)then; xb=x; fb=fx; end if
         else
            delta=0.5_dp*delta
         end if
         res%history(iter)=fb
         if(maxval(abs(delta))<=tol*(1.0_dp+maxval(abs(xb))))exit
      end do
      res%x=xb; res%value=fb; res%iterations=min(iter,niter)
      res%converged=iter<=niter; res%status=merge(pracma_ok,pracma_not_converged,res%converged)
      call shrink_history(res)
   end function hooke_jeeves

   function steep_descent(f,x0,tolerance,max_iter) result(res)
      procedure(objective_function)::f
      real(dp),intent(in)::x0(:)
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(optimization_result)::res
      real(dp),allocatable::x(:),g(:),xn(:)
      real(dp)::fx,fn,alpha,tol
      integer::iter,niter,istat
      tol=1.0e-8_dp; if(present(tolerance))tol=tolerance
      niter=5000; if(present(max_iter))niter=max_iter
      allocate(x(size(x0)),g(size(x0)),xn(size(x0)),res%x(size(x0)),res%history(niter))
      x=x0; fx=f(x); res%evaluations=1
      do iter=1,niter
         call grad(f,x,g,status=istat); res%evaluations=res%evaluations+2*size(x)
         if(sqrt(sum(g*g))<=tol)exit
         alpha=1.0_dp
         do
            xn=x-alpha*g; fn=f(xn); res%evaluations=res%evaluations+1
            if(fn<=fx-1.0e-4_dp*alpha*sum(g*g) .or. alpha<=1.0e-12_dp)exit
            alpha=0.5_dp*alpha
         end do
         x=xn; fx=fn; res%history(iter)=fx
         if(alpha<=1.0e-12_dp)exit
      end do
      res%x=x; res%value=fx; res%iterations=min(iter,niter)
      res%converged=sqrt(sum(g*g))<=tol; res%status=merge(pracma_ok,pracma_not_converged,res%converged)
      call shrink_history(res)
   end function steep_descent

   function fminunc(f,x0,tolerance,max_iter) result(res)
      procedure(objective_function)::f
      real(dp),intent(in)::x0(:)
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(optimization_result)::res
      real(dp),allocatable::x(:),g(:),gnew(:),h(:,:),p(:),xn(:),svec(:),yvec(:)
      real(dp)::fx,fn,alpha,tol,ys,rho
      integer::n,iter,niter,istat
      n=size(x0); tol=1.0e-8_dp; if(present(tolerance))tol=tolerance
      niter=max(500,100*n); if(present(max_iter))niter=max_iter
      allocate(x(n),g(n),gnew(n),h(n,n),p(n),xn(n),svec(n),yvec(n),res%x(n),res%history(niter))
      x=x0; h=identity_matrix(n); fx=f(x); call grad(f,x,g,status=istat); res%evaluations=1+2*n
      do iter=1,niter
         if(sqrt(sum(g*g))<=tol*(1.0_dp+abs(fx)))exit
         p=-matmul(h,g)
         if(dot_product(p,g)>=0.0_dp)p=-g
         alpha=1.0_dp
         do
            xn=x+alpha*p; fn=f(xn); res%evaluations=res%evaluations+1
            if(fn<=fx+1.0e-4_dp*alpha*dot_product(g,p) .or. alpha<=1.0e-12_dp)exit
            alpha=0.5_dp*alpha
         end do
         call grad(f,xn,gnew,status=istat); res%evaluations=res%evaluations+2*n
         svec=xn-x; yvec=gnew-g; ys=dot_product(yvec,svec)
         if(ys>sqrt(eps_dp)*sqrt(sum(yvec*yvec)*sum(svec*svec)))then
            rho=1.0_dp/ys
            h=matmul(identity_matrix(n)-rho*outer_product(svec,yvec), &
               matmul(h,identity_matrix(n)-rho*outer_product(yvec,svec)))+rho*outer_product(svec,svec)
         else
            h=identity_matrix(n)
         end if
         x=xn; fx=fn; g=gnew; res%history(iter)=fx
         if(alpha<=1.0e-12_dp)exit
      end do
      res%x=x; res%value=fx; res%iterations=min(iter,niter)
      res%converged=sqrt(sum(g*g))<=tol*(1.0_dp+abs(fx))
      res%status=merge(pracma_ok,pracma_not_converged,res%converged)
      call shrink_history(res)
   end function fminunc

   function fletcher_powell(f,x0,tolerance,max_iter) result(res)
      procedure(objective_function)::f
      real(dp),intent(in)::x0(:)
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(optimization_result)::res
      res=fminunc(f,x0,tolerance,max_iter)
   end function fletcher_powell

   function gaussNewton(residual,x0,m,tolerance,max_iter) result(res)
      procedure(vector_function)::residual
      real(dp),intent(in)::x0(:)
      integer,intent(in)::m
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(regression_result)::res
      real(dp),allocatable::x(:),r(:),rn(:),j(:,:),dx(:),xn(:)
      real(dp)::rss,rssn,alpha,tol
      integer::n,iter,niter,istat
      n=size(x0); tol=1.0e-8_dp; if(present(tolerance))tol=tolerance
      niter=500; if(present(max_iter))niter=max_iter
      allocate(x(n),r(m),rn(m),j(m,n),dx(n),xn(n))
      allocate(res%coefficients(n),res%fitted(m),res%residuals(m))
      x=x0; call residual(x,r); rss=sum(r*r)
      do iter=1,niter
         call jacobian(residual,x,m,j,status=istat)
         dx=mldivide(j,-r,istat)
         if(istat/=pracma_ok)exit
         alpha=1.0_dp
         do
            xn=x+alpha*dx; call residual(xn,rn); rssn=sum(rn*rn)
            if(rssn<rss .or. alpha<=1.0e-8_dp)exit
            alpha=0.5_dp*alpha
         end do
         x=xn; r=rn; rss=rssn
         if(sqrt(sum((alpha*dx)**2))<=tol*(1.0_dp+sqrt(sum(x*x))))exit
      end do
      res%coefficients=x; res%residuals=r; res%fitted=-r; res%rss=rss
      res%iterations=min(iter,niter); res%converged=iter<=niter
      res%status=merge(pracma_ok,pracma_not_converged,res%converged)
   end function gaussNewton

   function lsqnonlin(residual,x0,m,tolerance,max_iter) result(res)
      procedure(vector_function)::residual
      real(dp),intent(in)::x0(:)
      integer,intent(in)::m
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(regression_result)::res
      res=gaussNewton(residual,x0,m,tolerance,max_iter)
   end function lsqnonlin

   function curvefit(model,xdata,ydata,p0,tolerance,max_iter) result(res)
      procedure(regression_model)::model
      real(dp),intent(in)::xdata(:),ydata(:),p0(:)
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(regression_result)::res
      real(dp),allocatable::p(:),pred(:),r(:),j(:,:),pp(:),pm(:),predp(:),predm(:),dpv(:)
      real(dp)::h,rss,rssn,alpha,tol
      integer::n,m,i,iter,niter,istat
      n=size(p0); m=size(ydata); tol=1.0e-8_dp; if(present(tolerance))tol=tolerance
      niter=500; if(present(max_iter))niter=max_iter
      allocate(p(n),pred(m),r(m),j(m,n),pp(n),pm(n),predp(m),predm(m),dpv(n))
      allocate(res%coefficients(n),res%fitted(m),res%residuals(m))
      if(size(xdata)/=m)then
         res%status=pracma_dimension_mismatch; return
      end if
      p=p0; call model(p,xdata,pred); r=pred-ydata; rss=sum(r*r)
      do iter=1,niter
         do i=1,n
            h=eps_dp**(1.0_dp/3.0_dp)*max(1.0_dp,abs(p(i)))
            pp=p; pm=p; pp(i)=p(i)+h; pm(i)=p(i)-h
            call model(pp,xdata,predp); call model(pm,xdata,predm)
            j(:,i)=(predp-predm)/(2.0_dp*h)
         end do
         dpv=mldivide(j,-r,istat); if(istat/=pracma_ok)exit
         alpha=1.0_dp
         do
            pp=p+alpha*dpv; call model(pp,xdata,predp); rssn=sum((predp-ydata)**2)
            if(rssn<rss .or. alpha<=1.0e-8_dp)exit
            alpha=0.5_dp*alpha
         end do
         p=pp; pred=predp; r=pred-ydata; rss=rssn
         if(sqrt(sum((alpha*dpv)**2))<=tol*(1.0_dp+sqrt(sum(p*p))))exit
      end do
      res%coefficients=p; res%fitted=pred; res%residuals=r; res%rss=rss
      res%iterations=min(iter,niter); res%converged=iter<=niter
      res%status=merge(pracma_ok,pracma_not_converged,res%converged)
   end function curvefit

   function lsqcurvefit(model,xdata,ydata,p0,tolerance,max_iter) result(res)
      procedure(regression_model)::model
      real(dp),intent(in)::xdata(:),ydata(:),p0(:)
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(regression_result)::res
      res=curvefit(model,xdata,ydata,p0,tolerance,max_iter)
   end function lsqcurvefit

   function lsqnonneg(a,b,tolerance,max_iter,status) result(x)
      real(dp),intent(in)::a(:,:),b(:)
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      integer,intent(out),optional::status
      real(dp),allocatable::x(:),g(:),xn(:)
      real(dp)::tol,alpha,lipschitz
      integer::iter,niter,istat
      allocate(x(size(a,2)),g(size(a,2)),xn(size(a,2))); x=0.0_dp
      tol=1.0e-9_dp; if(present(tolerance))tol=tolerance
      niter=10000; if(present(max_iter))niter=max_iter
      lipschitz=max(1.0e-12_dp,maxval(sum(abs(matmul(transpose(a),a)),dim=2)))
      alpha=1.0_dp/lipschitz; istat=pracma_not_converged
      do iter=1,niter
         g=matmul(transpose(a),matmul(a,x)-b)
         xn=max(0.0_dp,x-alpha*g)
         if(sqrt(sum((xn-x)**2))<=tol*(1.0_dp+sqrt(sum(x*x))))then
            x=xn; istat=pracma_ok; exit
         end if
         x=xn
      end do
      if(present(status))status=istat
   end function lsqnonneg

   function L1linreg(a,b,tolerance,max_iter,status) result(x)
      real(dp),intent(in)::a(:,:),b(:)
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      integer,intent(out),optional::status
      real(dp),allocatable::x(:),w(:),aw(:,:),bw(:),r(:)
      real(dp)::tol,delta
      integer::iter,niter,i,istat
      tol=1.0e-8_dp; if(present(tolerance))tol=tolerance
      niter=200; if(present(max_iter))niter=max_iter
      x=mldivide(a,b,istat); allocate(w(size(b)),aw(size(a,1),size(a,2)),bw(size(b)),r(size(b)))
      do iter=1,niter
         r=b-matmul(a,x); w=1.0_dp/max(abs(r),1.0e-8_dp)
         do i=1,size(a,1); aw(i,:)=sqrt(w(i))*a(i,:); bw(i)=sqrt(w(i))*b(i); end do
         r=x; x=mldivide(aw,bw,istat); delta=sqrt(sum((x-r)**2))
         if(delta<=tol*(1.0_dp+sqrt(sum(x*x))))exit
      end do
      if(iter>niter)istat=pracma_not_converged
      if(present(status))status=istat
   end function L1linreg

   function geo_median(points,tolerance,max_iter,status) result(median)
      real(dp),intent(in)::points(:,:)
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      integer,intent(out),optional::status
      real(dp),allocatable::median(:),next(:),dist(:)
      real(dp)::tol,den
      integer::iter,niter,i,istat
      allocate(median(size(points,2)),next(size(points,2)),dist(size(points,1)))
      median=sum(points,dim=1)/real(size(points,1),dp)
      tol=1.0e-9_dp; if(present(tolerance))tol=tolerance
      niter=1000; if(present(max_iter))niter=max_iter
      istat=pracma_not_converged
      do iter=1,niter
         do i=1,size(points,1); dist(i)=sqrt(sum((points(i,:)-median)**2)); end do
         if(minval(dist)<=1.0e-14_dp)then; median=points(minloc(dist,dim=1),:); istat=pracma_ok; exit; end if
         next=0.0_dp; den=0.0_dp
         do i=1,size(points,1); next=next+points(i,:)/dist(i); den=den+1.0_dp/dist(i); end do
         next=next/den
         if(sqrt(sum((next-median)**2))<=tol*(1.0_dp+sqrt(sum(median*median))))then
            median=next; istat=pracma_ok; exit
         end if
         median=next
      end do
      if(present(status))status=istat
   end function geo_median

   function qpspecial(g,x0,max_iter) result(res)
      real(dp),intent(in)::g(:,:)
      real(dp),intent(in),optional::x0(:)
      integer,intent(in),optional::max_iter
      type(qpspecial_result)::res
      real(dp),allocatable::dmat(:,:),dvec(:),amat(:,:),bvec(:),xstart(:)
      type(quadprog_result)::qpres
      integer::n,m,niter,i
      m=size(g,1); n=size(g,2); niter=100
      if(present(max_iter))niter=max_iter
      allocate(dmat(n,n),dvec(n),amat(n,n+1),bvec(n+1),xstart(n))
      dmat=2.0_dp*matmul(transpose(g),g)+1.0e-12_dp*identity_matrix(n)
      dvec=0.0_dp; amat=0.0_dp; bvec=0.0_dp
      amat(:,1)=1.0_dp; bvec(1)=1.0_dp
      do i=1,n; amat(i,i+1)=1.0_dp; end do
      qpres=solve_qp(dmat,dvec,amat,bvec,meq=1)
      allocate(res%x(n),res%d(m)); res%x=qpres%solution; res%d=matmul(g,res%x)
      res%q=sum(res%d*res%d); res%iterations=qpres%iterations(1)
      res%converged=qpres%status==qp_success
      res%status=merge(pracma_ok,pracma_infeasible,res%converged)
      if(present(x0))then
         xstart=max(0.0_dp,x0)
      end if
      if(niter<0)res%status=pracma_invalid_argument
   end function qpspecial

   function qpsolve_projection(d,a,b,meq,tolerance) result(res)
      real(dp),intent(in)::d(:),a(:,:),b(:)
      integer,intent(in),optional::meq
      real(dp),intent(in),optional::tolerance
      type(quadprog_result)::res
      real(dp),allocatable::dmat(:,:)
      real(dp)::tol
      integer::neq
      neq=0; if(present(meq))neq=meq
      tol=sqrt(epsilon(1.0_dp)); if(present(tolerance))tol=max(tolerance,0.0_dp)
      allocate(dmat(size(d),size(d)))
      dmat=identity_matrix(size(d))
      if(tol<0.0_dp) dmat(1,1)=dmat(1,1)
      res=solve_qp(dmat,d,a,b,meq=neq)
   end function qpsolve_projection

   function qpsolve(d,a,b,meq,tolerance) result(res)
      real(dp),intent(in)::d(:),a(:,:),b(:)
      integer,intent(in),optional::meq
      real(dp),intent(in),optional::tolerance
      type(quadprog_result)::res
      res=qpsolve_projection(d,a,b,meq,tolerance)
   end function qpsolve

   function quadprog(c,d,a,b,aeq,beq,lb,ub) result(res)
      real(dp),intent(in)::c(:,:),d(:)
      real(dp),intent(in),optional::a(:,:),b(:),aeq(:,:),beq(:),lb(:),ub(:)
      type(quadprog_result)::res
      real(dp),allocatable::amat(:,:),bvec(:),dvec(:)
      integer::n,q,meq,pos,i
      n=size(d); meq=0; q=0
      if(present(aeq))q=q+size(aeq,1)
      if(present(a))q=q+size(a,1)
      if(present(lb))q=q+n
      if(present(ub))q=q+n
      allocate(amat(n,q),bvec(q),dvec(n)); amat=0.0_dp; bvec=0.0_dp; pos=0
      if(present(aeq))then
         meq=size(aeq,1)
         amat(:,1:meq)=transpose(aeq); bvec(1:meq)=beq; pos=meq
      end if
      if(present(a))then
         amat(:,pos+1:pos+size(a,1))=-transpose(a)
         bvec(pos+1:pos+size(a,1))=-b; pos=pos+size(a,1)
      end if
      if(present(lb))then
         do i=1,n; amat(i,pos+i)=1.0_dp; bvec(pos+i)=lb(i); end do
         pos=pos+n
      end if
      if(present(ub))then
         do i=1,n; amat(i,pos+i)=-1.0_dp; bvec(pos+i)=-ub(i); end do
      end if
      dvec=-d
      res=solve_qp(c,dvec,amat,bvec,meq=meq)
   end function quadprog

   function linprog(cc,a,b,aeq,beq,lb,ub,maximize) result(res)
      real(dp),intent(in)::cc(:)
      real(dp),intent(in),optional::a(:,:),b(:),aeq(:,:),beq(:),lb(:),ub(:)
      logical,intent(in),optional::maximize
      type(linprog_result)::res
      type(quadprog_result)::qpres
      real(dp),allocatable::c(:,:),d(:)
      logical::maxi
      integer::n
      n=size(cc); allocate(c(n,n),d(n)); c=1.0e-10_dp*identity_matrix(n)
      maxi=.false.; if(present(maximize))maxi=maximize
      if(maxi)then; d=-cc; else; d=cc; end if
      qpres=quadprog(c,d,a,b,aeq,beq,lb,ub)
      allocate(res%x(n),res%dual(size(qpres%lagrangian)))
      res%x=qpres%solution; res%dual=qpres%lagrangian; res%value=dot_product(cc,res%x)
      res%iterations=qpres%iterations(1); res%converged=qpres%status==qp_success
      res%status=merge(pracma_ok,pracma_infeasible,res%converged)
   end function linprog

   function fmincon(f,x0,lower,upper,ineq,mineq,eq,meq,tolerance,max_iter) result(res)
      procedure(objective_function)::f
      real(dp),intent(in)::x0(:)
      real(dp),intent(in),optional::lower(:),upper(:)
      procedure(vector_function),optional::ineq,eq
      integer,intent(in),optional::mineq,meq,max_iter
      real(dp),intent(in),optional::tolerance
      type(optimization_result)::res
      real(dp),allocatable::x(:),g(:),xn(:),ci(:),ce(:)
      real(dp)::fx,fn,penalty,alpha,tol,rho
      integer::iter,niter,mi,me
      mi=0; me=0; if(present(mineq))mi=mineq; if(present(meq))me=meq
      tol=1.0e-7_dp; if(present(tolerance))tol=tolerance
      niter=5000; if(present(max_iter))niter=max_iter
      allocate(x(size(x0)),g(size(x0)),xn(size(x0)),ci(mi),ce(me),res%x(size(x0)),res%history(niter))
      x=x0; call apply_bounds(x,lower,upper); rho=100.0_dp
      fx=penalized_value(f,x,ineq,mi,eq,me,rho,ci,ce)
      do iter=1,niter
         call penalty_gradient(f,x,g,ineq,mi,eq,me,rho,lower,upper)
         if(sqrt(sum(g*g))<=tol)exit
         alpha=1.0_dp
         do
            xn=x-alpha*g; call apply_bounds(xn,lower,upper)
            fn=penalized_value(f,xn,ineq,mi,eq,me,rho,ci,ce)
            if(fn<fx .or. alpha<=1.0e-12_dp)exit
            alpha=0.5_dp*alpha
         end do
         x=xn; fx=fn; res%history(iter)=fx
         penalty=constraint_violation(ineq,mi,eq,me,x,ci,ce)
         if(penalty>sqrt(tol) .and. modulo(iter,100)==0)rho=min(1.0e10_dp,10.0_dp*rho)
         if(alpha<=1.0e-12_dp)exit
      end do
      res%x=x; res%value=f(x); res%iterations=min(iter,niter)
      penalty=constraint_violation(ineq,mi,eq,me,x,ci,ce)
      res%converged=penalty<=sqrt(tol) .and. sqrt(sum(g*g))<=10.0_dp*sqrt(tol)
      res%status=merge(pracma_ok,pracma_not_converged,res%converged)
      call shrink_history(res)
   end function fmincon

   function penalized_value(f,x,ineq,mi,eq,me,rho,ci,ce) result(v)
      procedure(objective_function)::f
      real(dp),intent(in)::x(:),rho
      procedure(vector_function),optional::ineq,eq
      integer,intent(in)::mi,me
      real(dp),intent(inout)::ci(:),ce(:)
      real(dp)::v
      v=f(x)
      if(present(ineq).and.mi>0)then; call ineq(x,ci); v=v+rho*sum(max(ci,0.0_dp)**2); end if
      if(present(eq).and.me>0)then; call eq(x,ce); v=v+rho*sum(ce*ce); end if
   end function penalized_value

   subroutine penalty_gradient(f,x,g,ineq,mi,eq,me,rho,lower,upper)
      procedure(objective_function)::f
      real(dp),intent(in)::x(:),rho
      real(dp),intent(out)::g(:)
      procedure(vector_function),optional::ineq,eq
      integer,intent(in)::mi,me
      real(dp),intent(in),optional::lower(:),upper(:)
      real(dp),allocatable::xp(:),xm(:),ci(:),ce(:)
      real(dp)::h,fp,fm
      integer::i
      allocate(xp(size(x)),xm(size(x)),ci(mi),ce(me))
      do i=1,size(x)
         h=eps_dp**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x(i)))
         xp=x; xm=x; xp(i)=xp(i)+h; xm(i)=xm(i)-h
         call apply_bounds(xp,lower,upper); call apply_bounds(xm,lower,upper)
         fp=penalized_value(f,xp,ineq,mi,eq,me,rho,ci,ce)
         fm=penalized_value(f,xm,ineq,mi,eq,me,rho,ci,ce)
         g(i)=(fp-fm)/max(tiny(1.0_dp),xp(i)-xm(i))
      end do
   end subroutine penalty_gradient

   function constraint_violation(ineq,mi,eq,me,x,ci,ce) result(v)
      procedure(vector_function),optional::ineq,eq
      integer,intent(in)::mi,me
      real(dp),intent(in)::x(:)
      real(dp),intent(inout)::ci(:),ce(:)
      real(dp)::v
      v=0.0_dp
      if(present(ineq).and.mi>0)then; call ineq(x,ci); v=max(v,maxval(max(ci,0.0_dp))); end if
      if(present(eq).and.me>0)then; call eq(x,ce); v=max(v,maxval(abs(ce))); end if
   end function constraint_violation

   subroutine apply_bounds(x,lower,upper)
      real(dp),intent(inout)::x(:)
      real(dp),intent(in),optional::lower(:),upper(:)
      if(present(lower))x=max(x,lower)
      if(present(upper))x=min(x,upper)
   end subroutine apply_bounds

   subroutine shrink_history(res)
      type(optimization_result),intent(inout)::res
      real(dp),allocatable::tmp(:)
      integer::n
      if(.not.allocated(res%history))return
      n=max(0,res%iterations)
      allocate(tmp(n))
      if(n>0)tmp=res%history(1:n)
      call move_alloc(tmp,res%history)
   end subroutine shrink_history

end module pracma_optimization
