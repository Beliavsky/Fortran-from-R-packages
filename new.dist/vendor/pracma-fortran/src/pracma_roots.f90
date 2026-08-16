! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_roots
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use pracma_kinds, only : dp, eps_dp
   use pracma_status
   use pracma_types, only : root_result
   use pracma_callbacks
   use pracma_linalg, only : mldivide, solve_linear, identity_matrix
   use pracma_differentiation, only : numderiv, jacobian
   use pracma_polynomial, only : horner, hornerdefl
   implicit none
   private

   public :: bisect, secant, regulaFalsi, brentDekker, fzero
   public :: newtonRaphson, halley, newtonHorner, muller, ridders
   public :: findzeros, findmins, fsolve, fzsolve, newtonsys, broyden
   public :: itersolve, aitken

contains

   function bisect(f,a,b,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(root_result)::res
      real(dp)::left,right,mid,fl,fr,fm,tol
      integer::niter,iter
      tol=1.0e-12_dp; if(present(tolerance))tol=tolerance
      niter=200; if(present(max_iter))niter=max_iter
      left=a; right=b; fl=f(left); fr=f(right); res%evaluations=2
      if(.not.ieee_is_finite(fl) .or. .not.ieee_is_finite(fr))then
         res%status=pracma_nonfinite; return
      end if
      if(abs(fl)<=tol)then
         res%root=left; res%value=fl; res%converged=.true.; return
      else if(abs(fr)<=tol)then
         res%root=right; res%value=fr; res%converged=.true.; return
      else if(fl*fr>0.0_dp)then
         res%status=pracma_not_bracketed; return
      end if
      do iter=1,niter
         mid=left+0.5_dp*(right-left); fm=f(mid); res%evaluations=res%evaluations+1
         if(abs(fm)<=tol .or. abs(right-left)<=tol*(1.0_dp+abs(mid)))then
            res%root=mid; res%value=fm; res%iterations=iter
            res%status=pracma_ok; res%converged=.true.; return
         end if
         if(fl*fm<=0.0_dp)then
            right=mid; fr=fm
         else
            left=mid; fl=fm
         end if
      end do
      res%root=mid; res%value=fm; res%iterations=niter
      res%status=pracma_not_converged
   end function bisect

   function secant(f,x0,x1,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::x0,x1
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(root_result)::res
      real(dp)::a,b,fa,fb,c,fc,den,tol
      integer::iter,niter
      tol=1.0e-12_dp; if(present(tolerance))tol=tolerance
      niter=100; if(present(max_iter))niter=max_iter
      a=x0; b=x1; fa=f(a); fb=f(b); res%evaluations=2
      do iter=1,niter
         den=fb-fa
         if(abs(den)<=eps_dp*max(1.0_dp,abs(fa),abs(fb)))then
            res%status=pracma_singular; exit
         end if
         c=b-fb*(b-a)/den; fc=f(c); res%evaluations=res%evaluations+1
         if(abs(fc)<=tol .or. abs(c-b)<=tol*(1.0_dp+abs(c)))then
            res%root=c; res%value=fc; res%iterations=iter
            res%status=pracma_ok; res%converged=.true.; return
         end if
         a=b; fa=fb; b=c; fb=fc
      end do
      res%root=b; res%value=fb; res%iterations=min(iter,niter)
      if(res%status==pracma_ok)res%status=pracma_not_converged
   end function secant

   function regulaFalsi(f,a,b,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(root_result)::res
      real(dp)::left,right,fl,fr,c,fc,tol,old
      integer::iter,niter
      tol=1.0e-12_dp; if(present(tolerance))tol=tolerance
      niter=200; if(present(max_iter))niter=max_iter
      left=a; right=b; fl=f(left); fr=f(right); res%evaluations=2; old=huge(1.0_dp)
      if(fl*fr>0.0_dp)then
         res%status=pracma_not_bracketed; return
      end if
      do iter=1,niter
         c=(left*fr-right*fl)/(fr-fl); fc=f(c); res%evaluations=res%evaluations+1
         if(abs(fc)<=tol .or. abs(c-old)<=tol*(1.0_dp+abs(c)))then
            res%root=c; res%value=fc; res%iterations=iter
            res%status=pracma_ok; res%converged=.true.; return
         end if
         old=c
         if(fl*fc<0.0_dp)then
            right=c; fr=fc; fl=0.5_dp*fl
         else
            left=c; fl=fc; fr=0.5_dp*fr
         end if
      end do
      res%root=c; res%value=fc; res%iterations=niter; res%status=pracma_not_converged
   end function regulaFalsi

   function brentDekker(f,a,b,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(root_result)::res
      real(dp)::aa,bb,cc,dd,ee,fa,fb,fc,p,q,r,s,tol,tol1,xm,temp
      integer::iter,niter
      tol=1.0e-12_dp; if(present(tolerance))tol=tolerance
      niter=200; if(present(max_iter))niter=max_iter
      aa=a; bb=b; fa=f(aa); fb=f(bb); res%evaluations=2
      if(fa*fb>0.0_dp)then
         res%status=pracma_not_bracketed; return
      end if
      cc=bb; fc=fb; dd=bb-aa; ee=dd
      do iter=1,niter
         if(fb*fc>0.0_dp)then
            cc=aa; fc=fa; dd=bb-aa; ee=dd
         end if
         if(abs(fc)<abs(fb))then
            temp=aa; aa=bb; bb=cc; cc=temp
            temp=fa; fa=fb; fb=fc; fc=temp
         end if
         tol1=2.0_dp*eps_dp*abs(bb)+0.5_dp*tol
         xm=0.5_dp*(cc-bb)
         if(abs(xm)<=tol1 .or. abs(fb)<=tol)then
            res%root=bb; res%value=fb; res%iterations=iter
            res%status=pracma_ok; res%converged=.true.; return
         end if
         if(abs(ee)>=tol1 .and. abs(fa)>abs(fb))then
            s=fb/fa
            if(abs(aa-cc)<=tol1)then
               p=2.0_dp*xm*s; q=1.0_dp-s
            else
               q=fa/fc; r=fb/fc
               p=s*(2.0_dp*xm*q*(q-r)-(bb-aa)*(r-1.0_dp))
               q=(q-1.0_dp)*(r-1.0_dp)*(s-1.0_dp)
            end if
            if(p>0.0_dp)q=-q
            p=abs(p)
            if(2.0_dp*p<min(3.0_dp*xm*q-abs(tol1*q),abs(ee*q)))then
               ee=dd; dd=p/q
            else
               dd=xm; ee=dd
            end if
         else
            dd=xm; ee=dd
         end if
         aa=bb; fa=fb
         if(abs(dd)>tol1)then
            bb=bb+dd
         else
            bb=bb+sign(tol1,xm)
         end if
         fb=f(bb); res%evaluations=res%evaluations+1
      end do
      res%root=bb; res%value=fb; res%iterations=niter; res%status=pracma_not_converged
   end function brentDekker

   function fzero(f,a,b,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a
      real(dp),intent(in),optional::b,tolerance
      integer,intent(in),optional::max_iter
      type(root_result)::res
      real(dp)::left,right,step
      integer::k
      if(present(b))then
         res=brentDekker(f,a,b,tolerance,max_iter)
      else
         step=max(0.01_dp,0.01_dp*abs(a)); left=a-step; right=a+step
         do k=1,60
            if(f(left)*f(right)<=0.0_dp)exit
            step=step*1.6_dp; left=a-step; right=a+step
         end do
         if(k>60)then
            res%status=pracma_not_bracketed
         else
            res=brentDekker(f,left,right,tolerance,max_iter)
         end if
      end if
   end function fzero

   function newtonRaphson(f,df,x0,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      procedure(scalar_derivative),optional::df
      real(dp),intent(in)::x0
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(root_result)::res
      real(dp)::x,fx,der,dx,tol
      integer::iter,niter
      tol=1.0e-12_dp; if(present(tolerance))tol=tolerance
      niter=100; if(present(max_iter))niter=max_iter
      x=x0
      do iter=1,niter
         fx=f(x); res%evaluations=res%evaluations+1
         if(abs(fx)<=tol)then
            res%root=x; res%value=fx; res%iterations=iter-1
            res%status=pracma_ok; res%converged=.true.; return
         end if
         if(present(df))then
            der=df(x); res%evaluations=res%evaluations+1
         else
            der=numderiv(f,x); res%evaluations=res%evaluations+2
         end if
         if(abs(der)<=eps_dp*max(1.0_dp,abs(fx)))then
            res%status=pracma_singular; exit
         end if
         dx=fx/der; x=x-dx
         if(abs(dx)<=tol*(1.0_dp+abs(x)))then
            fx=f(x); res%evaluations=res%evaluations+1
            res%root=x; res%value=fx; res%iterations=iter
            res%status=pracma_ok; res%converged=.true.; return
         end if
      end do
      res%root=x; res%value=f(x); res%iterations=min(iter,niter)
      if(res%status==pracma_ok)res%status=pracma_not_converged
   end function newtonRaphson

   function halley(f,df,ddf,x0,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      procedure(scalar_derivative)::df,ddf
      real(dp),intent(in)::x0
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(root_result)::res
      real(dp)::x,fx,d1,d2,den,dx,tol
      integer::iter,niter
      tol=1.0e-12_dp; if(present(tolerance))tol=tolerance
      niter=100; if(present(max_iter))niter=max_iter
      x=x0
      do iter=1,niter
         fx=f(x); d1=df(x); d2=ddf(x); res%evaluations=res%evaluations+3
         den=2.0_dp*d1*d1-fx*d2
         if(abs(den)<=eps_dp*max(1.0_dp,abs(d1*d1)))then
            res%status=pracma_singular; exit
         end if
         dx=2.0_dp*fx*d1/den; x=x-dx
         if(abs(dx)<=tol*(1.0_dp+abs(x)))then
            res%root=x; res%value=f(x); res%iterations=iter
            res%status=pracma_ok; res%converged=.true.; return
         end if
      end do
      res%root=x; res%value=f(x); res%iterations=min(iter,niter)
      if(res%status==pracma_ok)res%status=pracma_not_converged
   end function halley

   function newtonHorner(p,x0,tolerance,max_iter) result(res)
      real(dp),intent(in)::p(:),x0
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(root_result)::res
      real(dp)::x,fx,dfx,dx,tol
      integer::iter,niter
      tol=1.0e-12_dp; if(present(tolerance))tol=tolerance
      niter=100; if(present(max_iter))niter=max_iter
      x=x0
      do iter=1,niter
         call horner(p,x,fx,dfx)
         if(abs(dfx)<=eps_dp*max(1.0_dp,abs(fx)))then
            res%status=pracma_singular; exit
         end if
         dx=fx/dfx; x=x-dx
         if(abs(dx)<=tol*(1.0_dp+abs(x)))then
            call horner(p,x,fx)
            res%root=x; res%value=fx; res%iterations=iter
            res%status=pracma_ok; res%converged=.true.; return
         end if
      end do
      call horner(p,x,fx); res%root=x; res%value=fx; res%iterations=min(iter,niter)
      if(res%status==pracma_ok)res%status=pracma_not_converged
   end function newtonHorner

   function muller(f,x0,x1,x2,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::x0,x1,x2
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(root_result)::res
      real(dp)::a,b,c,h1,h2,d1,d2,d,disc,e,dx,x3,f0,f1,f2,f3,tol
      integer::iter,niter
      tol=1.0e-12_dp; if(present(tolerance))tol=tolerance
      niter=100; if(present(max_iter))niter=max_iter
      a=x0; b=x1; c=x2; f0=f(a); f1=f(b); f2=f(c); res%evaluations=3
      do iter=1,niter
         h1=b-a; h2=c-b; d1=(f1-f0)/h1; d2=(f2-f1)/h2
         d=(d2-d1)/(h2+h1); b=d2+h2*d; disc=b*b-4.0_dp*f2*d
         if(disc<0.0_dp)then
            res%status=pracma_unsupported; exit
         end if
         disc=sqrt(disc)
         if(abs(b+disc)>abs(b-disc))then; e=b+disc; else; e=b-disc; end if
         if(abs(e)<=tiny(1.0_dp))then; res%status=pracma_singular; exit; end if
         dx=-2.0_dp*f2/e; x3=c+dx; f3=f(x3); res%evaluations=res%evaluations+1
         if(abs(dx)<=tol*(1.0_dp+abs(x3)) .or. abs(f3)<=tol)then
            res%root=x3; res%value=f3; res%iterations=iter
            res%status=pracma_ok; res%converged=.true.; return
         end if
         a=b; f0=f1; b=c; f1=f2; c=x3; f2=f3
      end do
      res%root=c; res%value=f2; res%iterations=min(iter,niter)
      if(res%status==pracma_ok)res%status=pracma_not_converged
   end function muller

   function ridders(f,a,b,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(root_result)::res
      real(dp)::x1,x2,f1,f2,xm,fm,s,xnew,fnew,tol
      integer::iter,niter
      tol=1.0e-12_dp; if(present(tolerance))tol=tolerance
      niter=100; if(present(max_iter))niter=max_iter
      x1=a; x2=b; f1=f(x1); f2=f(x2); res%evaluations=2
      if(f1*f2>0.0_dp)then; res%status=pracma_not_bracketed; return; end if
      do iter=1,niter
         xm=0.5_dp*(x1+x2); fm=f(xm); res%evaluations=res%evaluations+1
         s=sqrt(max(0.0_dp,fm*fm-f1*f2))
         if(s<=tiny(1.0_dp))then; res%status=pracma_singular; exit; end if
         xnew=xm+(xm-x1)*sign(1.0_dp,f1-f2)*fm/s
         fnew=f(xnew); res%evaluations=res%evaluations+1
         if(abs(fnew)<=tol .or. abs(xnew-xm)<=tol*(1.0_dp+abs(xnew)))then
            res%root=xnew; res%value=fnew; res%iterations=iter
            res%status=pracma_ok; res%converged=.true.; return
         end if
         if(fm*fnew<0.0_dp)then
            x1=xm; f1=fm; x2=xnew; f2=fnew
         else if(f1*fnew<0.0_dp)then
            x2=xnew; f2=fnew
         else
            x1=xnew; f1=fnew
         end if
      end do
      res%root=xnew; res%value=fnew; res%iterations=min(iter,niter)
      if(res%status==pracma_ok)res%status=pracma_not_converged
   end function ridders

   function findzeros(f,a,b,n,tolerance) result(z)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      integer,intent(in),optional::n
      real(dp),intent(in),optional::tolerance
      real(dp),allocatable::z(:),tmp(:)
      real(dp)::x1,x2,f1,f2
      type(root_result)::r
      integer::m,i,k
      m=1000; if(present(n))m=n
      allocate(tmp(m)); k=0; x1=a; f1=f(x1)
      do i=1,m
         x2=a+real(i,dp)*(b-a)/real(m,dp); f2=f(x2)
         if(abs(f1)<=merge(tolerance,1.0e-10_dp,present(tolerance)))then
            if(k==0 .or. abs(x1-tmp(k))>1.0e-8_dp)then; k=k+1; tmp(k)=x1; end if
         else if(f1*f2<0.0_dp)then
            r=brentDekker(f,x1,x2,tolerance)
            if(r%converged)then
               if(k==0 .or. abs(r%root-tmp(k))>1.0e-8_dp)then; k=k+1; tmp(k)=r%root; end if
            end if
         end if
         x1=x2; f1=f2
      end do
      allocate(z(k)); if(k>0)z=tmp(1:k)
   end function findzeros

   function findmins(f,a,b,n) result(xmin)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b
      integer,intent(in),optional::n
      real(dp),allocatable::xmin(:),tmp(:)
      real(dp)::x0,x1,x2,f0,f1,f2
      integer::m,i,k
      m=1000; if(present(n))m=n
      allocate(tmp(m)); k=0
      x0=a; f0=f(x0); x1=a+(b-a)/real(m,dp); f1=f(x1)
      do i=2,m
         x2=a+real(i,dp)*(b-a)/real(m,dp); f2=f(x2)
         if(f1<f0 .and. f1<=f2)then; k=k+1; tmp(k)=x1; end if
         x0=x1; f0=f1; x1=x2; f1=f2
      end do
      allocate(xmin(k)); if(k>0)xmin=tmp(1:k)
   end function findmins

   subroutine newtonsys(f,x,status,tolerance,max_iter)
      procedure(vector_function)::f
      real(dp),intent(inout)::x(:)
      integer,intent(out),optional::status
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      real(dp),allocatable::fx(:),j(:,:),dx(:)
      real(dp)::tol
      integer::n,iter,niter,istat
      n=size(x); allocate(fx(n),j(n,n),dx(n))
      tol=1.0e-10_dp; if(present(tolerance))tol=tolerance
      niter=100; if(present(max_iter))niter=max_iter
      istat=pracma_not_converged
      do iter=1,niter
         call f(x,fx)
         if(sqrt(sum(fx*fx))<=tol)then; istat=pracma_ok; exit; end if
         call jacobian(f,x,n,j,status=istat)
         dx=mldivide(j,-fx,istat)
         if(istat/=pracma_ok)exit
         x=x+dx
         if(sqrt(sum(dx*dx))<=tol*(1.0_dp+sqrt(sum(x*x))))then
            call f(x,fx); istat=merge(pracma_ok,pracma_not_converged,sqrt(sum(fx*fx))<=sqrt(tol))
            exit
         end if
      end do
      if(present(status))status=istat
   end subroutine newtonsys

   function fsolve(f,x0,tolerance,max_iter,status) result(x)
      procedure(vector_function)::f
      real(dp),intent(in)::x0(:)
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      integer,intent(out),optional::status
      real(dp),allocatable::x(:)
      integer::istat
      allocate(x(size(x0))); x=x0
      call newtonsys(f,x,istat,tolerance,max_iter)
      if(present(status))status=istat
   end function fsolve

   function fzsolve(f,x0,tolerance,max_iter,status) result(x)
      procedure(vector_function)::f
      real(dp),intent(in)::x0(:)
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      integer,intent(out),optional::status
      real(dp),allocatable::x(:)
      x=fsolve(f,x0,tolerance,max_iter,status)
   end function fzsolve

   subroutine broyden(f,x,status,tolerance,max_iter)
      procedure(vector_function)::f
      real(dp),intent(inout)::x(:)
      integer,intent(out),optional::status
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      real(dp),allocatable::fx(:),fnew(:),b(:,:),s(:),y(:),xnew(:)
      real(dp)::tol,den,alpha
      integer::n,iter,niter,istat
      n=size(x); allocate(fx(n),fnew(n),b(n,n),s(n),y(n),xnew(n))
      tol=1.0e-10_dp; if(present(tolerance))tol=tolerance
      niter=200; if(present(max_iter))niter=max_iter
      b=identity_matrix(n); call f(x,fx); istat=pracma_not_converged
      do iter=1,niter
         s=mldivide(b,-fx,istat)
         if(istat/=pracma_ok)exit
         alpha=1.0_dp
         do
            xnew=x+alpha*s; call f(xnew,fnew)
            if(sqrt(sum(fnew*fnew))<sqrt(sum(fx*fx)) .or. alpha<1.0e-6_dp)exit
            alpha=0.5_dp*alpha
         end do
         s=xnew-x; y=fnew-fx; den=dot_product(s,s)
         if(den>tiny(1.0_dp))b=b+spread(y-matmul(b,s),2,n)*spread(s,1,n)/den
         x=xnew; fx=fnew
         if(sqrt(sum(fx*fx))<=tol)then; istat=pracma_ok; exit; end if
      end do
      if(present(status))status=istat
   end subroutine broyden

   function itersolve(g,x0,tolerance,max_iter,status) result(x)
      procedure(scalar_function)::g
      real(dp),intent(in)::x0
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      integer,intent(out),optional::status
      real(dp)::x,xnew,tol
      integer::iter,niter,istat
      tol=1.0e-12_dp; if(present(tolerance))tol=tolerance
      niter=1000; if(present(max_iter))niter=max_iter
      x=x0; istat=pracma_not_converged
      do iter=1,niter
         xnew=g(x)
         if(abs(xnew-x)<=tol*(1.0_dp+abs(xnew)))then; x=xnew; istat=pracma_ok; exit; end if
         x=xnew
      end do
      if(present(status))status=istat
   end function itersolve

   function aitken(sequence) result(accelerated)
      real(dp),intent(in)::sequence(:)
      real(dp),allocatable::accelerated(:)
      real(dp)::den
      integer::i,n
      n=max(0,size(sequence)-2); allocate(accelerated(n))
      do i=1,n
         den=sequence(i+2)-2.0_dp*sequence(i+1)+sequence(i)
         if(abs(den)<=eps_dp*max(1.0_dp,maxval(abs(sequence(i:i+2)))))then
            accelerated(i)=sequence(i+2)
         else
            accelerated(i)=sequence(i)-(sequence(i+1)-sequence(i))**2/den
         end if
      end do
   end function aitken

end module pracma_roots
