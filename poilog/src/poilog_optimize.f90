! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
module poilog_optimize
   use poilog_kinds, only : dp
   use poilog_math, only : is_finite_dp
   implicit none
   private
   public :: optim_result, minimize

   type :: optim_result
      real(dp), allocatable :: par(:)
      real(dp) :: value = huge(1.0_dp)
      integer :: convergence = 1
      integer :: iterations = 0
   end type optim_result

   abstract interface
      function objective_fn(x) result(f)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp) :: f
      end function objective_fn
   end interface

contains

   function minimize(f,x0,method,maxit,tol) result(res)
      procedure(objective_fn) :: f
      real(dp), intent(in) :: x0(:)
      character(len=*), intent(in), optional :: method
      integer, intent(in), optional :: maxit
      real(dp), intent(in), optional :: tol
      type(optim_result) :: res
      character(len=:), allocatable :: m
      integer :: mi
      real(dp) :: tt
      mi=1000; if(present(maxit)) mi=maxit
      tt=1.0e-6_dp; if(present(tol)) tt=tol
      m='bfgs'; if(present(method)) m=lower_string(trim(method))
      if(index(m,'nelder')>0) then
         res=nelder_mead(f,x0,mi,tt)
      else
         res=bfgs(f,x0,mi,tt)
      end if
   end function minimize

   function bfgs(f,x0,maxit,tol) result(res)
      procedure(objective_fn) :: f
      real(dp), intent(in) :: x0(:),tol
      integer, intent(in) :: maxit
      type(optim_result) :: res
      integer :: n,i,it,ls
      real(dp), allocatable :: x(:),xn(:),g(:),gn(:),p(:),h(:,:),s(:),y(:),hy(:)
      real(dp) :: fx,fn,alpha,gp,ys,yhy,gnorm

      n=size(x0)
      allocate(x(n),xn(n),g(n),gn(n),p(n),h(n,n),s(n),y(n),hy(n))
      x=x0; h=0.0_dp
      do i=1,n; h(i,i)=1.0_dp; end do
      fx=f(x); call gradient(f,x,g)
      if(.not.is_finite_dp(fx)) then
         allocate(res%par(n)); res%par=x; res%value=fx; res%convergence=2; return
      end if

      do it=1,maxit
         gnorm=maxval(abs(g))
         if(gnorm <= tol*(1.0_dp+abs(fx))) then
            res%convergence=0; exit
         end if
         p=-matmul(h,g)
         gp=dot_product(g,p)
         if(gp >= -epsilon(1.0_dp)) then
            p=-g; gp=-dot_product(g,g); h=0.0_dp
            do i=1,n; h(i,i)=1.0_dp; end do
         end if
         alpha=1.0_dp
         fn=huge(1.0_dp)
         do ls=1,50
            xn=x+alpha*p
            fn=f(xn)
            if(is_finite_dp(fn)) then
               if(fn <= fx+1.0e-4_dp*alpha*gp) exit
            end if
            alpha=0.5_dp*alpha
         end do
         if(ls>50 .or. alpha < 1.0e-12_dp) then
            res%convergence=3; exit
         end if
         call gradient(f,xn,gn)
         s=xn-x; y=gn-g; ys=dot_product(y,s)
         if(ys > sqrt(epsilon(1.0_dp))*sqrt(dot_product(s,s)*dot_product(y,y))) then
            hy=matmul(h,y); yhy=dot_product(y,hy)
            h=h+((ys+yhy)/(ys*ys))*outer(s,s)-(outer(hy,s)+outer(s,hy))/ys
         else
            h=0.0_dp; do i=1,n; h(i,i)=1.0_dp; end do
         end if
         x=xn; g=gn; fx=fn
      end do
      if(it>maxit) res%convergence=1
      allocate(res%par(n)); res%par=x; res%value=fx; res%iterations=min(it,maxit)
   end function bfgs

   function nelder_mead(f,x0,maxit,tol) result(res)
      procedure(objective_fn) :: f
      real(dp), intent(in) :: x0(:),tol
      integer, intent(in) :: maxit
      type(optim_result) :: res
      integer :: n,i,it,ilo,ihi,inhi
      real(dp), allocatable :: simp(:,:),fv(:),cent(:),xr(:),xe(:),xc(:)
      real(dp) :: spread
      n=size(x0)
      allocate(simp(n,n+1),fv(n+1),cent(n),xr(n),xe(n),xc(n))
      simp(:,1)=x0
      do i=1,n
         simp(:,i+1)=x0
         simp(i,i+1)=x0(i)+0.05_dp*max(1.0_dp,abs(x0(i)))
      end do
      do i=1,n+1; fv(i)=f(simp(:,i)); end do
      do it=1,maxit
         call order_simplex(fv,ilo,ihi,inhi)
         spread=maxval(abs(fv-fv(ilo)))
         if(spread <= tol*(1.0_dp+abs(fv(ilo)))) then
            res%convergence=0; exit
         end if
         cent=0.0_dp
         do i=1,n+1; if(i/=ihi) cent=cent+simp(:,i); end do
         cent=cent/real(n,dp)
         xr=cent+(cent-simp(:,ihi))
         if(f(xr)<fv(ilo)) then
            xe=cent+2.0_dp*(xr-cent)
            if(f(xe)<f(xr)) then; simp(:,ihi)=xe; else; simp(:,ihi)=xr; end if
         else if(f(xr)<fv(inhi)) then
            simp(:,ihi)=xr
         else
            xc=cent+0.5_dp*(simp(:,ihi)-cent)
            if(f(xc)<fv(ihi)) then
               simp(:,ihi)=xc
            else
               do i=1,n+1
                  if(i/=ilo) simp(:,i)=simp(:,ilo)+0.5_dp*(simp(:,i)-simp(:,ilo))
               end do
            end if
         end if
         do i=1,n+1; fv(i)=f(simp(:,i)); end do
      end do
      call order_simplex(fv,ilo,ihi,inhi)
      allocate(res%par(n)); res%par=simp(:,ilo); res%value=fv(ilo); res%iterations=min(it,maxit)
      if(it>maxit) res%convergence=1
   end function nelder_mead

   subroutine gradient(f,x,g)
      procedure(objective_fn) :: f
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      real(dp), allocatable :: xp(:),xm(:)
      real(dp) :: h,fp,fm
      integer :: i
      allocate(xp(size(x)),xm(size(x)))
      do i=1,size(x)
         h=epsilon(1.0_dp)**(1.0_dp/3.0_dp)*(1.0_dp+abs(x(i)))
         xp=x; xm=x; xp(i)=xp(i)+h; xm(i)=xm(i)-h
         fp=f(xp); fm=f(xm)
         if(is_finite_dp(fp).and.is_finite_dp(fm)) then
            g(i)=(fp-fm)/(2.0_dp*h)
         else
            g(i)=0.0_dp
         end if
      end do
   end subroutine gradient

   pure function outer(a,b) result(c)
      real(dp), intent(in) :: a(:),b(:)
      real(dp) :: c(size(a),size(b))
      integer :: i,j
      do j=1,size(b); do i=1,size(a); c(i,j)=a(i)*b(j); end do; end do
   end function outer

   subroutine order_simplex(fv,ilo,ihi,inhi)
      real(dp), intent(in) :: fv(:)
      integer, intent(out) :: ilo,ihi,inhi
      integer :: i
      ilo=minloc(fv,dim=1); ihi=maxloc(fv,dim=1); inhi=ilo
      do i=1,size(fv)
         if(i/=ihi) then
            if(inhi==ilo .or. fv(i)>fv(inhi)) inhi=i
         end if
      end do
   end subroutine order_simplex

   pure function lower_string(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      integer :: i,k
      t=s
      do i=1,len(s)
         k=iachar(t(i:i)); if(k>=iachar('A').and.k<=iachar('Z')) t(i:i)=achar(k+32)
      end do
   end function lower_string

end module poilog_optimize
