! SPDX-License-Identifier: GPL-2.0-only
module streg_optimize
   use streg_kinds, only : dp
   use streg_linalg, only : identity_matrix
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: bfgs_minimize, nelder_mead, numerical_hessian

   abstract interface
      function objective_function(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_function
   end interface

contains

   subroutine numerical_gradient(objective,x,g)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(size(x))
      real(dp) :: h,fp,fm
      real(dp) :: xp(size(x)),xm(size(x))
      integer :: i
      do i=1,size(x)
         h=epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x(i)))
         xp=x; xm=x; xp(i)=xp(i)+h; xm(i)=xm(i)-h
         fp=objective(xp); fm=objective(xm)
         if(ieee_is_finite(fp).and.ieee_is_finite(fm))then
            g(i)=(fp-fm)/(2.0_dp*h)
         else
            g(i)=0.0_dp
         end if
      end do
   end subroutine numerical_gradient

   subroutine bfgs_minimize(objective,start,optimum,value,converged,iterations,max_iter,tolerance)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: start(:)
      real(dp), allocatable, intent(out) :: optimum(:)
      real(dp), intent(out) :: value
      logical, intent(out) :: converged
      integer, intent(out) :: iterations
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: x(:),g(:),gnew(:),hmat(:,:),p(:),xnew(:),s(:),y(:),vmat(:,:)
      real(dp) :: f,fnew,alpha,tol,ys,gd,rho,rel
      integer :: n,maxit,ls
      n=size(start); maxit=1000; if(present(max_iter))maxit=max(1,max_iter)
      tol=1.0e-8_dp; if(present(tolerance))tol=max(tolerance,epsilon(1.0_dp))
      allocate(x(n),g(n),gnew(n),hmat(n,n),p(n),xnew(n),s(n),y(n),vmat(n,n),optimum(n))
      x=start; f=objective(x); hmat=identity_matrix(n)
      call numerical_gradient(objective,x,g)
      converged=.false.; iterations=0
      do iterations=1,maxit
         if(maxval(abs(g))<=sqrt(tol)*(1.0_dp+abs(f)))then
            converged=.true.; exit
         end if
         p=-matmul(hmat,g); gd=dot_product(g,p)
         if(gd>=-epsilon(1.0_dp))then
            p=-g; gd=-dot_product(g,g); hmat=identity_matrix(n)
         end if
         alpha=1.0_dp; fnew=huge(1.0_dp)
         do ls=1,60
            xnew=x+alpha*p; fnew=objective(xnew)
            if(ieee_is_finite(fnew))then
               if(fnew<=f+1.0e-4_dp*alpha*gd)exit
            end if
            alpha=0.5_dp*alpha
         end do
         if(alpha<1.0e-16_dp .or. .not.ieee_is_finite(fnew))exit
         call numerical_gradient(objective,xnew,gnew)
         s=xnew-x; y=gnew-g; ys=dot_product(y,s)
         rel=abs(fnew-f)/(1.0_dp+abs(f))
         if(ys>sqrt(epsilon(1.0_dp))*sqrt(max(dot_product(s,s)*dot_product(y,y),tiny(1.0_dp))))then
            rho=1.0_dp/ys
            vmat=identity_matrix(n)-rho*spread(s,2,n)*spread(y,1,n)
            hmat=matmul(vmat,matmul(hmat,transpose(vmat)))+rho*spread(s,2,n)*spread(s,1,n)
            hmat=0.5_dp*(hmat+transpose(hmat))
         else
            hmat=identity_matrix(n)
         end if
         x=xnew; f=fnew; g=gnew
         if(rel<=tol .and. maxval(abs(s))<=sqrt(tol)*(1.0_dp+maxval(abs(x))))then
            converged=.true.; exit
         end if
      end do
      optimum=x; value=f
   end subroutine bfgs_minimize

   subroutine nelder_mead(objective,start,optimum,value,converged,iterations,max_iter,tolerance)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: start(:)
      real(dp), allocatable, intent(out) :: optimum(:)
      real(dp), intent(out) :: value
      logical, intent(out) :: converged
      integer, intent(out) :: iterations
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: simplex(:,:),f(:),centroid(:),xr(:),xe(:),xc(:)
      real(dp) :: fr,fe,fc,tol,scale
      integer :: n,i,j,best,worst,second,maxit
      n=size(start); maxit=1500; if(present(max_iter))maxit=max(1,max_iter)
      tol=1.0e-8_dp; if(present(tolerance))tol=max(tolerance,epsilon(1.0_dp))
      allocate(simplex(n,n+1),f(n+1),centroid(n),xr(n),xe(n),xc(n),optimum(n))
      simplex(:,1)=start
      do j=1,n
         simplex(:,j+1)=start
         scale=0.05_dp*max(1.0_dp,abs(start(j)))
         simplex(j,j+1)=simplex(j,j+1)+scale
      end do
      do i=1,n+1; f(i)=objective(simplex(:,i)); end do
      converged=.false.
      do iterations=1,maxit
         best=minloc(f,dim=1); worst=maxloc(f,dim=1); second=0
         do i=1,n+1
            if(i==worst)cycle
            if(second==0 .or. f(i)>f(second))second=i
         end do
         if(maxval(abs(f-f(best)))<=tol*(1.0_dp+abs(f(best))) .and. &
            maxval(abs(simplex-spread(simplex(:,best),2,n+1)))<=sqrt(tol)*(1.0_dp+maxval(abs(simplex(:,best)))))then
            converged=.true.; exit
         end if
         centroid=0.0_dp
         do i=1,n+1; if(i/=worst)centroid=centroid+simplex(:,i); end do
         centroid=centroid/real(n,dp)
         xr=2.0_dp*centroid-simplex(:,worst); fr=objective(xr)
         if(fr<f(best))then
            xe=centroid+2.0_dp*(xr-centroid); fe=objective(xe)
            if(fe<fr)then; simplex(:,worst)=xe; f(worst)=fe
            else; simplex(:,worst)=xr; f(worst)=fr; end if
         else if(fr<f(second))then
            simplex(:,worst)=xr; f(worst)=fr
         else
            if(fr<f(worst))then; xc=centroid+0.5_dp*(xr-centroid)
            else; xc=centroid+0.5_dp*(simplex(:,worst)-centroid); end if
            fc=objective(xc)
            if(fc<min(fr,f(worst)))then
               simplex(:,worst)=xc; f(worst)=fc
            else
               do i=1,n+1
                  if(i/=best)then
                     simplex(:,i)=simplex(:,best)+0.5_dp*(simplex(:,i)-simplex(:,best))
                     f(i)=objective(simplex(:,i))
                  end if
               end do
            end if
         end if
      end do
      best=minloc(f,dim=1); optimum=simplex(:,best); value=f(best)
   end subroutine nelder_mead

   subroutine numerical_hessian(objective,x,hessian)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: hessian(:,:)
      real(dp) :: xp(size(x)),xm(size(x)),xpp(size(x)),xpm(size(x)),xmp(size(x)),xmm(size(x))
      real(dp) :: hi,hj,f0
      integer :: i,j,n
      n=size(x); allocate(hessian(n,n)); f0=objective(x)
      do i=1,n
         hi=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(i)))
         xp=x; xm=x; xp(i)=xp(i)+hi; xm(i)=xm(i)-hi
         hessian(i,i)=(objective(xp)-2.0_dp*f0+objective(xm))/(hi*hi)
         do j=i+1,n
            hj=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(j)))
            xpp=x; xpm=x; xmp=x; xmm=x
            xpp(i)=xpp(i)+hi; xpp(j)=xpp(j)+hj
            xpm(i)=xpm(i)+hi; xpm(j)=xpm(j)-hj
            xmp(i)=xmp(i)-hi; xmp(j)=xmp(j)+hj
            xmm(i)=xmm(i)-hi; xmm(j)=xmm(j)-hj
            hessian(i,j)=(objective(xpp)-objective(xpm)-objective(xmp)+objective(xmm))/(4.0_dp*hi*hj)
            hessian(j,i)=hessian(i,j)
         end do
      end do
   end subroutine numerical_hessian

end module streg_optimize
