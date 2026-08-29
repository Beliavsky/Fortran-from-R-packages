module spatialextremes_utils
   use spatialextremes_base, only: dp
   implicit none
   private
   public :: stirling_second_kind,bell_number,vandercorput,penalization,finite_difference_hessian
contains
   pure real(dp) function stirling_second_kind(n,k) result(s)
      integer,intent(in)::n,k
      real(dp)::a(0:n,0:k)
      integer::i,j
      a=0.0_dp
      a(0,0)=1.0_dp
      do i=1,n
         do j=1,min(i,k)
         a(i,j)=a(i-1,j-1)+real(j,dp)*a(i-1,j)
         end do
      end do
      s=a(n,k)
   end function

   pure integer function bell_number(n) result(b)
      integer,intent(in)::n
      integer::k
      real(dp)::s
      s=0
      do k=0,n
      s=s+stirling_second_kind(n,k)
      end do
      if(s>real(huge(b),dp))then
      b=huge(b)
      else
      b=nint(s)
      end if
   end function

   pure real(dp) function vandercorput(index,base) result(x)
      integer,intent(in)::index,base
      integer::n
      real(dp)::den
      n=index
      x=0
      den=1
      do while(n>0)
      den=den*base
      x=x+mod(n,base)/den
      n=n/base
      end do
   end function

   pure real(dp) function penalization(pmat,coef,lambda,nparametric) result(pen)
      real(dp),intent(in)::pmat(:,:),coef(:),lambda
      integer,intent(in),optional::nparametric
      integer::np
      real(dp),allocatable::c(:)
      np=0
      if(present(nparametric))np=nparametric
      if(np>=size(coef))then
      pen=0
      return
      end if
      c=coef(np+1:)
      if(size(pmat,1)==size(c))then
      pen=lambda*dot_product(c,matmul(pmat,c))
      else
      pen=0
      end if
   end function

   subroutine finite_difference_hessian(fn,x,hess,step)
      abstract interface
         function f_iface(z) result(v)
            import dp
            real(dp),intent(in)::z(:)
            real(dp)::v
         end function
      end interface
      procedure(f_iface)::fn
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::hess(size(x),size(x))
      real(dp),intent(in),optional::step
      real(dp)::h,fi,fpp,fpm,fmp,fmm
      real(dp)::xp(size(x)),xm(size(x)),xpp(size(x)),xpm(size(x)),xmp(size(x)),xmm(size(x))
      integer::i,j
      h=1.0e-5_dp
      if(present(step))h=step
      fi=fn(x)
      hess=0
      do i=1,size(x)
         xp=x
         xm=x
         xp(i)=xp(i)+h
         xm(i)=xm(i)-h
         hess(i,i)=(fn(xp)-2*fi+fn(xm))/(h*h)
         do j=i+1,size(x)
            xpp=x
            xpm=x
            xmp=x
            xmm=x
            xpp(i)=xpp(i)+h
            xpp(j)=xpp(j)+h
            xpm(i)=xpm(i)+h
            xpm(j)=xpm(j)-h
            xmp(i)=xmp(i)-h
            xmp(j)=xmp(j)+h
            xmm(i)=xmm(i)-h
            xmm(j)=xmm(j)-h
            fpp=fn(xpp)
            fpm=fn(xpm)
            fmp=fn(xmp)
            fmm=fn(xmm)
            hess(i,j)=(fpp-fpm-fmp+fmm)/(4*h*h)
            hess(j,i)=hess(i,j)
         end do
      end do
   end subroutine
end module spatialextremes_utils
