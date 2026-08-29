module spatialextremes_copula
   use spatialextremes_base, only: dp,pi,chol_upper,logdet_spd,nan_dp,chisq_rand
   use r_compat, only: dnorm,dt,qnorm,qt,rnorm1
   implicit none
   private
   public :: gaussian_copula_loglik,student_copula_loglik,simulate_gaussian_copula,simulate_student_copula
contains
   real(dp) function gaussian_copula_loglik(u,cor) result(ll)
      real(dp),intent(in)::u(:,:),cor(:,:)
      ll=gaussian_copula_loglik_core(u,cor)
   end function gaussian_copula_loglik

   real(dp) function gaussian_copula_loglik_core(u,cor) result(ll)
      real(dp),intent(in)::u(:,:),cor(:,:)
      real(dp)::r(size(cor,1),size(cor,2)),z(size(cor,1)),z0(size(cor,1)),ldet
      integer::info,i,j,n,p
      n=size(u,1)
      p=size(u,2)
      call chol_upper(cor,r,info)
      if(info/=0)then
      ll=-huge(1.0_dp)
      return
      end if
      ldet=0
      do j=1,p
      ldet=ldet+2*log(r(j,j))
      end do
      ll=0
      do i=1,n
         do j=1,p
         z0(j)=qnorm(u(i,j),0.0_dp,1.0_dp,.true.)
         end do
         z=z0
         call triangular_whiten(r,z)
         ll=ll-0.5_dp*ldet-0.5_dp*(sum(z*z)-sum(z0*z0))
      end do
   end function gaussian_copula_loglik_core

   real(dp) function student_copula_loglik(u,nu,cor) result(ll)
      real(dp),intent(in)::u(:,:),nu,cor(:,:)
      real(dp)::r(size(cor,1),size(cor,2)),z(size(cor,1)),z0(size(cor,1)),ldet,q,c
      integer::info,i,j,n,p
      n=size(u,1)
      p=size(u,2)
      call chol_upper(cor,r,info)
      if(info/=0)then
      ll=-huge(1.0_dp)
      return
      end if
      ldet=0
      do j=1,p
      ldet=ldet+2*log(r(j,j))
      end do
      c=log_gamma(0.5_dp*(nu+p))-log_gamma(0.5_dp*nu)-0.5_dp*p*log(nu*pi)
      ll=0
      do i=1,n
         do j=1,p
         z0(j)=qt(u(i,j),nu)
         end do
         z=z0
         call triangular_whiten(r,z)
         q=sum(z*z)
         ll=ll+c-0.5_dp*ldet-0.5_dp*(nu+p)*log(1.0_dp+q/nu)
         do j=1,p
         ll=ll-log(dt(z0(j),nu,.false.))
         end do
      end do
   end function student_copula_loglik

   function simulate_gaussian_copula(n,cor) result(u)
      real(dp),intent(in)::cor(:,:)
      integer,intent(in)::n
      real(dp)::u(n,size(cor,1)),r(size(cor,1),size(cor,2)),z(size(cor,1))
      integer::info,i,j
      call chol_upper(cor,r,info)
      if(info/=0)then
      u=nan_dp()
      return
      end if
      do i=1,n
         do j=1,size(cor,1)
         z(j)=rnorm1()
         end do
         z=matmul(transpose(r),z)
         do j=1,size(cor,1)
         u(i,j)=0.5_dp*erfc(-z(j)/sqrt(2.0_dp))
         end do
      end do
   end function simulate_gaussian_copula

   function simulate_student_copula(n,nu,cor) result(u)
      real(dp),intent(in)::cor(:,:),nu
      integer,intent(in)::n
      real(dp)::u(n,size(cor,1)),r(size(cor,1),size(cor,2)),z(size(cor,1)),s
      integer::info,i,j
      call chol_upper(cor,r,info)
      if(info/=0)then
      u=nan_dp()
      return
      end if
      do i=1,n
         do j=1,size(cor,1)
         z(j)=rnorm1()
         end do
         z=matmul(transpose(r),z)
         s=sqrt(chisq_rand(nu)/nu)
         z=z/s
         do j=1,size(cor,1)
         u(i,j)=student_cdf(z(j),nu)
         end do
      end do
   end function simulate_student_copula

   subroutine triangular_whiten(r,z)
      real(dp),intent(in)::r(:,:)
      real(dp),intent(inout)::z(:)
      integer::i,j,n
      n=size(z)
      ! Solve R^T x = z (R upper)
      do i=1,n
         if(i>1)z(i)=z(i)-dot_product(r(1:i-1,i),z(1:i-1))
         z(i)=z(i)/r(i,i)
      end do
   end subroutine triangular_whiten

   pure real(dp) function student_cdf(x,nu) result(p)
      use r_compat, only: pt
      real(dp),intent(in)::x,nu
      p=pt(x,nu)
   end function
end module spatialextremes_copula
