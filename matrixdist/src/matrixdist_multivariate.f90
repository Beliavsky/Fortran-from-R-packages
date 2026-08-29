! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_multivariate
   use r_compat, only: dp
   use matrixdist_linalg
   use matrixdist_ph
   use matrixdist_dph
   implicit none
   private
   public :: bivph_density, bivph_tail, bivph_laplace, bivph_moment, bivph_mean, bivph_cov
   public :: bivdph_density, bivdph_tail, bivdph_pgf, bivdph_moment, bivdph_mean, bivdph_cov
   public :: mph_density_point, mph_cdf_point, mph_laplace_point, mph_moment, mph_mean, mph_cov
   public :: mdph_density_point, mdph_pgf_point, mdph_factorial_moment, mdph_mean, mdph_cov
   public :: mphstar_mean, mphstar_cov

contains

   function factorial_real(k) result(v)
      integer,intent(in)::k
      real(dp)::v
      integer::i
      v=1.0_dp
      do i=2,k
      v=v*real(i,dp)
      end do
   end function

   function bivph_density(x1,x2,alpha,s11,s12,s22) result(f)
      real(dp),intent(in)::x1,x2,alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp)::f
      real(dp),allocatable::m1(:,:),m2(:,:),t2(:)
      if(x1<0.0_dp .or. x2<0.0_dp)then
      f=0.0_dp
      return
      end if
      m1=matrix_exponential(s11*x1)
      m2=matrix_exponential(s22*x2)
      t2=ph_exit_rates(s22)
      f=dot_product(alpha,matmul(m1,matmul(s12,matmul(m2,t2))))
   end function

   function bivph_tail(x1,x2,alpha,s11,s12,s22) result(q)
      real(dp),intent(in)::x1,x2,alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp)::q
      real(dp),allocatable::u1(:,:),m1(:,:),m2(:,:),e(:)
      allocate(e(size(s22,1)))
      e=1.0_dp
      u1=matrix_inverse(-s11)
      m1=matrix_exponential(s11*x1)
      m2=matrix_exponential(s22*x2)
      q=dot_product(alpha,matmul(u1,matmul(m1,matmul(s12,matmul(m2,e)))))
   end function

   function bivph_laplace(r1,r2,alpha,s11,s12,s22) result(v)
      real(dp),intent(in)::r1,r2,alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp)::v
      real(dp),allocatable::u1(:,:),u2(:,:),t2(:)
      u1=matrix_inverse(r1*eye_matrix(size(s11,1))-s11)
      u2=matrix_inverse(r2*eye_matrix(size(s22,1))-s22)
      t2=ph_exit_rates(s22)
      v=dot_product(alpha,matmul(u1,matmul(s12,matmul(u2,t2))))
   end function

   function bivph_moment(k1,k2,alpha,s11,s12,s22) result(v)
      integer,intent(in)::k1,k2
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp)::v
      real(dp),allocatable::u1(:,:),u2(:,:),p1(:,:),p2(:,:),e(:)
      if(k1<0 .or. k2<0 .or. k1+k2==0) error stop "bivph_moment: invalid order"
      u1=matrix_inverse(-s11)
      u2=matrix_inverse(-s22)
      p1=matrix_power(k1+1,u1)
      p2=matrix_power(k2,u2)
      allocate(e(size(s22,1)))
      e=1.0_dp
      v=factorial_real(k1)*factorial_real(k2)* &
        dot_product(alpha,matmul(p1,matmul(s12,matmul(p2,e))))
   end function

   function bivph_mean(alpha,s11,s12,s22) result(m)
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp)::m(2)
      m(1)=bivph_moment(1,0,alpha,s11,s12,s22)
      m(2)=bivph_moment(0,1,alpha,s11,s12,s22)
   end function

   function bivph_cov(alpha,s11,s12,s22) result(c)
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp)::c(2,2),m(2)
      m=bivph_mean(alpha,s11,s12,s22)
      c(1,1)=bivph_moment(2,0,alpha,s11,s12,s22)-m(1)**2
      c(2,2)=bivph_moment(0,2,alpha,s11,s12,s22)-m(2)**2
      c(1,2)=bivph_moment(1,1,alpha,s11,s12,s22)-m(1)*m(2)
      c(2,1)=c(1,2)
   end function

   function bivdph_density(k1,k2,alpha,s11,s12,s22) result(f)
      integer,intent(in)::k1,k2
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp)::f
      real(dp),allocatable::p1(:,:),p2(:,:),t2(:)
      if(k1<1 .or. k2<1)then
      f=0.0_dp
      return
      end if
      p1=matrix_power(k1-1,s11)
      p2=matrix_power(k2-1,s22)
      t2=dph_exit_probs(s22)
      f=dot_product(alpha,matmul(p1,matmul(s12,matmul(p2,t2))))
   end function

   function bivdph_tail(k1,k2,alpha,s11,s12,s22) result(q)
      integer,intent(in)::k1,k2
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp)::q
      real(dp),allocatable::p1(:,:),p2(:,:),e(:)
      p1=matrix_power(max(0,k1),s11)
      p2=matrix_power(max(0,k2),s22)
      allocate(e(size(s22,1)))
      e=1.0_dp
      q=dot_product(alpha,matmul(p1,matmul(s12,matmul(p2,e))))
   end function

   function bivdph_pgf(z1,z2,alpha,s11,s12,s22) result(v)
      real(dp),intent(in)::z1,z2,alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp)::v
      real(dp),allocatable::u1(:,:),u2(:,:),t2(:)
      u1=matrix_inverse(eye_matrix(size(s11,1))-z1*s11)
      u2=matrix_inverse(eye_matrix(size(s22,1))-z2*s22)
      t2=dph_exit_probs(s22)
      v=z1*z2*dot_product(alpha,matmul(u1,matmul(s12,matmul(u2,t2))))
   end function

   function bivdph_moment(k1,k2,alpha,s11,s12,s22) result(v)
      ! Joint factorial moment, matching matrixdist's bivdph moment method.
      integer,intent(in)::k1,k2
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp)::v
      real(dp),allocatable::a1(:,:),a2(:,:),u1(:,:),u2(:,:),b1(:,:),b2(:,:),e(:)
      if(k1<=0 .or. k2<=0) error stop "bivdph_moment: positive orders required"
      a1=matrix_power(k1-1,s11)
      u1=matrix_inverse(eye_matrix(size(s11,1))-s11)
      b1=matrix_power(k1+1,u1)
      a2=matrix_power(k2-1,s22)
      u2=matrix_inverse(eye_matrix(size(s22,1))-s22)
      b2=matrix_power(k2,u2)
      allocate(e(size(s22,1)))
      e=1.0_dp
      v=factorial_real(k1)*factorial_real(k2)* &
        dot_product(alpha,matmul(a1,matmul(b1,matmul(s12,matmul(a2,matmul(b2,e))))))
   end function

   function bivdph_mean(alpha,s11,s12,s22) result(m)
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp)::m(2)
      real(dp),allocatable::alpha2(:),u1(:,:)
      m(1)=dph_mean(alpha,s11)
      u1=matrix_inverse(eye_matrix(size(s11,1))-s11)
      alpha2=matmul(transpose(s12),matmul(transpose(u1),alpha))
      m(2)=dph_mean(alpha2,s22)
   end function

   function bivdph_cov(alpha,s11,s12,s22) result(c)
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp)::c(2,2),m(2),cross
      real(dp),allocatable::alpha2(:),u1(:,:)
      m=bivdph_mean(alpha,s11,s12,s22)
      c(1,1)=dph_variance(alpha,s11)
      u1=matrix_inverse(eye_matrix(size(s11,1))-s11)
      alpha2=matmul(transpose(s12),matmul(transpose(u1),alpha))
      c(2,2)=dph_variance(alpha2,s22)
      cross=bivdph_moment(1,1,alpha,s11,s12,s22)
      c(1,2)=cross-m(1)*m(2)
      c(2,1)=c(1,2)
   end function

   function mph_density_point(y,alpha,s,delta) result(v)
      real(dp),intent(in)::y(:),alpha(:),s(:,:,:)
      logical,intent(in),optional::delta(:)
      real(dp)::v,prodv,f
      real(dp),allocatable::unit(:)
      integer::p,d,j,i
      logical::unc
      p=size(alpha)
      d=size(s,3)
      if(size(y)/=d) error stop "mph_density_point: dimension mismatch"
      v=0.0_dp
      allocate(unit(p))
      do j=1,p
         unit=0.0_dp
         unit(j)=1.0_dp
         prodv=1.0_dp
         do i=1,d
            unc=.true.
            if(present(delta))unc=delta(i)
            if(unc)then
            f=ph_density(y(i),unit,s(:,:,i))
            else
            f=ph_survival(y(i),unit,s(:,:,i))
            end if
            prodv=prodv*f
         end do
         v=v+alpha(j)*prodv
      end do
   end function

   function mph_cdf_point(y,alpha,s,lower_tail) result(v)
      real(dp),intent(in)::y(:),alpha(:),s(:,:,:)
      logical,intent(in),optional::lower_tail
      real(dp)::v,prodv
      real(dp),allocatable::unit(:)
      integer::p,d,j,i
      logical::lower
      lower=.true.
      if(present(lower_tail))lower=lower_tail
      p=size(alpha)
      d=size(s,3)
      v=0.0_dp
      allocate(unit(p))
      do j=1,p
         unit=0.0_dp
         unit(j)=1.0_dp
         prodv=1.0_dp
         do i=1,d
         prodv=prodv*ph_cdf(y(i),unit,s(:,:,i),lower)
         end do
         v=v+alpha(j)*prodv
      end do
   end function

   function mph_laplace_point(r,alpha,s) result(v)
      real(dp),intent(in)::r(:),alpha(:),s(:,:,:)
      real(dp)::v,prodv
      real(dp),allocatable::unit(:)
      integer::p,d,j,i
      p=size(alpha)
      d=size(s,3)
      v=0.0_dp
      allocate(unit(p))
      do j=1,p
         unit=0.0_dp
         unit(j)=1.0_dp
         prodv=1.0_dp
         do i=1,d
         prodv=prodv*ph_laplace(r(i),unit,s(:,:,i))
         end do
         v=v+alpha(j)*prodv
      end do
   end function

   function mph_moment(k,alpha,s) result(v)
      integer,intent(in)::k(:)
      real(dp),intent(in)::alpha(:),s(:,:,:)
      real(dp)::v,prodv
      real(dp),allocatable::unit(:),u(:,:),pk(:,:),e(:)
      integer::p,d,j,i
      p=size(alpha)
      d=size(s,3)
      if(size(k)/=d)error stop "mph_moment: dimension mismatch"
      allocate(unit(p),e(p))
      e=1.0_dp
      v=0.0_dp
      do j=1,p
         unit=0.0_dp
         unit(j)=1.0_dp
         prodv=1.0_dp
         do i=1,d
            if(k(i)==0) cycle
            u=matrix_inverse(-s(:,:,i))
            pk=matrix_power(k(i),u)
            prodv=prodv*factorial_real(k(i))*dot_product(unit,matmul(pk,e))
         end do
         v=v+alpha(j)*prodv
      end do
   end function

   function mph_mean(alpha,s) result(m)
      real(dp),intent(in)::alpha(:),s(:,:,:)
      real(dp)::m(size(s,3))
      integer::i,d
      integer,allocatable::k(:)
      d=size(s,3)
      allocate(k(d))
      k=0
      do i=1,d
      k=0
      k(i)=1
      m(i)=mph_moment(k,alpha,s)
      end do
   end function

   function mph_cov(alpha,s) result(c)
      real(dp),intent(in)::alpha(:),s(:,:,:)
      real(dp)::c(size(s,3),size(s,3)),m(size(s,3))
      integer,allocatable::k(:)
      integer::i,j,d
      d=size(s,3)
      allocate(k(d))
      m=mph_mean(alpha,s)
      do i=1,d
         do j=i,d
            k=0
            k(i)=k(i)+1
            k(j)=k(j)+1
            c(i,j)=mph_moment(k,alpha,s)-m(i)*m(j)
            c(j,i)=c(i,j)
         end do
      end do
   end function

   function mdph_density_point(k,alpha,s) result(v)
      integer,intent(in)::k(:)
      real(dp),intent(in)::alpha(:),s(:,:,:)
      real(dp)::v,prodv
      real(dp),allocatable::unit(:)
      integer::p,d,j,i
      p=size(alpha)
      d=size(s,3)
      v=0.0_dp
      allocate(unit(p))
      do j=1,p
         unit=0.0_dp
         unit(j)=1.0_dp
         prodv=1.0_dp
         do i=1,d
         prodv=prodv*dph_density(k(i),unit,s(:,:,i))
         end do
         v=v+alpha(j)*prodv
      end do
   end function

   function mdph_pgf_point(z,alpha,s) result(v)
      real(dp),intent(in)::z(:),alpha(:),s(:,:,:)
      real(dp)::v,prodv
      real(dp),allocatable::unit(:)
      integer::p,d,j,i
      p=size(alpha)
      d=size(s,3)
      v=0.0_dp
      allocate(unit(p))
      do j=1,p
         unit=0.0_dp
         unit(j)=1.0_dp
         prodv=1.0_dp
         do i=1,d
         prodv=prodv*dph_pgf(z(i),unit,s(:,:,i))
         end do
         v=v+alpha(j)*prodv
      end do
   end function

   function mdph_factorial_moment(k,alpha,s) result(v)
      integer,intent(in)::k(:)
      real(dp),intent(in)::alpha(:),s(:,:,:)
      real(dp)::v,prodv
      real(dp),allocatable::unit(:)
      integer::p,d,j,i
      p=size(alpha)
      d=size(s,3)
      v=0.0_dp
      allocate(unit(p))
      do j=1,p
         unit=0.0_dp
         unit(j)=1.0_dp
         prodv=1.0_dp
         do i=1,d
         prodv=prodv*dph_factorial_moment(k(i),unit,s(:,:,i))
         end do
         v=v+alpha(j)*prodv
      end do
   end function

   function mdph_mean(alpha,s) result(m)
      real(dp),intent(in)::alpha(:),s(:,:,:)
      real(dp)::m(size(s,3))
      integer::i
      do i=1,size(s,3)
      m(i)=dph_mean(alpha,s(:,:,i))
      end do
   end function

   function mdph_cov(alpha,s) result(c)
      real(dp),intent(in)::alpha(:),s(:,:,:)
      real(dp)::c(size(s,3),size(s,3)),m(size(s,3)),cross,prodv
      real(dp),allocatable::unit(:)
      integer::i,j,h,p,d
      p=size(alpha)
      d=size(s,3)
      m=mdph_mean(alpha,s)
      allocate(unit(p))
      do i=1,d
         c(i,i)=dph_variance(alpha,s(:,:,i))
         do j=i+1,d
            cross=0.0_dp
            do h=1,p
               unit=0.0_dp
               unit(h)=1.0_dp
               prodv=dph_mean(unit,s(:,:,i))*dph_mean(unit,s(:,:,j))
               cross=cross+alpha(h)*prodv
            end do
            c(i,j)=cross-m(i)*m(j)
            c(j,i)=c(i,j)
         end do
      end do
   end function

   function mphstar_mean(alpha,s,reward) result(m)
      real(dp),intent(in)::alpha(:),s(:,:),reward(:,:)
      real(dp)::m(size(reward,2))
      real(dp),allocatable::u(:,:)
      integer::j
      u=matrix_inverse(-s)
      do j=1,size(reward,2)
      m(j)=dot_product(alpha,matmul(u,reward(:,j)))
      end do
   end function

   function mphstar_cov(alpha,s,reward) result(c)
      real(dp),intent(in)::alpha(:),s(:,:),reward(:,:)
      real(dp)::c(size(reward,2),size(reward,2)),m(size(reward,2)),cross
      real(dp),allocatable::u(:,:),ri(:,:),rj(:,:)
      integer::i,j,p,d,k
      p=size(alpha)
      d=size(reward,2)
      u=matrix_inverse(-s)
      m=mphstar_mean(alpha,s,reward)
      allocate(ri(p,p),rj(p,p))
      ri=0.0_dp
      rj=0.0_dp
      do i=1,d
         do j=i,d
            ri=0.0_dp
            rj=0.0_dp
            do k=1,p
            ri(k,k)=reward(k,i)
            rj(k,k)=reward(k,j)
            end do
            cross=dot_product(alpha,matmul(u,matmul(ri,matmul(u,reward(:,j))))) + &
                  dot_product(alpha,matmul(u,matmul(rj,matmul(u,reward(:,i)))))
            c(i,j)=cross-m(i)*m(j)
            c(j,i)=c(i,j)
         end do
      end do
   end function

end module matrixdist_multivariate
