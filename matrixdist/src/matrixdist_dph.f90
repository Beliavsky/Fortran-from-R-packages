! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_dph
   use r_compat, only: dp
   use matrixdist_linalg
   use matrixdist_types, only: dph_type
   implicit none
   private
   public :: make_dph, dph_density, dph_cdf, dph_survival, dph_pgf
   public :: dph_factorial_moment, dph_mean, dph_variance, dph_loglik
   public :: dph_sum, dph_mixture, dph_minimum, dph_maximum, dph_exit_probs

contains

   function make_dph(alpha,s) result(x)
      real(dp), intent(in) :: alpha(:),s(:,:)
      type(dph_type) :: x
      if (size(s,1)/=size(s,2) .or. size(alpha)/=size(s,1)) error stop "make_dph: dimensions"
      x%alpha=alpha
      x%s=s
   end function make_dph

   function dph_exit_probs(s) result(t)
      real(dp), intent(in) :: s(:,:)
      real(dp) :: t(size(s,1))
      t=1.0_dp-sum(s,dim=2)
   end function dph_exit_probs

   function dph_density(k,alpha,s) result(f)
      integer, intent(in) :: k
      real(dp), intent(in) :: alpha(:),s(:,:)
      real(dp) :: f
      real(dp), allocatable :: t(:),p(:,:)
      if (k < 1) then
         f=0.0_dp
         return
      end if
      t=dph_exit_probs(s)
      p=matrix_power(k-1,s)
      f=dot_product(alpha,matmul(p,t))
      if (f < 0.0_dp .and. abs(f)<100.0_dp*epsilon(f)) f=0.0_dp
   end function dph_density

   function dph_survival(k,alpha,s) result(q)
      integer, intent(in) :: k
      real(dp), intent(in) :: alpha(:),s(:,:)
      real(dp) :: q
      real(dp), allocatable :: p(:,:),e(:)
      if (k < 0) then
         q=1.0_dp
         return
      end if
      allocate(e(size(alpha)))
      e=1.0_dp
      p=matrix_power(k,s)
      q=dot_product(alpha,matmul(p,e))
      q=min(1.0_dp,max(0.0_dp,q))
   end function dph_survival

   function dph_cdf(k,alpha,s,lower_tail) result(pv)
      integer, intent(in) :: k
      real(dp), intent(in) :: alpha(:),s(:,:)
      logical, intent(in), optional :: lower_tail
      real(dp) :: pv,q
      logical :: lower
      lower=.true.
      if(present(lower_tail)) lower=lower_tail
      q=dph_survival(k,alpha,s)
      if(lower) then
      pv=1.0_dp-q
      else
      pv=q
      end if
   end function dph_cdf

   function dph_pgf(z,alpha,s) result(g)
      real(dp), intent(in) :: z,alpha(:),s(:,:)
      real(dp) :: g
      real(dp), allocatable :: a(:,:),t(:),v(:)
      a=eye_matrix(size(alpha))-z*s
      t=dph_exit_probs(s)
      v=solve_vector(a,t)
      g=z*dot_product(alpha,v)
   end function dph_pgf

   function dph_factorial_moment(k,alpha,s) result(mom)
      integer,intent(in)::k
      real(dp),intent(in)::alpha(:),s(:,:)
      real(dp)::mom,fact
      real(dp),allocatable::m1(:,:),m2(:,:),u(:,:),e(:)
      integer::j,p
      if(k<=0) error stop "dph_factorial_moment: k must be positive"
      p=size(alpha)
      allocate(e(p))
      e=1.0_dp
      m1=matrix_power(k-1,s)
      u=matrix_inverse(eye_matrix(p)-s)
      m2=matrix_power(k,u)
      fact=1.0_dp
      do j=2,k
      fact=fact*real(j,dp)
      end do
      mom=fact*dot_product(alpha,matmul(m1,matmul(m2,e)))
   end function dph_factorial_moment

   function dph_mean(alpha,s) result(m)
      real(dp),intent(in)::alpha(:),s(:,:)
      real(dp)::m
      real(dp),allocatable::u(:,:),e(:)
      allocate(e(size(alpha)))
      e=1.0_dp
      u=matrix_inverse(eye_matrix(size(alpha))-s)
      m=dot_product(alpha,matmul(u,e))
   end function dph_mean

   function dph_variance(alpha,s) result(v)
      real(dp),intent(in)::alpha(:),s(:,:)
      real(dp)::v,m,sm
      real(dp),allocatable::u(:,:),u2(:,:),e(:)
      allocate(e(size(alpha)))
      e=1.0_dp
      u=matrix_inverse(eye_matrix(size(alpha))-s)
      u2=matmul(u,u)
      m=dot_product(alpha,matmul(u,e))
      sm=2.0_dp*dot_product(alpha,matmul(s,matmul(u2,e)))
      v=sm+m-m*m
   end function dph_variance

   function dph_loglik(alpha,s,obs,weight) result(ll)
      real(dp),intent(in)::alpha(:),s(:,:)
      integer,intent(in)::obs(:)
      real(dp),intent(in),optional::weight(:)
      real(dp)::ll,f,w
      integer::i
      ll=0.0_dp
      do i=1,size(obs)
         w=1.0_dp
         if(present(weight))w=weight(i)
         f=dph_density(obs(i),alpha,s)
         if(f<=0.0_dp) then
         ll=-huge(1.0_dp)
         return
         end if
         ll=ll+w*log(f)
      end do
   end function dph_loglik

   function dph_sum(x1,x2) result(z)
      type(dph_type),intent(in)::x1,x2
      type(dph_type)::z
      integer::p1,p2
      real(dp),allocatable::t1(:)
      p1=size(x1%alpha)
      p2=size(x2%alpha)
      allocate(z%alpha(p1+p2),z%s(p1+p2,p1+p2))
      z%s=0.0_dp
      z%alpha(1:p1)=x1%alpha
      z%alpha(p1+1:)=0.0_dp
      z%s(1:p1,1:p1)=x1%s
      z%s(p1+1:,p1+1:)=x2%s
      t1=dph_exit_probs(x1%s)
      z%s(1:p1,p1+1:)=spread(t1,2,p2)*spread(x2%alpha,1,p1)
   end function dph_sum

   function dph_mixture(x1,x2,prob) result(z)
      type(dph_type),intent(in)::x1,x2
      real(dp),intent(in)::prob
      type(dph_type)::z
      integer::p1,p2
      p1=size(x1%alpha)
      p2=size(x2%alpha)
      allocate(z%alpha(p1+p2),z%s(p1+p2,p1+p2))
      z%s=0.0_dp
      z%alpha=[prob*x1%alpha,(1.0_dp-prob)*x2%alpha]
      z%s(1:p1,1:p1)=x1%s
      z%s(p1+1:,p1+1:)=x2%s
   end function dph_mixture

   function dph_minimum(x1,x2) result(z)
      type(dph_type),intent(in)::x1,x2
      type(dph_type)::z
      real(dp),allocatable::a1(:,:),a2(:,:),aa(:,:)
      allocate(a1(1,size(x1%alpha)),a2(1,size(x2%alpha)))
      a1(1,:)=x1%alpha
      a2(1,:)=x2%alpha
      aa=kronecker(a1,a2)
      z%alpha=aa(1,:)
      z%s=kronecker(x1%s,x2%s)
   end function dph_minimum

   function dph_maximum(x1,x2) result(z)
      type(dph_type),intent(in)::x1,x2
      type(dph_type)::z
      integer::p1,p2,n12,n
      real(dp),allocatable::a1(:,:),a2(:,:),aa(:,:),t1(:),t2(:)
      p1=size(x1%alpha)
      p2=size(x2%alpha)
      n12=p1*p2
      n=n12+p1+p2
      allocate(a1(1,p1),a2(1,p2))
      a1(1,:)=x1%alpha
      a2(1,:)=x2%alpha
      aa=kronecker(a1,a2)
      allocate(z%alpha(n),z%s(n,n))
      z%alpha=0.0_dp
      z%s=0.0_dp
      z%alpha(1:n12)=aa(1,:)
      z%s(1:n12,1:n12)=kronecker(x1%s,x2%s)
      t1=dph_exit_probs(x1%s)
      t2=dph_exit_probs(x2%s)
      z%s(1:n12,n12+1:n12+p1)=kronecker(x1%s,reshape(t2,[p2,1]))
      z%s(1:n12,n12+p1+1:n)=kronecker(reshape(t1,[p1,1]),x2%s)
      z%s(n12+1:n12+p1,n12+1:n12+p1)=x1%s
      z%s(n12+p1+1:n,n12+p1+1:n)=x2%s
   end function dph_maximum

end module matrixdist_dph
