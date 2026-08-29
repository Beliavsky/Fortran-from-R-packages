module lavaan_linalg
   use lavaan_kinds, only : dp
   implicit none
   private
   public :: chol_lower, inverse_spd, inverse_general, logdet_spd, solve_linear
   public :: sym_eigen_jacobi, symmetric_sqrt, trace_matrix, sample_mean_cov
   public :: vec, vech, vech_reverse, commutation_matrix, duplication_matrix
   public :: orthogonal_complement
contains
   subroutine chol_lower(a, l, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: l(:, :)
      integer, intent(out) :: info
      integer :: n, i, j, k
      real(dp) :: s
      n = size(a,1)
      allocate(l(n,n))
      l = 0.0_dp
      info = 0
      if (size(a,2) /= n) then
      info=-1
      return
      end if
      do i=1,n
         do j=1,i
            s=a(i,j)
            do k=1,j-1
               s=s-l(i,k)*l(j,k)
            end do
            if (i==j) then
               if (s <= 0.0_dp) then
               info=i
               return
               end if
               l(i,j)=sqrt(s)
            else
               l(i,j)=s/l(j,j)
            end if
         end do
      end do
   end subroutine chol_lower

   subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: aa(:, :), bb(:), row(:)
      real(dp) :: fac, piv
      integer :: n, i, j, k, p
      n=size(b)
      allocate(aa(n,n),bb(n),x(n),row(n))
      aa=a
      bb=b
      info=0
      if(size(a,1)/=n .or. size(a,2)/=n) then
      info=-1
      x=0
      return
      end if
      do k=1,n-1
         p=k
         do i=k+1,n
            if(abs(aa(i,k))>abs(aa(p,k))) p=i
         end do
         if(abs(aa(p,k))<tiny(1.0_dp)) then
         info=k
         x=0
         return
         end if
         if(p/=k) then
            row=aa(k,:)
            aa(k,:)=aa(p,:)
            aa(p,:)=row
            piv=bb(k)
            bb(k)=bb(p)
            bb(p)=piv
         end if
         do i=k+1,n
            fac=aa(i,k)/aa(k,k)
            aa(i,k)=0.0_dp
            do j=k+1,n
            aa(i,j)=aa(i,j)-fac*aa(k,j)
            end do
            bb(i)=bb(i)-fac*bb(k)
         end do
      end do
      if(abs(aa(n,n))<tiny(1.0_dp)) then
      info=n
      x=0
      return
      end if
      do i=n,1,-1
         x(i)=bb(i)
         do j=i+1,n
         x(i)=x(i)-aa(i,j)*x(j)
         end do
         x(i)=x(i)/aa(i,i)
      end do
   end subroutine solve_linear

   subroutine inverse_general(a, ainv, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: ainv(:, :)
      integer, intent(out) :: info
      integer :: n,j,istat
      real(dp), allocatable :: e(:),x(:)
      n=size(a,1)
      allocate(ainv(n,n),e(n))
      ainv=0
      info=0
      do j=1,n
         e=0
         e(j)=1
         call solve_linear(a,e,x,istat)
         if(istat/=0) then
         info=istat
         return
         end if
         ainv(:,j)=x
      end do
   end subroutine inverse_general

   subroutine inverse_spd(a, ainv, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: ainv(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: l(:, :), y(:), x(:)
      integer :: n,i,j,k
      call chol_lower(a,l,info)
      n=size(a,1)
      allocate(ainv(n,n),y(n),x(n))
      ainv=0
      if(info/=0) return
      do j=1,n
         y=0
         x=0
         do i=1,n
            y(i)=merge(1.0_dp,0.0_dp,i==j)
            do k=1,i-1
            y(i)=y(i)-l(i,k)*y(k)
            end do
            y(i)=y(i)/l(i,i)
         end do
         do i=n,1,-1
            x(i)=y(i)
            do k=i+1,n
            x(i)=x(i)-l(k,i)*x(k)
            end do
            x(i)=x(i)/l(i,i)
         end do
         ainv(:,j)=x
      end do
   end subroutine inverse_spd

   function logdet_spd(a, info) result(v)
      real(dp), intent(in) :: a(:, :)
      integer, intent(out) :: info
      real(dp) :: v
      real(dp), allocatable :: l(:, :)
      integer :: i
      call chol_lower(a,l,info)
      if(info/=0) then
      v=huge(1.0_dp)
      return
      end if
      v=0
      do i=1,size(a,1)
      v=v+2.0_dp*log(l(i,i))
      end do
   end function logdet_spd

   subroutine sym_eigen_jacobi(a, values, vectors, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: b(:, :)
      real(dp) :: app,aqq,apq,tau,t,c,s,tmp,maxoff
      integer :: n,i,j,p,q,iter
      n=size(a,1)
      allocate(b(n,n),values(n),vectors(n,n))
      b=a
      vectors=0
      maxoff=0.0_dp
      do i=1,n
      vectors(i,i)=1
      end do
      info=0
      do iter=1,100*n*n
         maxoff=0
         p=1
         q=min(2,n)
         do i=1,n-1
         do j=i+1,n
            if(abs(b(i,j))>maxoff) then
            maxoff=abs(b(i,j))
            p=i
            q=j
            end if
         end do
         end do
         if(maxoff<1.0e-13_dp) exit
         app=b(p,p)
         aqq=b(q,q)
         apq=b(p,q)
         tau=(aqq-app)/(2*apq)
         if(tau>=0) then
         t=1/(tau+sqrt(1+tau*tau))
         else
         t=-1/(-tau+sqrt(1+tau*tau))
         end if
         c=1/sqrt(1+t*t)
         s=t*c
         do j=1,n
            if(j/=p .and. j/=q) then
               tmp=b(j,p)
               b(j,p)=c*tmp-s*b(j,q)
               b(p,j)=b(j,p)
               b(j,q)=s*tmp+c*b(j,q)
               b(q,j)=b(j,q)
            end if
         end do
         b(p,p)=c*c*app-2*s*c*apq+s*s*aqq
         b(q,q)=s*s*app+2*s*c*apq+c*c*aqq
         b(p,q)=0
         b(q,p)=0
         do j=1,n
            tmp=vectors(j,p)
            vectors(j,p)=c*tmp-s*vectors(j,q)
            vectors(j,q)=s*tmp+c*vectors(j,q)
         end do
      end do
      if(maxoff>=1.0e-10_dp) info=1
      do i=1,n
      values(i)=b(i,i)
      end do
      do i=1,n-1
      do j=i+1,n
         if(values(j)<values(i)) then
            tmp=values(i)
            values(i)=values(j)
            values(j)=tmp
            b(:,1)=vectors(:,i)
            vectors(:,i)=vectors(:,j)
            vectors(:,j)=b(:,1)
         end if
      end do
      end do
   end subroutine sym_eigen_jacobi

   subroutine symmetric_sqrt(a, root, inverse, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: root(:, :)
      logical, intent(in), optional :: inverse
      integer, intent(out) :: info
      real(dp), allocatable :: val(:),vec(:, :),d(:, :)
      logical :: inv
      integer :: i,n
      inv=.false.
      if(present(inverse)) inv=inverse
      call sym_eigen_jacobi(a,val,vec,info)
      n=size(a,1)
      allocate(d(n,n))
      d=0
      if(info/=0) then
      allocate(root(n,n))
      root=0
      return
      end if
      do i=1,n
         if(val(i)<=0) then
         info=i
         allocate(root(n,n))
         root=0
         return
         end if
         if(inv) then
         d(i,i)=1/sqrt(val(i))
         else
         d(i,i)=sqrt(val(i))
         end if
      end do
      root=matmul(vec,matmul(d,transpose(vec)))
   end subroutine symmetric_sqrt

   pure function trace_matrix(a) result(v)
      real(dp), intent(in) :: a(:, :)
      real(dp) :: v
      integer :: i
      v=0
      do i=1,min(size(a,1),size(a,2))
      v=v+a(i,i)
      end do
   end function trace_matrix

   subroutine sample_mean_cov(x, mean, cov, unbiased)
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: mean(:), cov(:, :)
      logical, intent(in), optional :: unbiased
      logical :: ub
      integer :: n,p,i
      real(dp) :: den
      real(dp), allocatable :: z(:)
      n=size(x,1)
      p=size(x,2)
      ub=.false.
      if(present(unbiased)) ub=unbiased
      allocate(mean(p),cov(p,p),z(p))
      mean=sum(x,dim=1)/real(n,dp)
      cov=0
      do i=1,n
      z=x(i,:)-mean
      cov=cov+spread(z,2,p)*spread(z,1,p)
      end do
      den=real(n,dp)
      if(ub .and. n>1) den=real(n-1,dp)
      cov=cov/den
   end subroutine sample_mean_cov

   function vec(a) result(v)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable :: v(:)
      integer :: i,j,k
      allocate(v(size(a)))
      k=0
      do j=1,size(a,2)
      do i=1,size(a,1)
      k=k+1
      v(k)=a(i,j)
      end do
      end do
   end function vec

   function vech(a) result(v)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable :: v(:)
      integer :: n,i,j,k
      n=size(a,1)
      allocate(v(n*(n+1)/2))
      k=0
      do j=1,n
      do i=j,n
      k=k+1
      v(k)=a(i,j)
      end do
      end do
   end function vech

   function vech_reverse(v,n) result(a)
      real(dp), intent(in) :: v(:)
      integer,intent(in)::n
      real(dp),allocatable::a(:,:)
      integer::i,j,k
      allocate(a(n,n))
      a=0
      k=0
      do j=1,n
      do i=j,n
      k=k+1
      a(i,j)=v(k)
      a(j,i)=v(k)
      end do
      end do
   end function vech_reverse

   function commutation_matrix(m,n) result(kmat)
      integer,intent(in)::m,n
      real(dp),allocatable::kmat(:,:)
      integer::i,j,r,c
      allocate(kmat(m*n,m*n))
      kmat=0
      do j=1,n
      do i=1,m
      c=i+(j-1)*m
      r=j+(i-1)*n
      kmat(r,c)=1
      end do
      end do
   end function commutation_matrix

   function duplication_matrix(n) result(dmat)
      integer,intent(in)::n
      real(dp),allocatable::dmat(:,:)
      integer::i,j,k,r1,r2
      allocate(dmat(n*n,n*(n+1)/2))
      dmat=0
      k=0
      do j=1,n
      do i=j,n
         k=k+1
         r1=i+(j-1)*n
         dmat(r1,k)=1
         if(i/=j) then
         r2=j+(i-1)*n
         dmat(r2,k)=1
         end if
      end do
      end do
   end function duplication_matrix

   subroutine orthogonal_complement(a, q2, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: q2(:, :)
      integer,intent(out)::info
      integer :: m,n,i,j,k,r
      real(dp),allocatable::q(:,:),v(:)
      real(dp)::nr
      m=size(a,1)
      n=size(a,2)
      allocate(q(m,m),v(m))
      q=0
      k=0
      info=0
      do j=1,n
         v=a(:,j)
         do i=1,k
         v=v-dot_product(q(:,i),v)*q(:,i)
         end do
         nr=sqrt(dot_product(v,v))
         if(nr>1e-12_dp) then
         k=k+1
         q(:,k)=v/nr
         end if
      end do
      r=k
      do j=1,m
         v=0
         v(j)=1
         do i=1,k
         v=v-dot_product(q(:,i),v)*q(:,i)
         end do
         nr=sqrt(dot_product(v,v))
         if(nr>1e-10_dp) then
         k=k+1
         q(:,k)=v/nr
         end if
         if(k==m) exit
      end do
      if(k<m) then
      info=1
      allocate(q2(m,0))
      return
      end if
      allocate(q2(m,m-r))
      if(m-r>0) q2=q(:,r+1:m)
   end subroutine orthogonal_complement
end module lavaan_linalg
