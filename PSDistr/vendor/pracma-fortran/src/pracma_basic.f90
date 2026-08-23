! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_basic
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use pracma_kinds, only : dp, pi_dp, eps_dp, huge_dp
   implicit none
   private

   public :: linspace, logspace, logseq, eye, ones, zeros, Diag, blkdiag
   public :: hilb, pascal, vander, hankel, Toeplitz, flipud, fliplr, rot90
   public :: circshift, repmat, meshgrid, size2, numel, nnz, isempty
   public :: dot, cross, crossn, Norm, fnorm, Trace, harmmean, geomean
   public :: trimmean, Mode, std, std_err, ceil, Fix, mod, rem, idivide
   public :: gcd, Lcm, nextpow2, pow2, distmat, pdist, pdist2, hausdorff_dist
   public :: accumarray, bsxfun, uniq, find, finds, findintervals
   public :: deg2rad, rad2deg, sind, cosd, tand, cotd, asind, acosd
   public :: atand, acotd, atan2d, secd, cscd, asecd, acscd
   public :: cot, csc, sec, acot, acsc, asec, coth, csch, sech, acoth
   public :: acsch, asech, sigmoid, logit, hypot_pracma, eps, real_part, imag_part, angle
   public :: sort_real

   interface Diag
      module procedure diag_from_vector
      module procedure diagonal_from_matrix
   end interface Diag
   interface zeros
      module procedure zeros_vector
      module procedure zeros_matrix
   end interface zeros
   interface ones
      module procedure ones_vector
      module procedure ones_matrix
   end interface ones
   interface Norm
      module procedure norm_vector
      module procedure norm_matrix
   end interface Norm
   interface std
      module procedure std_vector
   end interface std
   interface ceil
      module procedure ceil_real
   end interface ceil
   interface Fix
      module procedure fix_real
   end interface Fix
   interface mod
      module procedure mod_real
   end interface mod
   interface rem
      module procedure rem_real
   end interface rem
   interface pow2
      module procedure pow2_scalar
      module procedure pow2_pair
   end interface pow2
   interface eps
      module procedure eps_scalar
   end interface eps

contains

   function linspace(a, b, n) result(x)
      real(dp), intent(in) :: a, b
      integer, intent(in), optional :: n
      real(dp), allocatable :: x(:)
      integer :: m, i
      m = 100
      if (present(n)) m = n
      if (m <= 0) then
         allocate(x(0))
      else if (m == 1) then
         allocate(x(1)); x = b
      else
         allocate(x(m))
         do i = 1, m
            x(i) = a + real(i - 1, dp)*(b - a)/real(m - 1, dp)
         end do
         x(m) = b
      end if
   end function linspace

   function logspace(a, b, n, base) result(x)
      real(dp), intent(in) :: a, b
      integer, intent(in), optional :: n
      real(dp), intent(in), optional :: base
      real(dp), allocatable :: x(:), e(:)
      real(dp) :: radix
      radix = 10.0_dp
      if (present(base)) radix = base
      if (present(n)) then
         e = linspace(a, b, n)
      else
         e = linspace(a, b)
      end if
      allocate(x(size(e))); x = radix**e
   end function logspace

   function logseq(a, b, n, base) result(x)
      real(dp), intent(in) :: a, b
      integer, intent(in), optional :: n
      real(dp), intent(in), optional :: base
      real(dp), allocatable :: x(:)
      real(dp) :: radix
      integer :: m
      radix = exp(1.0_dp)
      if (present(base)) radix = base
      m = 50
      if (present(n)) m = n
      if (a <= 0.0_dp .or. b <= 0.0_dp .or. radix <= 0.0_dp .or. &
          abs(radix - 1.0_dp) <= eps_dp) then
         allocate(x(0)); return
      end if
      x = logspace(log(a)/log(radix), log(b)/log(radix), m, radix)
   end function logseq

   function eye(n, m) result(a)
      integer, intent(in) :: n
      integer, intent(in), optional :: m
      real(dp), allocatable :: a(:, :)
      integer :: j, nc
      nc = n
      if (present(m)) nc = m
      allocate(a(max(0,n), max(0,nc))); a = 0.0_dp
      do j = 1, min(n, nc)
         a(j,j) = 1.0_dp
      end do
   end function eye

   function zeros_vector(n) result(x)
      integer, intent(in) :: n
      real(dp), allocatable :: x(:)
      allocate(x(max(0,n))); x = 0.0_dp
   end function zeros_vector

   function zeros_matrix(n, m) result(a)
      integer, intent(in) :: n, m
      real(dp), allocatable :: a(:, :)
      allocate(a(max(0,n),max(0,m))); a = 0.0_dp
   end function zeros_matrix

   function ones_vector(n) result(x)
      integer, intent(in) :: n
      real(dp), allocatable :: x(:)
      allocate(x(max(0,n))); x = 1.0_dp
   end function ones_vector

   function ones_matrix(n, m) result(a)
      integer, intent(in) :: n, m
      real(dp), allocatable :: a(:, :)
      allocate(a(max(0,n),max(0,m))); a = 1.0_dp
   end function ones_matrix

   function diag_from_vector(x, k) result(a)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: k
      real(dp), allocatable :: a(:, :)
      integer :: d, i, r, c
      d = 0
      if (present(k)) d = k
      allocate(a(size(x)+abs(d), size(x)+abs(d))); a = 0.0_dp
      do i = 1, size(x)
         if (d >= 0) then
            r = i; c = i+d
         else
            r = i-d; c = i
         end if
         a(r,c) = x(i)
      end do
   end function diag_from_vector

   function diagonal_from_matrix(a, k) result(x)
      real(dp), intent(in) :: a(:, :)
      integer, intent(in), optional :: k
      real(dp), allocatable :: x(:)
      integer :: d, n, i, r, c
      d = 0
      if (present(k)) d = k
      if (d >= 0) then
         n = min(size(a,1), size(a,2)-d)
      else
         n = min(size(a,1)+d, size(a,2))
      end if
      n = max(0,n); allocate(x(n))
      do i = 1, n
         if (d >= 0) then
            r=i; c=i+d
         else
            r=i-d; c=i
         end if
         x(i)=a(r,c)
      end do
   end function diagonal_from_matrix

   function blkdiag(a, b) result(c)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), allocatable :: c(:, :)
      integer :: m1,n1,m2,n2
      m1=size(a,1); n1=size(a,2); m2=size(b,1); n2=size(b,2)
      allocate(c(m1+m2,n1+n2)); c=0.0_dp
      c(1:m1,1:n1)=a; c(m1+1:,n1+1:)=b
   end function blkdiag

   function hilb(n) result(a)
      integer, intent(in) :: n
      real(dp), allocatable :: a(:, :)
      integer :: i,j
      allocate(a(max(0,n),max(0,n)))
      do j=1,n
         do i=1,n
            a(i,j)=1.0_dp/real(i+j-1,dp)
         end do
      end do
   end function hilb

   function pascal(n, kind_value) result(a)
      integer, intent(in) :: n
      integer, intent(in), optional :: kind_value
      real(dp), allocatable :: a(:, :)
      integer :: i,j,knd
      knd=0
      if (present(kind_value)) knd=kind_value
      allocate(a(max(0,n),max(0,n))); a=0.0_dp
      if (n<=0) return
      do i=1,n
         a(i,1)=1.0_dp; a(1,i)=1.0_dp
      end do
      do i=2,n
         do j=2,n
            a(i,j)=a(i-1,j)+a(i,j-1)
         end do
      end do
      if (knd==1) then
         do j=1,n
            if (modulo(j,2)==0) a(:,j)=-a(:,j)
         end do
      else if (knd==2) then
         a=transpose(a)
      end if
   end function pascal

   function vander(x, increasing) result(a)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: increasing
      real(dp), allocatable :: a(:, :)
      logical :: inc
      integer :: n,i,j,p
      n=size(x); inc=.false.
      if (present(increasing)) inc=increasing
      allocate(a(n,n))
      do i=1,n
         do j=1,n
            if (inc) then
               p=j-1
            else
               p=n-j
            end if
            a(i,j)=x(i)**p
         end do
      end do
   end function vander

   function hankel(c, r) result(a)
      real(dp), intent(in) :: c(:)
      real(dp), intent(in), optional :: r(:)
      real(dp), allocatable :: a(:, :), rr(:)
      integer :: m,n,i,j,k
      m=size(c)
      if (present(r)) then
         n=size(r); allocate(rr(n)); rr=r
      else
         n=m; allocate(rr(n)); rr=0.0_dp
         if(n>0) rr(1)=c(m)
      end if
      allocate(a(m,n))
      do i=1,m
         do j=1,n
            k=i+j-1
            if (k<=m) then
               a(i,j)=c(k)
            else
               a(i,j)=rr(k-m+1)
            end if
         end do
      end do
   end function hankel

   function Toeplitz(c, r) result(a)
      real(dp), intent(in) :: c(:)
      real(dp), intent(in), optional :: r(:)
      real(dp), allocatable :: a(:, :), rr(:)
      integer :: m,n,i,j
      m=size(c)
      if (present(r)) then
         n=size(r); allocate(rr(n)); rr=r
      else
         n=m; allocate(rr(n)); rr=c
      end if
      if (m>0 .and. n>0) rr(1)=c(1)
      allocate(a(m,n))
      do i=1,m
         do j=1,n
            if (i>=j) then
               a(i,j)=c(i-j+1)
            else
               a(i,j)=rr(j-i+1)
            end if
         end do
      end do
   end function Toeplitz

   function flipud(a) result(b)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable :: b(:, :)
      integer :: i,m
      m=size(a,1); allocate(b(m,size(a,2)))
      do i=1,m
         b(i,:)=a(m-i+1,:)
      end do
   end function flipud

   function fliplr(a) result(b)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable :: b(:, :)
      integer :: j,n
      n=size(a,2); allocate(b(size(a,1),n))
      do j=1,n
         b(:,j)=a(:,n-j+1)
      end do
   end function fliplr

   function rot90(a, k) result(b)
      real(dp), intent(in) :: a(:, :)
      integer, intent(in), optional :: k
      real(dp), allocatable :: b(:, :), t(:, :)
      integer :: q,iter
      q=1
      if(present(k)) q=modulo(k,4)
      allocate(t(size(a,1),size(a,2))); t=a
      do iter=1,q
         b=flipud(transpose(t)); call move_alloc(b,t)
      end do
      call move_alloc(t,b)
   end function rot90

   function circshift(x, k) result(y)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: k
      real(dp), allocatable :: y(:)
      integer :: n,s,i
      n=size(x); allocate(y(n))
      if(n==0) return
      s=modulo(k,n)
      do i=1,n
         y(modulo(i-1+s,n)+1)=x(i)
      end do
   end function circshift

   function repmat(a, nr, nc) result(b)
      real(dp), intent(in) :: a(:, :)
      integer, intent(in) :: nr, nc
      real(dp), allocatable :: b(:, :)
      integer :: i,j,m,n
      m=size(a,1); n=size(a,2)
      allocate(b(max(0,nr)*m,max(0,nc)*n))
      do i=0,nr-1
         do j=0,nc-1
            b(i*m+1:(i+1)*m,j*n+1:(j+1)*n)=a
         end do
      end do
   end function repmat

   subroutine meshgrid(x, y, xx, yy)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), allocatable, intent(out) :: xx(:, :), yy(:, :)
      integer :: i,j
      allocate(xx(size(y),size(x)),yy(size(y),size(x)))
      do i=1,size(y)
         do j=1,size(x)
            xx(i,j)=x(j); yy(i,j)=y(i)
         end do
      end do
   end subroutine meshgrid

   pure function size2(a) result(s)
      real(dp), intent(in) :: a(:, :)
      integer :: s(2)
      s=[size(a,1),size(a,2)]
   end function size2

   pure integer function numel(a)
      real(dp), intent(in) :: a(:, :)
      numel=size(a)
   end function numel

   pure integer function nnz(a, tolerance)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in), optional :: tolerance
      real(dp) :: tol
      tol=0.0_dp
      if(present(tolerance)) tol=tolerance
      nnz=count(abs(a)>tol)
   end function nnz

   pure logical function isempty(a)
      real(dp), intent(in) :: a(:)
      isempty=size(a)==0
   end function isempty

   pure real(dp) function dot(x,y)
      real(dp), intent(in) :: x(:),y(:)
      if(size(x)/=size(y)) then
         dot=ieee_value(0.0_dp,ieee_quiet_nan)
      else
         dot=dot_product(x,y)
      end if
   end function dot

   pure function cross(x,y) result(z)
      real(dp), intent(in) :: x(3),y(3)
      real(dp) :: z(3)
      z=[x(2)*y(3)-x(3)*y(2),x(3)*y(1)-x(1)*y(3),x(1)*y(2)-x(2)*y(1)]
   end function cross

   function crossn(x) result(z)
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable :: z(:)
      integer :: n,j
      real(dp), allocatable :: minor(:, :)
      n=size(x,1); allocate(z(n)); z=0.0_dp
      if(size(x,2)/=n-1 .or. n<2) return
      do j=1,n
         minor=remove_row(x,j)
         if (modulo(j,2)==1) then
            z(j)=det_small(minor)
         else
            z(j)=-det_small(minor)
         end if
      end do
   end function crossn

   function remove_row(a,r) result(b)
      real(dp),intent(in)::a(:,:)
      integer,intent(in)::r
      real(dp),allocatable::b(:,:)
      integer::m
      m=size(a,1); allocate(b(m-1,size(a,2)))
      if(r>1)b(1:r-1,:)=a(1:r-1,:)
      if(r<m)b(r:,:)=a(r+1:m,:)
   end function remove_row

   function remove_col(a,c) result(b)
      real(dp),intent(in)::a(:,:)
      integer,intent(in)::c
      real(dp),allocatable::b(:,:)
      integer::n
      n=size(a,2); allocate(b(size(a,1),n-1))
      if(c>1)b(:,1:c-1)=a(:,1:c-1)
      if(c<n)b(:,c:)=a(:,c+1:n)
   end function remove_col

   recursive function det_small(a) result(d)
      real(dp),intent(in)::a(:,:)
      real(dp)::d
      real(dp),allocatable::m(:,:)
      integer::j,n
      n=size(a,1)
      if(n/=size(a,2))then
         d=0.0_dp; return
      end if
      if(n==0)then
         d=1.0_dp; return
      else if(n==1)then
         d=a(1,1); return
      else if(n==2)then
         d=a(1,1)*a(2,2)-a(1,2)*a(2,1); return
      end if
      d=0.0_dp
      do j=1,n
         m=remove_col(remove_row(a,1),j)
         if (modulo(j,2)==1) then
            d=d+a(1,j)*det_small(m)
         else
            d=d-a(1,j)*det_small(m)
         end if
      end do
   end function det_small

   pure real(dp) function norm_vector(x,p)
      real(dp),intent(in)::x(:)
      real(dp),intent(in),optional::p
      real(dp)::q
      q=2.0_dp
      if(present(p))q=p
      if(size(x)==0)then
         norm_vector=0.0_dp
      else if(q<0.0_dp)then
         norm_vector=minval(abs(x))
      else if(q==0.0_dp)then
         norm_vector=real(count(x/=0.0_dp),dp)
      else if(q>=huge_dp/2.0_dp)then
         norm_vector=maxval(abs(x))
      else
         norm_vector=sum(abs(x)**q)**(1.0_dp/q)
      end if
   end function norm_vector

   pure real(dp) function norm_matrix(a,p)
      real(dp),intent(in)::a(:,:)
      character(len=*),intent(in),optional::p
      if(.not.present(p))then
         norm_matrix=sqrt(sum(a*a))
      else
         select case(trim(p))
         case('F','f','fro')
            norm_matrix=sqrt(sum(a*a))
         case('1')
            norm_matrix=maxval(sum(abs(a),dim=1))
         case('I','i','inf')
            norm_matrix=maxval(sum(abs(a),dim=2))
         case default
            norm_matrix=sqrt(sum(a*a))
         end select
      end if
   end function norm_matrix

   pure real(dp) function fnorm(a)
      real(dp),intent(in)::a(:,:)
      fnorm=sqrt(sum(a*a))
   end function fnorm

   pure real(dp) function Trace(a)
      real(dp),intent(in)::a(:,:)
      integer::i
      Trace=0.0_dp
      do i=1,min(size(a,1),size(a,2))
         Trace=Trace+a(i,i)
      end do
   end function Trace

   pure real(dp) function harmmean(x)
      real(dp),intent(in)::x(:)
      if(size(x)==0 .or. any(x==0.0_dp))then
         harmmean=ieee_value(0.0_dp,ieee_quiet_nan)
      else
         harmmean=real(size(x),dp)/sum(1.0_dp/x)
      end if
   end function harmmean

   pure real(dp) function geomean(x)
      real(dp),intent(in)::x(:)
      if(size(x)==0 .or. any(x<0.0_dp))then
         geomean=ieee_value(0.0_dp,ieee_quiet_nan)
      else if(any(x==0.0_dp))then
         geomean=0.0_dp
      else
         geomean=exp(sum(log(x))/real(size(x),dp))
      end if
   end function geomean

   function trimmean(x, percent) result(m)
      real(dp),intent(in)::x(:),percent
      real(dp)::m
      real(dp),allocatable::y(:)
      integer::k,n
      y=sort_real(x); n=size(y)
      k=int(floor(0.5_dp*percent/100.0_dp*real(n,dp)))
      if(n-2*k<=0)then
         m=ieee_value(0.0_dp,ieee_quiet_nan)
      else
         m=sum(y(k+1:n-k))/real(n-2*k,dp)
      end if
   end function trimmean

   function Mode(x) result(v)
      real(dp),intent(in)::x(:)
      real(dp)::v
      real(dp),allocatable::y(:)
      integer::i,best,count_now,best_count
      if(size(x)==0)then
         v=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      y=sort_real(x); best=1; count_now=1; best_count=1
      do i=2,size(y)
         if(y(i)==y(i-1))then
            count_now=count_now+1
         else
            count_now=1
         end if
         if(count_now>best_count)then
            best_count=count_now; best=i
         end if
      end do
      v=y(best)
   end function Mode

   pure real(dp) function std_vector(x, population)
      real(dp),intent(in)::x(:)
      logical,intent(in),optional::population
      logical::pop
      real(dp)::m
      integer::den
      pop=.false.
      if(present(population))pop=population
      if(size(x)==0 .or. (size(x)==1 .and. .not.pop))then
         std_vector=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      m=sum(x)/real(size(x),dp)
      if (pop) then
         den=size(x)
      else
         den=size(x)-1
      end if
      std_vector=sqrt(sum((x-m)**2)/real(den,dp))
   end function std_vector

   pure real(dp) function std_err(x)
      real(dp),intent(in)::x(:)
      if(size(x)==0)then
         std_err=ieee_value(0.0_dp,ieee_quiet_nan)
      else
         std_err=std_vector(x)/sqrt(real(size(x),dp))
      end if
   end function std_err

   elemental real(dp) function ceil_real(x)
      real(dp),intent(in)::x
      ceil_real=real(ceiling(x),dp)
   end function ceil_real

   elemental real(dp) function fix_real(x)
      real(dp),intent(in)::x
      fix_real=real(int(x),dp)
   end function fix_real

   elemental real(dp) function mod_real(x,y)
      real(dp),intent(in)::x,y
      if(y==0.0_dp)then
         mod_real=ieee_value(0.0_dp,ieee_quiet_nan)
      else
         mod_real=modulo(x,y)
      end if
   end function mod_real

   elemental real(dp) function rem_real(x,y)
      real(dp),intent(in)::x,y
      if(y==0.0_dp)then
         rem_real=ieee_value(0.0_dp,ieee_quiet_nan)
      else
         rem_real=x-real(int(x/y),dp)*y
      end if
   end function rem_real

   elemental integer function idivide(a,b,mode)
      integer,intent(in)::a,b
      character(len=*),intent(in),optional::mode
      real(dp)::q
      if(b==0)then
         idivide=0; return
      end if
      q=real(a,dp)/real(b,dp)
      if(.not.present(mode))then
         idivide=int(q)
      else
         select case(trim(mode))
         case('fix'); idivide=int(q)
         case('floor'); idivide=floor(q)
         case('ceil'); idivide=ceiling(q)
         case('round'); idivide=nint(q)
         case default; idivide=int(q)
         end select
      end if
   end function idivide

   pure integer function gcd(a,b)
      integer,intent(in)::a,b
      integer::x,y,t
      x=abs(a); y=abs(b)
      do while(y/=0)
         t=modulo(x,y); x=y; y=t
      end do
      gcd=x
   end function gcd

   pure integer function Lcm(a,b)
      integer,intent(in)::a,b
      if(a==0 .or. b==0)then
         Lcm=0
      else
         Lcm=abs((a/gcd(a,b))*b)
      end if
   end function Lcm

   elemental integer function nextpow2(x)
      real(dp),intent(in)::x
      if(x<=0.0_dp)then
         nextpow2=0
      else
         nextpow2=ceiling(log(x)/log(2.0_dp)-8.0_dp*eps_dp)
      end if
   end function nextpow2

   elemental real(dp) function pow2_scalar(x)
      real(dp),intent(in)::x
      pow2_scalar=2.0_dp**x
   end function pow2_scalar

   elemental real(dp) function pow2_pair(f,e)
      real(dp),intent(in)::f
      integer,intent(in)::e
      pow2_pair=f*2.0_dp**e
   end function pow2_pair

   function distmat(x,y,p) result(d)
      real(dp),intent(in)::x(:,:),y(:,:)
      real(dp),intent(in),optional::p
      real(dp),allocatable::d(:,:)
      real(dp)::q
      integer::i,j
      q=2.0_dp
      if(present(p))q=p
      allocate(d(size(x,1),size(y,1)))
      if(size(x,2)/=size(y,2))then
         d=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      do i=1,size(x,1)
         do j=1,size(y,1)
            if(q>=huge_dp/2.0_dp)then
               d(i,j)=maxval(abs(x(i,:)-y(j,:)))
            else
               d(i,j)=sum(abs(x(i,:)-y(j,:))**q)**(1.0_dp/q)
            end if
         end do
      end do
   end function distmat

   function pdist(x,p) result(d)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(in),optional::p
      real(dp),allocatable::d(:)
      real(dp)::q
      integer::i,j,k,n
      q=2.0_dp
      if(present(p))q=p
      n=size(x,1); allocate(d(n*(n-1)/2)); k=0
      do j=1,n-1
         do i=j+1,n
            k=k+1; d(k)=sum(abs(x(i,:)-x(j,:))**q)**(1.0_dp/q)
         end do
      end do
   end function pdist

   function pdist2(x,y,p) result(d)
      real(dp),intent(in)::x(:,:),y(:,:)
      real(dp),intent(in),optional::p
      real(dp),allocatable::d(:,:)
      if(present(p))then
         d=distmat(x,y,p)
      else
         d=distmat(x,y)
      end if
   end function pdist2

   function hausdorff_dist(x,y,p) result(h)
      real(dp),intent(in)::x(:,:),y(:,:)
      real(dp),intent(in),optional::p
      real(dp)::h
      real(dp),allocatable::d(:,:)
      if(present(p))then
         d=distmat(x,y,p)
      else
         d=distmat(x,y)
      end if
      h=max(maxval(minval(d,dim=2)),maxval(minval(d,dim=1)))
   end function hausdorff_dist

   function accumarray(indices,values,n,fill_value) result(out)
      integer,intent(in)::indices(:)
      real(dp),intent(in)::values(:)
      integer,intent(in),optional::n
      real(dp),intent(in),optional::fill_value
      real(dp),allocatable::out(:)
      logical,allocatable::seen(:)
      real(dp)::fill
      integer::m,i
      fill=0.0_dp
      if(present(fill_value))fill=fill_value
      if(size(indices)/=size(values))then
         allocate(out(0)); return
      end if
      m=0
      if(size(indices)>0)m=maxval(indices)
      if(present(n))m=n
      allocate(out(max(0,m)),seen(max(0,m)))
      out=0.0_dp; seen=.false.
      do i=1,size(indices)
         if(indices(i)>=1 .and. indices(i)<=m)then
            out(indices(i))=out(indices(i))+values(i)
            seen(indices(i))=.true.
         end if
      end do
      where(.not.seen) out=fill
   end function accumarray

   function bsxfun(a,b,op) result(c)
      real(dp),intent(in)::a(:,:),b(:,:)
      character(len=*),intent(in)::op
      real(dp),allocatable::c(:,:)
      integer::m,n,i,j,ia,ib,ja,jb
      if(.not.(size(a,1)==size(b,1).or.size(a,1)==1.or.size(b,1)==1) .or. &
         .not.(size(a,2)==size(b,2).or.size(a,2)==1.or.size(b,2)==1))then
         allocate(c(0,0)); return
      end if
      m=max(size(a,1),size(b,1)); n=max(size(a,2),size(b,2)); allocate(c(m,n))
      do i=1,m
         do j=1,n
            if (size(a,1)==1) then; ia=1; else; ia=i; end if
            if (size(b,1)==1) then; ib=1; else; ib=i; end if
            if (size(a,2)==1) then; ja=1; else; ja=j; end if
            if (size(b,2)==1) then; jb=1; else; jb=j; end if
            select case(trim(op))
            case('+','plus'); c(i,j)=a(ia,ja)+b(ib,jb)
            case('-','minus'); c(i,j)=a(ia,ja)-b(ib,jb)
            case('*','times'); c(i,j)=a(ia,ja)*b(ib,jb)
            case('/','rdivide'); c(i,j)=a(ia,ja)/b(ib,jb)
            case('max'); c(i,j)=max(a(ia,ja),b(ib,jb))
            case('min'); c(i,j)=min(a(ia,ja),b(ib,jb))
            case default; c(i,j)=ieee_value(0.0_dp,ieee_quiet_nan)
            end select
         end do
      end do
   end function bsxfun

   function uniq(x,tolerance) result(y)
      real(dp),intent(in)::x(:)
      real(dp),intent(in),optional::tolerance
      real(dp),allocatable::y(:),tmp(:),s(:)
      real(dp)::tol
      integer::i,k
      tol=0.0_dp
      if(present(tolerance))tol=tolerance
      s=sort_real(x); allocate(tmp(size(s))); k=0
      do i=1,size(s)
         if(i==1)then
            k=k+1; tmp(k)=s(i)
         else if(abs(s(i)-s(i-1))>tol)then
            k=k+1; tmp(k)=s(i)
         end if
      end do
      allocate(y(k))
      if(k>0)y=tmp(1:k)
   end function uniq

   function find(mask) result(idx)
      logical,intent(in)::mask(:)
      integer,allocatable::idx(:)
      integer::i,k
      allocate(idx(count(mask))); k=0
      do i=1,size(mask)
         if(mask(i))then
            k=k+1; idx(k)=i
         end if
      end do
   end function find

   function finds(x,value,tolerance) result(idx)
      real(dp),intent(in)::x(:),value
      real(dp),intent(in),optional::tolerance
      integer,allocatable::idx(:)
      real(dp)::tol
      tol=0.0_dp
      if(present(tolerance))tol=tolerance
      idx=find(abs(x-value)<=tol)
   end function finds

   function findintervals(x,breaks) result(idx)
      real(dp),intent(in)::x(:),breaks(:)
      integer,allocatable::idx(:)
      integer::i,j
      allocate(idx(size(x)))
      do i=1,size(x)
         idx(i)=0
         do j=1,size(breaks)
            if(x(i)>=breaks(j))idx(i)=j
         end do
      end do
   end function findintervals

   elemental real(dp) function deg2rad(x)
      real(dp),intent(in)::x
      deg2rad=x*pi_dp/180.0_dp
   end function deg2rad
   elemental real(dp) function rad2deg(x)
      real(dp),intent(in)::x
      rad2deg=x*180.0_dp/pi_dp
   end function rad2deg
   elemental real(dp) function sind(x)
      real(dp),intent(in)::x
      sind=sin(deg2rad(x))
   end function sind
   elemental real(dp) function cosd(x)
      real(dp),intent(in)::x
      cosd=cos(deg2rad(x))
   end function cosd
   elemental real(dp) function tand(x)
      real(dp),intent(in)::x
      tand=tan(deg2rad(x))
   end function tand
   elemental real(dp) function cotd(x)
      real(dp),intent(in)::x
      cotd=1.0_dp/tand(x)
   end function cotd
   elemental real(dp) function asind(x)
      real(dp),intent(in)::x
      asind=rad2deg(asin(x))
   end function asind
   elemental real(dp) function acosd(x)
      real(dp),intent(in)::x
      acosd=rad2deg(acos(x))
   end function acosd
   elemental real(dp) function atand(x)
      real(dp),intent(in)::x
      atand=rad2deg(atan(x))
   end function atand
   elemental real(dp) function acotd(x)
      real(dp),intent(in)::x
      acotd=rad2deg(atan2(1.0_dp,x))
   end function acotd
   elemental real(dp) function atan2d(y,x)
      real(dp),intent(in)::y,x
      atan2d=rad2deg(atan2(y,x))
   end function atan2d
   elemental real(dp) function secd(x)
      real(dp),intent(in)::x
      secd=1.0_dp/cosd(x)
   end function secd
   elemental real(dp) function cscd(x)
      real(dp),intent(in)::x
      cscd=1.0_dp/sind(x)
   end function cscd
   elemental real(dp) function asecd(x)
      real(dp),intent(in)::x
      asecd=acosd(1.0_dp/x)
   end function asecd
   elemental real(dp) function acscd(x)
      real(dp),intent(in)::x
      acscd=asind(1.0_dp/x)
   end function acscd
   elemental real(dp) function cot(x)
      real(dp),intent(in)::x
      cot=1.0_dp/tan(x)
   end function cot
   elemental real(dp) function csc(x)
      real(dp),intent(in)::x
      csc=1.0_dp/sin(x)
   end function csc
   elemental real(dp) function sec(x)
      real(dp),intent(in)::x
      sec=1.0_dp/cos(x)
   end function sec
   elemental real(dp) function acot(x)
      real(dp),intent(in)::x
      acot=atan2(1.0_dp,x)
   end function acot
   elemental real(dp) function acsc(x)
      real(dp),intent(in)::x
      acsc=asin(1.0_dp/x)
   end function acsc
   elemental real(dp) function asec(x)
      real(dp),intent(in)::x
      asec=acos(1.0_dp/x)
   end function asec
   elemental real(dp) function coth(x)
      real(dp),intent(in)::x
      coth=1.0_dp/tanh(x)
   end function coth
   elemental real(dp) function csch(x)
      real(dp),intent(in)::x
      csch=1.0_dp/sinh(x)
   end function csch
   elemental real(dp) function sech(x)
      real(dp),intent(in)::x
      sech=1.0_dp/cosh(x)
   end function sech
   elemental real(dp) function acoth(x)
      real(dp),intent(in)::x
      acoth=0.5_dp*log((x+1.0_dp)/(x-1.0_dp))
   end function acoth
   elemental real(dp) function acsch(x)
      real(dp),intent(in)::x
      acsch=asinh(1.0_dp/x)
   end function acsch
   elemental real(dp) function asech(x)
      real(dp),intent(in)::x
      asech=acosh(1.0_dp/x)
   end function asech
   elemental real(dp) function sigmoid(x)
      real(dp),intent(in)::x
      sigmoid=1.0_dp/(1.0_dp+exp(-x))
   end function sigmoid
   elemental real(dp) function logit(x)
      real(dp),intent(in)::x
      logit=log(x/(1.0_dp-x))
   end function logit
   elemental real(dp) function hypot_pracma(x,y)
      real(dp),intent(in)::x,y
      real(dp)::s,a,b
      a=abs(x); b=abs(y); s=max(a,b)
      if(s==0.0_dp)then
         hypot_pracma=0.0_dp
      else
         hypot_pracma=s*sqrt((a/s)**2+(b/s)**2)
      end if
   end function hypot_pracma
   elemental real(dp) function eps_scalar(x)
      real(dp),intent(in)::x
      if(x==0.0_dp)then
         eps_scalar=tiny(1.0_dp)
      else
         eps_scalar=spacing(x)
      end if
   end function eps_scalar
   elemental real(dp) function real_part(z)
      complex(dp),intent(in)::z
      real_part=real(z,dp)
   end function real_part
   elemental real(dp) function imag_part(z)
      complex(dp),intent(in)::z
      imag_part=aimag(z)
   end function imag_part
   elemental real(dp) function angle(z)
      complex(dp),intent(in)::z
      angle=atan2(aimag(z),real(z,dp))
   end function angle

   function sort_real(x) result(y)
      real(dp),intent(in)::x(:)
      real(dp),allocatable::y(:)
      integer::i,j
      real(dp)::v
      allocate(y(size(x))); y=x
      do i=2,size(y)
         v=y(i); j=i-1
         do while(j>=1)
            if(y(j)<=v)exit
            y(j+1)=y(j); j=j-1
         end do
         y(j+1)=v
      end do
   end function sort_real

end module pracma_basic
