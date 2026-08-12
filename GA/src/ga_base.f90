! GA-fortran: translation of computational algorithms from R package GA.
! Copyright (C) 2012-2026 Luca Scrucca (original GA package)
! Fortran translation: 2026
! License: GPL-2.0-or-later. See LICENSE and original/.

module ga_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
end module ga_kinds

module ga_random
  use ga_kinds, only : dp
  implicit none
  private
  public :: ga_seed, runif, randint, sample_without_replacement, shuffle_int
contains
  subroutine ga_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    integer(kind=8) :: z
    call random_seed(size=n)
    allocate(put(n))
    z = int(seed,8)
    if (z == 0_8) z = 104729_8
    do i = 1, n
      z = modulo(1103515245_8*z + 12345_8 + 7919_8*i, 2147483647_8)
      if (z <= 0_8) z = z + 2147483646_8
      put(i) = int(z)
    end do
    call random_seed(put=put)
  end subroutine ga_seed

  real(dp) function runif(a,b) result(x)
    real(dp), intent(in), optional :: a,b
    real(dp) :: u, lo, hi
    call random_number(u)
    lo = 0.0_dp; hi = 1.0_dp
    if (present(a)) lo = a
    if (present(b)) hi = b
    x = lo + (hi-lo)*u
  end function runif

  integer function randint(lo,hi) result(k)
    integer, intent(in) :: lo, hi
    real(dp) :: u
    if (hi < lo) error stop "randint: invalid bounds"
    call random_number(u)
    k = lo + int(u*real(hi-lo+1,dp))
    if (k > hi) k = hi
  end function randint

  subroutine shuffle_int(x)
    integer, intent(inout) :: x(:)
    integer :: i,j,t
    do i = size(x), 2, -1
      j = randint(1,i)
      t=x(i); x(i)=x(j); x(j)=t
    end do
  end subroutine shuffle_int

  subroutine sample_without_replacement(n,k,out)
    integer, intent(in) :: n,k
    integer, intent(out) :: out(k)
    integer, allocatable :: a(:)
    integer :: i
    if (k < 0 .or. k > n) error stop "sample_without_replacement: invalid k"
    allocate(a(n))
    a = [(i,i=1,n)]
    call shuffle_int(a)
    out = a(1:k)
  end subroutine sample_without_replacement
end module ga_random

module ga_utils
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use ga_kinds, only : dp
  use ga_random, only : runif
  implicit none
  private
  public :: decimal2binary, binary2decimal, binary2gray, gray2binary
  public :: rank_decreasing, weighted_index, ga_pmutation, garun, optim_probsel
  public :: mean_sd, fitness_summary, unique_rows_real, unique_rows_int
  public :: repair_solution, reflect_solution
  public :: argsort_desc, argsort_asc, clamp
contains
  pure real(dp) function clamp(x,a,b) result(y)
    real(dp), intent(in) :: x,a,b
    y = min(max(x,a),b)
  end function clamp

  subroutine decimal2binary(x,b)
    integer, intent(in) :: x
    integer, intent(out) :: b(:)
    integer :: v,i
    v=x; b=0
    do i=size(b),1,-1
      b(i)=mod(v,2)
      v=v/2
    end do
  end subroutine decimal2binary

  integer function binary2decimal(b) result(x)
    integer, intent(in) :: b(:)
    integer :: i
    x=0
    do i=1,size(b)
      x=2*x+b(i)
    end do
  end function binary2decimal

  subroutine binary2gray(b,g)
    integer, intent(in) :: b(:)
    integer, intent(out) :: g(size(b))
    integer :: i
    if(size(b)==0) return
    g(1)=b(1)
    do i=2,size(b)
      g(i)=ieor(b(i-1),b(i))
    end do
  end subroutine binary2gray

  subroutine gray2binary(g,b)
    integer, intent(in) :: g(:)
    integer, intent(out) :: b(size(g))
    integer :: i
    if(size(g)==0) return
    b(1)=g(1)
    do i=2,size(g)
      b(i)=ieor(b(i-1),g(i))
    end do
  end subroutine gray2binary

  subroutine rank_decreasing(x,rank)
    real(dp), intent(in) :: x(:)
    integer, intent(out) :: rank(size(x))
    integer, allocatable :: idx(:)
    integer :: i,j,t,n
    n=size(x); allocate(idx(n)); idx=[(i,i=1,n)]
    ! Stable insertion sort, finite values first in decreasing order.
    do i=2,n
      t=idx(i); j=i-1
      do while(j>=1)
        if (.not. greater_idx(t,idx(j),x)) exit
        idx(j+1)=idx(j); j=j-1
      end do
      idx(j+1)=t
    end do
    do i=1,n
      rank(idx(i))=i
    end do
  contains
    logical function greater_idx(i,j,v)
      integer,intent(in)::i,j
      real(dp),intent(in)::v(:)
      if (.not.ieee_is_finite(v(i))) then
        greater_idx=.false.
      else if (.not.ieee_is_finite(v(j))) then
        greater_idx=.true.
      else if (v(i)>v(j)) then
        greater_idx=.true.
      else if (v(i)<v(j)) then
        greater_idx=.false.
      else
        greater_idx=i<j
      end if
    end function greater_idx
  end subroutine rank_decreasing

  integer function weighted_index(prob) result(idx)
    real(dp), intent(in) :: prob(:)
    real(dp) :: s,u,c
    integer :: i,n
    n=size(prob); s=sum(max(prob,0.0_dp))
    if (s<=0.0_dp .or. .not.ieee_is_finite(s)) then
      idx=1+int(runif()*real(n,dp)); if(idx>n) idx=n; return
    end if
    u=runif()*s; c=0.0_dp
    do i=1,n
      c=c+max(prob(i),0.0_dp)
      if(u<=c) then; idx=i; return; end if
    end do
    idx=n
  end function weighted_index

  subroutine optim_probsel(x,prob,q)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::prob(size(x))
    real(dp),intent(in),optional::q
    integer,allocatable::r(:)
    real(dp)::qv,eps,s
    integer::i
    eps=sqrt(epsilon(1.0_dp));qv=0.25_dp;if(present(q))qv=q
    qv=min(max(eps,qv),1.0_dp-eps);allocate(r(size(x)));call rank_decreasing(x,r)
    do i=1,size(x)
      prob(i)=exp(log(qv)+real(r(i)-1,dp)*log(1.0_dp-qv))
      if(.not.ieee_is_finite(x(i)))prob(i)=0.0_dp
    end do
    s=sum(prob);if(s>0.0_dp)prob=prob/s
  end subroutine optim_probsel

  pure real(dp) function ga_pmutation(iter,maxiter,p0,p,t) result(pm)
    integer, intent(in) :: iter,maxiter
    real(dp), intent(in), optional :: p0,p,t
    real(dp) :: p0v,pv,tv
    p0v=0.5_dp; pv=0.01_dp; tv=real(nint(real(maxiter,dp)/2.0_dp),dp)
    if(present(p0)) p0v=p0
    if(present(p)) pv=p
    if(present(t)) tv=t
    if(tv<=0.0_dp) then
      pm=pv
    else
      pm=(p0v-pv)*exp(-2.0_dp*real(iter-1,dp)/tv)+pv
    end if
  end function ga_pmutation

  integer function garun(best,eps) result(nrun)
    real(dp), intent(in) :: best(:)
    real(dp), intent(in), optional :: eps
    real(dp) :: e,m
    e=sqrt(epsilon(1.0_dp)); if(present(eps)) e=eps
    if(size(best)==0) then; nrun=0; return; end if
    m=maxval(best)
    nrun=count(best>=m-e)
  end function garun

  subroutine mean_sd(x,mean,sd)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::mean,sd
    integer::n
    n=size(x)
    if(n==0) then; mean=0; sd=0; return; end if
    mean=sum(x)/real(n,dp)
    if(n>1) then
      sd=sqrt(sum((x-mean)**2)/real(n-1,dp))
    else
      sd=0.0_dp
    end if
  end subroutine mean_sd

  subroutine fitness_summary(f,s)
    real(dp),intent(in)::f(:)
    real(dp),intent(out)::s(6)
    real(dp),allocatable::a(:)
    real(dp)::m,d(5),n4
    integer::i,n,k
    n=count([(ieee_is_finite(f(i)),i=1,size(f))])
    if(n==0) then
      s=0.0_dp
      return
    end if
    allocate(a(n)); k=0
    do i=1,size(f)
      if(ieee_is_finite(f(i))) then
        k=k+1
        a(k)=f(i)
      end if
    end do
    call sort_real(a)
    m=sum(a)/real(n,dp)
    n4=real(floor(real(n+3,dp)/2.0_dp),dp)/2.0_dp
    d=[1.0_dp,n4,real(n+1,dp)/2.0_dp,real(n+1,dp)-n4,real(n,dp)]
    s(1)=a(n)
    s(2)=m
    s(3)=hinge(a,d(4))
    s(4)=hinge(a,d(3))
    s(5)=hinge(a,d(2))
    s(6)=a(1)
  contains
    subroutine sort_real(v)
      real(dp),intent(inout)::v(:)
      integer::ii,jj
      real(dp)::t
      do ii=2,size(v)
        t=v(ii); jj=ii-1
        do while(jj>=1)
          if(v(jj)<=t) exit
          v(jj+1)=v(jj)
          jj=jj-1
        end do
        v(jj+1)=t
      end do
    end subroutine sort_real
    real(dp) function hinge(v,z)
      real(dp),intent(in)::v(:),z
      integer::il,ih
      il=max(1,int(floor(z))); ih=min(size(v),int(ceiling(z)))
      hinge=0.5_dp*(v(il)+v(ih))
    end function hinge
  end subroutine fitness_summary

  subroutine repair_solution(x,lo,up)
    real(dp),intent(inout)::x(:)
    real(dp),intent(in)::lo(:),up(:)
    real(dp),allocatable::xl(:),xu(:)
    allocate(xl(size(x)),xu(size(x)))
    xl=lo-x; xl=xl+abs(xl)
    xu=x-up; xu=xu+abs(xu)
    x=x-(xu-xl)/2.0_dp
  end subroutine repair_solution

  subroutine reflect_solution(x,lo,up,tol)
    real(dp),intent(inout)::x(:)
    real(dp),intent(in)::lo(:),up(:)
    real(dp),intent(in),optional::tol
    real(dp),allocatable::xl(:),xu(:),r(:)
    real(dp)::e,t
    t=sqrt(epsilon(1.0_dp)); if(present(tol))t=tol
    allocate(xl(size(x)),xu(size(x)),r(size(x))); r=up-lo
    do
      e=sum(x-up+abs(x-up)+lo-x+abs(lo-x))
      if(e<=t)exit
      xu=x-up; xu=xu+abs(xu); xu=xu+r-abs(xu-r)
      xl=lo-x; xl=xl+abs(xl); xl=xl+r-abs(xl-r)
      x=x-(xu-xl)/2.0_dp
    end do
  end subroutine reflect_solution

  subroutine argsort_desc(x,idx)
    real(dp),intent(in)::x(:); integer,intent(out)::idx(size(x)); integer::i,j,t
    idx=[(i,i=1,size(x))]
    do i=2,size(idx)
    t=idx(i)
    j=i-1
    do while(j>=1)
      if(x(idx(j))>=x(t)) exit
      idx(j+1)=idx(j)
      j=j-1
    end do
    idx(j+1)=t
    end do
  end subroutine argsort_desc

  subroutine argsort_asc(x,idx)
    real(dp),intent(in)::x(:); integer,intent(out)::idx(size(x)); integer::i,j,t
    idx=[(i,i=1,size(x))]
    do i=2,size(idx)
    t=idx(i)
    j=i-1
    do while(j>=1)
      if(x(idx(j))<=x(t)) exit
      idx(j+1)=idx(j)
      j=j-1
    end do
    idx(j+1)=t
    end do
  end subroutine argsort_asc

  integer function unique_rows_real(a,tol) result(nu)
    real(dp),intent(in)::a(:,:); real(dp),intent(in),optional::tol
    real(dp)::t; integer::i,j; logical::dup
    t=sqrt(epsilon(1.0_dp)); if(present(tol)) t=tol
    nu=0
    do i=1,size(a,1)
    dup=.false.
    do j=1,i-1
    if(maxval(abs(a(i,:)-a(j,:)))<=t) then
    dup=.true.
    exit
    end if
    end do
    if(.not.dup) nu=nu+1
    end do
  end function
  integer function unique_rows_int(a) result(nu)
    integer,intent(in)::a(:,:); integer::i,j; logical::dup
    nu=0
    do i=1,size(a,1)
    dup=.false.
    do j=1,i-1
    if(all(a(i,:)==a(j,:))) then
    dup=.true.
    exit
    end if
    end do
    if(.not.dup) nu=nu+1
    end do
  end function
end module ga_utils
