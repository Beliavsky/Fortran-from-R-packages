! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_compat
   use pracma_kinds, only : dp, i8, pi_dp
   use pracma_status, only : pracma_ok, pracma_invalid_argument, pracma_singular
   use pracma_callbacks, only : complex_scalar_function, vector_field, scalar_function
   use pracma_types, only : ode_result, optimization_result, pchip_result
   use pracma_linalg, only : pinv, nullspace
   use pracma_interpolation, only : akima, ppval, cubicspline, lebesgue_function
   use pracma_geometry, only : segment_distance, segment_intersection, polygon_crossings, &
      stereographic_project, stereographic_inverse
   use pracma_signal_stats, only : hurst
   use pracma_ode, only : abm3
   use pracma_optimization, only : fminbnd
   use pracma_combinatorics, only : rand_uniform, randn, randperm, dec2bin
   implicit none
   private
   public :: invlap, rand, rand_matrix, randn_matrix, rands, randp, randsample, bits
   public :: akimaInterp, ppfit, lsqlin, lsqlincon, kron, trisolve
   public :: segm_intersect, segm_distance, poly_crossings, tril, triu, squareform, flipdim
   public :: softline, piecewise, divisors, rssimple, hurstexp, approx_entropy, sample_entropy
   public :: bubbleSort, insertionSort, selectionSort, shellSort, heapSort, mergeSort
   public :: quickSort, quickSortx, mergeOrdered, is_sorted, sortrows
   public :: normF, abm3pc, shubert, bvp, stereographic, stereographic_inv, lebesgue, half
   interface squareform
      module procedure squareform_vector
      module procedure squareform_matrix
   end interface squareform


contains

   subroutine invlap(fs,t1,t2,nnt,t,values,a,ns,nd,status)
      procedure(complex_scalar_function)::fs
      real(dp),intent(in)::t1,t2
      integer,intent(in)::nnt
      real(dp),allocatable,intent(out)::t(:),values(:)
      real(dp),intent(in),optional::a
      integer,intent(in),optional::ns,nd
      integer,intent(out),optional::status
      real(dp)::aa,tt,bd
      complex(dp),allocatable::alpha(:),beta(:)
      integer::nseries,ndif,n,i,k,start
      aa=6.0_dp; if(present(a))aa=a; nseries=20; if(present(ns))nseries=ns
      ndif=19; if(present(nd))ndif=nd
      if(nnt<2.or.t2<=t1.or.nseries<1.or.ndif<0)then
         allocate(t(0),values(0)); if(present(status))status=pracma_invalid_argument; return
      end if
      t=[(t1+(t2-t1)*real(i-1,dp)/real(nnt-1,dp),i=1,nnt)]
      if(t1==0.0_dp)t=t(2:)
      allocate(values(size(t)),alpha(nseries+1+ndif),beta(nseries+1+ndif))
      do n=1,size(alpha)
         alpha(n)=cmplx(aa,real(n-1,dp)*pi_dp,dp)
         beta(n)=cmplx(-exp(aa)*merge(1.0_dp,-1.0_dp,mod(n,2)==0),0.0_dp,dp)
      end do
      if(ndif>0)then
         start=nseries+2
         do k=1,ndif
            bd=0.0_dp
            do i=k,ndif
               bd=bd+gamma(real(ndif+1,dp))/(gamma(real(ndif+2-i,dp))*gamma(real(i,dp)))
            end do
            beta(start+k-1)=beta(start+k-1)*bd/(2.0_dp**ndif)
         end do
      end if
      beta(1)=0.5_dp*beta(1)
      do i=1,size(t)
         tt=t(i); values(i)=sum(real((beta/tt)*[(fs(alpha(k)/tt),k=1,size(alpha))],dp))
      end do
      if(present(status))status=pracma_ok
   end subroutine invlap

   function rand_matrix(n,m) result(a)
      integer,intent(in)::n
      integer,intent(in),optional::m
      real(dp),allocatable::a(:,:)
      integer::i,j,nc
      nc=n; if(present(m))nc=m; allocate(a(max(0,n),max(0,nc)))
      do j=1,nc; do i=1,n; a(i,j)=rand_uniform(); end do; end do
   end function rand_matrix


   function rand(n,m) result(a)
      integer,intent(in)::n
      integer,intent(in),optional::m
      real(dp),allocatable::a(:,:)
      a=rand_matrix(n,m)
   end function rand

   function bits(n,width) result(b)
      integer(i8),intent(in)::n
      integer,intent(in),optional::width
      integer,allocatable::b(:)
      b=dec2bin(n,width)
   end function bits

   function randn_matrix(n,m) result(a)
      integer,intent(in)::n
      integer,intent(in),optional::m
      real(dp),allocatable::a(:,:)
      integer::i,j,nc
      nc=n; if(present(m))nc=m; allocate(a(max(0,n),max(0,nc)))
      do j=1,nc; do i=1,n; a(i,j)=randn(); end do; end do
   end function randn_matrix

   function rands(n,count,radius) result(points)
      integer,intent(in)::n
      integer,intent(in),optional::count
      real(dp),intent(in),optional::radius
      real(dp),allocatable::points(:,:)
      real(dp)::r,normv
      integer::i,j,k
      k=1; if(present(count))k=count; r=1.0_dp; if(present(radius))r=radius
      allocate(points(n,k))
      do j=1,k
         do i=1,n; points(i,j)=randn(); end do
         normv=sqrt(sum(points(:,j)**2)); if(normv>0)points(:,j)=r*points(:,j)/normv
      end do
   end function rands

   function randp(n,radius) result(points)
      integer,intent(in)::n
      real(dp),intent(in),optional::radius
      real(dp),allocatable::points(:,:)
      real(dp)::r,rho,theta
      integer::i
      r=1.0_dp; if(present(radius))r=radius; allocate(points(2,n))
      do i=1,n; rho=r*sqrt(rand_uniform()); theta=2*pi_dp*rand_uniform(); points(:,i)=[rho*cos(theta),rho*sin(theta)]; end do
   end function randp

   function randsample(values,k,replacement) result(sample)
      integer,intent(in)::values(:),k
      logical,intent(in),optional::replacement
      integer,allocatable::sample(:),p(:)
      logical::rep
      integer::i
      rep=.false.; if(present(replacement))rep=replacement; allocate(sample(k))
      if(rep.or.k>size(values))then
         do i=1,k; sample(i)=values(1+min(size(values)-1,int(rand_uniform()*size(values)))); end do
      else
         p=randperm(size(values)); sample=values(p(:k))
      end if
   end function randsample

   function akimaInterp(x,y,xi) result(yi)
      real(dp),intent(in)::x(:),y(:),xi(:)
      real(dp),allocatable::yi(:)
      type(pchip_result)::pp
      pp=akima(x,y); yi=ppval(pp,xi)
   end function akimaInterp

   function ppfit(x,y) result(pp)
      real(dp),intent(in)::x(:),y(:)
      type(pchip_result)::pp
      pp=cubicspline(x,y)
   end function ppfit

   function lsqlin(a,b,c,d,tolerance,status) result(x)
      real(dp),intent(in)::a(:,:),b(:)
      real(dp),intent(in),optional::c(:,:),d(:),tolerance
      integer,intent(out),optional::status
      real(dp),allocatable::x(:),xc(:),z(:,:),xn(:)
      real(dp)::tol
      integer::st
      tol=1e-12_dp; if(present(tolerance))tol=tolerance; st=pracma_ok
      if(.not.present(c).and..not.present(d))then; x=matmul(pinv(a),b)
      else if(present(c).and.present(d))then
         xc=matmul(pinv(c),d)
         if(maxval(abs(matmul(c,xc)-d))>tol)then; allocate(x(0)); st=pracma_invalid_argument
         else
            z=nullspace(c)
            if(size(z,2)==0)then; x=xc
            else; xn=matmul(pinv(matmul(a,z)),b-matmul(a,xc)); x=xc+matmul(z,xn); end if
         end if
      else; allocate(x(0)); st=pracma_invalid_argument; end if
      if(present(status))status=st
   end function lsqlin

   function lsqlincon(a,b,c,d,tolerance,status) result(x)
      real(dp),intent(in)::a(:,:),b(:),c(:,:),d(:)
      real(dp),intent(in),optional::tolerance
      integer,intent(out),optional::status
      real(dp),allocatable::x(:)
      x=lsqlin(a,b,c,d,tolerance,status)
   end function lsqlincon

   pure function kron(a,b) result(c)
      real(dp),intent(in)::a(:,:),b(:,:)
      real(dp)::c(size(a,1)*size(b,1),size(a,2)*size(b,2))
      integer::i,j,mb,nb
      mb=size(b,1); nb=size(b,2)
      do j=1,size(a,2); do i=1,size(a,1); c((i-1)*mb+1:i*mb,(j-1)*nb+1:j*nb)=a(i,j)*b; end do; end do
   end function kron

   function trisolve(diag,upper,lower,rhs,status) result(x)
      real(dp),intent(in)::diag(:),upper(:),lower(:),rhs(:)
      integer,intent(out),optional::status
      real(dp),allocatable::x(:),a(:),b(:),c(:),d(:)
      real(dp)::m
      integer::i,n,st
      n=size(diag); st=pracma_ok; allocate(x(n),a(n),b(n-1),c(n-1),d(n)); a=diag; b=upper; c=lower; d=rhs
      do i=2,n
         if(abs(a(i-1))<=tiny(1.0_dp))then; st=pracma_singular; x=huge(1.0_dp); if(present(status))status=st; return; end if
         m=c(i-1)/a(i-1); a(i)=a(i)-m*b(i-1); d(i)=d(i)-m*d(i-1)
      end do
      x(n)=d(n)/a(n); do i=n-1,1,-1; x(i)=(d(i)-b(i)*x(i+1))/a(i); end do
      if(present(status))status=st
   end function trisolve

   logical function segm_intersect(s1,s2) result(hit)
      real(dp),intent(in)::s1(2,2),s2(2,2)
      real(dp)::p(2)
      p=segment_intersection(s1(:,1),s1(:,2),s2(:,1),s2(:,2),hit)
   end function segm_intersect

   real(dp) function segm_distance(p1,p2,p3,p4) result(d)
      real(dp),intent(in)::p1(2),p2(2),p3(2)
      real(dp),intent(in),optional::p4(2)
      if(present(p4))then; d=segment_distance(p1,p2,p3,p4); else; d=point_seg(p3,p1,p2); end if
   end function segm_distance

   function poly_crossings(x,y,level) result(points)
      real(dp),intent(in)::x(:),y(:),level
      real(dp),allocatable::points(:,:)
      points=polygon_crossings(x,y,level)
   end function poly_crossings

   pure function tril(a,k) result(l)
      real(dp),intent(in)::a(:,:); integer,intent(in),optional::k
      real(dp)::l(size(a,1),size(a,2)); integer::i,j,kk
      kk=0; if(present(k))kk=k; l=0
      do j=1,size(a,2); do i=1,size(a,1); if(j-i<=kk)l(i,j)=a(i,j); end do; end do
   end function tril

   pure function triu(a,k) result(u)
      real(dp),intent(in)::a(:,:); integer,intent(in),optional::k
      real(dp)::u(size(a,1),size(a,2)); integer::i,j,kk
      kk=0; if(present(k))kk=k; u=0
      do j=1,size(a,2); do i=1,size(a,1); if(j-i>=kk)u(i,j)=a(i,j); end do; end do
   end function triu

   function squareform_vector(v) result(a)
      real(dp),intent(in)::v(:)
      real(dp),allocatable::a(:,:)
      integer::n,i,j,k
      n=int((1+sqrt(1.0_dp+8.0_dp*size(v)))/2); allocate(a(n,n)); a=0; k=0
      do i=1,n-1; do j=i+1,n; k=k+1; a(i,j)=v(k); a(j,i)=v(k); end do; end do
   end function squareform_vector

   function squareform_matrix(a) result(v)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable::v(:)
      integer::n,i,j,k
      n=size(a,1); allocate(v(n*(n-1)/2)); k=0
      do i=1,n-1; do j=i+1,n; k=k+1; v(k)=a(i,j); end do; end do
   end function squareform_matrix

   function flipdim(a,dim) result(b)
      real(dp),intent(in)::a(:,:); integer,intent(in)::dim
      real(dp),allocatable::b(:,:)
      integer::i
      allocate(b(size(a,1),size(a,2)))
      if(dim==1)then; do i=1,size(a,1); b(i,:)=a(size(a,1)-i+1,:); end do
      else; do i=1,size(a,2); b(:,i)=a(:,size(a,2)-i+1); end do; end if
   end function flipdim

   pure elemental real(dp) function softline(x,alpha) result(y)
      real(dp),intent(in)::x
      real(dp),intent(in),optional::alpha
      real(dp)::a
      a=1.0_dp; if(present(alpha))a=alpha; y=log(1.0_dp+exp(-abs(a*x)))/a+max(x,0.0_dp)
   end function softline

   function piecewise(x,breaks,values) result(y)
      real(dp),intent(in)::x(:),breaks(:),values(:)
      real(dp),allocatable::y(:)
      integer::i,j
      allocate(y(size(x)))
      do i=1,size(x); j=1; do while(j<size(breaks).and.x(i)>=breaks(j+1)); j=j+1; end do; y(i)=values(min(j,size(values))); end do
   end function piecewise

   function divisors(n,n0) result(d)
      integer,intent(in)::n
      integer,intent(in),optional::n0
      integer,allocatable::d(:),tmp(:)
      integer::lo,i,k
      lo=1; if(present(n0))lo=n0; allocate(tmp(max(0,n/2))); k=0
      do i=max(1,lo),n/2; if(mod(n,i)==0)then; k=k+1; tmp(k)=i; end if; end do
      allocate(d(k)); d=tmp(:k)
   end function divisors

   real(dp) function rssimple(x) result(h)
      real(dp),intent(in)::x(:)
      real(dp),allocatable::s(:)
      real(dp)::sd
      integer::i,n
      n=size(x); allocate(s(n)); do i=1,n; s(i)=sum(x(:i)-sum(x)/n); end do
      sd=sqrt(sum((x-sum(x)/n)**2)/max(1,n-1)); h=log((maxval(s)-minval(s))/sd)/log(real(n,dp))
   end function rssimple

   real(dp) function hurstexp(x) result(h)
      real(dp),intent(in)::x(:); h=hurst(x)
   end function hurstexp

   real(dp) function approx_entropy(ts,edim,r) result(apen)
      real(dp),intent(in)::ts(:)
      integer,intent(in),optional::edim
      real(dp),intent(in),optional::r
      integer::m,n,i,j,k,mm
      real(dp)::tol,phi(2),cnt
      n=size(ts); mm=2; if(present(edim))mm=edim; tol=0.2_dp*std_local(ts); if(present(r))tol=r
      do k=1,2
         m=mm+k-1; phi(k)=0
         do i=1,n-m+1
            cnt=0
            do j=1,n-m+1; if(maxval(abs(ts(i:i+m-1)-ts(j:j+m-1)))<=tol)cnt=cnt+1; end do
            phi(k)=phi(k)+log(cnt/real(n-m+1,dp))
         end do
         phi(k)=phi(k)/real(n-m+1,dp)
      end do
      apen=phi(1)-phi(2)
   end function approx_entropy

   real(dp) function sample_entropy(ts,edim,r,tau) result(sampen)
      real(dp),intent(in)::ts(:)
      integer,intent(in),optional::edim,tau
      real(dp),intent(in),optional::r
      real(dp),allocatable::x(:)
      integer::m,t,n,i,j,cm,cm1
      real(dp)::tol
      m=2; if(present(edim))m=edim; t=1; if(present(tau))t=tau; x=ts(1:size(ts):t); n=size(x)
      tol=0.2_dp*std_local(x); if(present(r))tol=r; cm=0; cm1=0
      do i=1,n-m-1; do j=i+1,n-m
         if(maxval(abs(x(i:i+m-1)-x(j:j+m-1)))<=tol)then; cm=cm+1; if(abs(x(i+m)-x(j+m))<=tol)cm1=cm1+1; end if
      end do; end do
      if(cm1==0.or.cm==0)then; sampen=huge(1.0_dp); else; sampen=-log(real(cm1,dp)/real(cm,dp)); end if
   end function sample_entropy

   function bubbleSort(x) result(y)
      real(dp),intent(in)::x(:); real(dp),allocatable::y(:); integer::i,j; real(dp)::t
      y=x; do i=1,size(y)-1; do j=1,size(y)-i; if(y(j)>y(j+1))then; t=y(j); y(j)=y(j+1); y(j+1)=t; end if; end do; end do
   end function bubbleSort
   function insertionSort(x) result(y)
      real(dp),intent(in)::x(:); real(dp),allocatable::y(:); integer::i,j; real(dp)::t
      y=x; do i=2,size(y); t=y(i); j=i-1; do while(j>=1); if(y(j)<=t)exit; y(j+1)=y(j); j=j-1; end do; y(j+1)=t; end do
   end function insertionSort
   function selectionSort(x) result(y)
      real(dp),intent(in)::x(:); real(dp),allocatable::y(:); integer::i,j,k; real(dp)::t
      y=x; do i=1,size(y)-1; k=i; do j=i+1,size(y); if(y(j)<y(k))k=j; end do; t=y(i); y(i)=y(k); y(k)=t; end do
   end function selectionSort
   function shellSort(x) result(y)
      real(dp),intent(in)::x(:); real(dp),allocatable::y(:); integer::gap,i,j; real(dp)::t
      y=x; gap=size(y)/2
      do while(gap>0)
         do i=gap+1,size(y)
            t=y(i); j=i
            do while(j>gap.and.y(j-gap)>t)
               y(j)=y(j-gap); j=j-gap
            end do
            y(j)=t
         end do
         gap=gap/2
      end do
   end function shellSort
   function heapSort(x) result(y)
      real(dp),intent(in)::x(:); real(dp),allocatable::y(:); y=insertionSort(x)
   end function heapSort
   function mergeSort(x) result(y)
      real(dp),intent(in)::x(:); real(dp),allocatable::y(:); y=insertionSort(x)
   end function mergeSort
   function quickSort(x) result(y)
      real(dp),intent(in)::x(:); real(dp),allocatable::y(:); y=insertionSort(x)
   end function quickSort

   pure real(dp) function normF(a) result(v)
      real(dp),intent(in)::a(:,:); v=sqrt(sum(a*a))
   end function normF

   function abm3pc(f,tspan,y0,nsteps) result(res)
      procedure(vector_field)::f; real(dp),intent(in)::tspan(2),y0(:); integer,intent(in)::nsteps; type(ode_result)::res
      res=abm3(f,tspan,y0,nsteps)
   end function abm3pc

   function shubert(f,a,b,lipschitz,tolerance,max_iter) result(res)
      procedure(scalar_function)::f
      real(dp),intent(in)::a,b,lipschitz
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter
      type(optimization_result)::res
      real(dp),allocatable::xs(:),zs(:),xt(:),zt(:)
      real(dp)::tol,x0,y0,xbest,ybest,fupper,xn,zn,yn,zl,zr,xl,xr
      integer::n,nmax,iter,i,j
      tol=1.0e-4_dp; if(present(tolerance))tol=tolerance
      nmax=1000; if(present(max_iter))nmax=max_iter
      allocate(xs(nmax+2),zs(nmax+2),xt(nmax+2),zt(nmax+2),res%x(1),res%history(nmax))
      x0=0.5_dp*(a+b); y0=f(x0); xbest=x0; ybest=y0
      fupper=y0+lipschitz*(b-a)/2.0_dp
      xs(1)=b; zs(1)=y0+fupper; xs(2)=a; zs(2)=y0+fupper; n=2
      do iter=1,nmax
         if(fupper-ybest<=tol)exit
         xn=xs(n); zn=zs(n); yn=f(xn)
         if(yn>ybest)then; xbest=xn; ybest=yn; end if
         zl=0.5_dp*(zn+yn); zr=zl
         xl=xn-(zn-yn)/(2.0_dp*lipschitz); xr=xn+(zn-yn)/(2.0_dp*lipschitz)
         n=n-1
         if(xl>=a.and.xl<=b)then; n=n+1; xs(n)=xl; zs(n)=zl; end if
         if(xr>=a.and.xr<=b)then; n=n+1; xs(n)=xr; zs(n)=zr; end if
         ! Stable insertion sort by increasing upper bound.
         do i=2,n
            xn=xs(i); zn=zs(i); j=i-1
            do while(j>=1)
               if(zs(j)<=zn)exit
               xs(j+1)=xs(j); zs(j+1)=zs(j); j=j-1
            end do
            xs(j+1)=xn; zs(j+1)=zn
         end do
         fupper=zs(n); res%history(iter)=ybest
      end do
      res%x(1)=xbest; res%value=ybest; res%iterations=min(iter,nmax)
      res%evaluations=res%iterations+1
      res%converged=fupper-ybest<=tol; res%status=merge(pracma_ok,pracma_invalid_argument,res%converged)
      if(res%iterations<size(res%history))res%history=res%history(:res%iterations)
   end function shubert

   function quickSortx(x) result(y)
      real(dp),intent(in)::x(:)
      real(dp),allocatable::y(:)
      y=quickSort(x)
   end function quickSortx

   function mergeOrdered(a,b) result(y)
      real(dp),intent(in)::a(:),b(:)
      real(dp),allocatable::y(:)
      integer::i,j,k
      allocate(y(size(a)+size(b))); i=1; j=1; k=0
      do while(i<=size(a).and.j<=size(b))
         k=k+1
         if(a(i)<=b(j))then; y(k)=a(i); i=i+1; else; y(k)=b(j); j=j+1; end if
      end do
      if(i<=size(a))y(k+1:)=a(i:)
      if(j<=size(b))y(k+1:)=b(j:)
   end function mergeOrdered

   pure logical function is_sorted(x) result(ok)
      real(dp),intent(in)::x(:)
      ok=size(x)<2.or.all(x(2:)>=x(:size(x)-1))
   end function is_sorted

   function sortrows(a,column) result(b)
      real(dp),intent(in)::a(:,:)
      integer,intent(in),optional::column
      real(dp),allocatable::b(:,:),row(:)
      integer::c,i,j
      c=1; if(present(column))c=column; b=a; allocate(row(size(a,2)))
      do i=2,size(b,1)
         row=b(i,:); j=i-1
         do while(j>=1)
            if(b(j,c)<=row(c))exit
            b(j+1,:)=b(j,:); j=j-1
         end do
         b(j+1,:)=row
      end do
   end function sortrows

   subroutine stereographic(points,projected,pole)
      real(dp),intent(in)::points(:,:)
      real(dp),allocatable,intent(out)::projected(:,:)
      integer,intent(in),optional::pole
      real(dp),allocatable::u(:),v(:)
      call stereographic_project(points(1,:),points(2,:),points(3,:),u,v,pole)
      allocate(projected(2,size(points,2))); projected(1,:)=u; projected(2,:)=v
   end subroutine stereographic

   subroutine stereographic_inv(projected,points,pole)
      real(dp),intent(in)::projected(:,:)
      real(dp),allocatable,intent(out)::points(:,:)
      integer,intent(in),optional::pole
      real(dp),allocatable::x(:),y(:),z(:)
      call stereographic_inverse(projected(1,:),projected(2,:),x,y,z,pole)
      allocate(points(3,size(projected,2))); points(1,:)=x; points(2,:)=y; points(3,:)=z
   end subroutine stereographic_inv

   function lebesgue(x,xi) result(values)
      real(dp),intent(in)::x(:),xi(:)
      real(dp),allocatable::values(:)
      values=lebesgue_function(x,xi)
   end function lebesgue

   function half(indices) result(out)
      integer,intent(in)::indices(:)
      integer,allocatable::out(:)
      integer::i,k,n
      n=size(indices); allocate(out(2*n-1)); k=0
      do i=1,n-1
         k=k+1; out(k)=indices(i)
         k=k+1; out(k)=indices(i)+(indices(i+1)-indices(i)+1)/2
      end do
      out(2*n-1)=indices(n)
   end function half

   function bvp(ffun,gfun,hfun,xspan,ybound,n) result(res)
      procedure(scalar_function)::ffun,gfun,hfun
      real(dp),intent(in)::xspan(2),ybound(2)
      integer,intent(in),optional::n
      type(ode_result)::res
      real(dp),allocatable::diag(:),upper(:),lower(:),rhs(:),interior(:)
      real(dp)::dt,x
      integer::nn,i,st
      nn=50; if(present(n))nn=n
      if(nn<2.or.xspan(2)<=xspan(1))then; res%status=pracma_invalid_argument; return; end if
      dt=(xspan(2)-xspan(1))/real(nn+1,dp)
      allocate(res%t(nn+2),res%y(1,nn+2),diag(nn),upper(nn-1),lower(nn-1),rhs(nn))
      res%t=[(xspan(1)+dt*real(i,dp),i=0,nn+1)]
      do i=1,nn
         x=res%t(i+1); diag(i)=-2.0_dp-dt*dt*gfun(x); rhs(i)=dt*dt*hfun(x)
         if(i<nn)upper(i)=1.0_dp-0.5_dp*dt*ffun(x)
         if(i>1)lower(i-1)=1.0_dp+0.5_dp*dt*ffun(x)
      end do
      rhs(1)=rhs(1)-ybound(1)*(1.0_dp+0.5_dp*dt*ffun(res%t(2)))
      rhs(nn)=rhs(nn)-ybound(2)*(1.0_dp-0.5_dp*dt*ffun(res%t(nn+1)))
      interior=trisolve(diag,upper,lower,rhs,st)
      res%y(1,1)=ybound(1); res%y(1,2:nn+1)=interior; res%y(1,nn+2)=ybound(2)
      res%accepted_steps=nn+1; res%converged=st==pracma_ok; res%status=st
   end function bvp

   pure real(dp) function point_seg(p,a,b) result(d)
      real(dp),intent(in)::p(2),a(2),b(2)
      real(dp)::ab(2),t
      ab=b-a
      if(sum(ab*ab)==0.0_dp)then
         d=sqrt(sum((p-a)**2))
      else
         t=max(0.0_dp,min(1.0_dp,dot_product(p-a,ab)/sum(ab*ab)))
         d=sqrt(sum((p-a-t*ab)**2))
      end if
   end function point_seg

   real(dp) function std_local(x) result(s)
      real(dp),intent(in)::x(:); real(dp)::m
      m=sum(x)/size(x); s=sqrt(sum((x-m)**2)/max(1,size(x)-1))
   end function std_local

end module pracma_compat
