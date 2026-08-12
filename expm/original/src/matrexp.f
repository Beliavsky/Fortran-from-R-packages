cccc-*- mode: fortran; kept-old-versions: 12;  kept-new-versions: 20; -*-

c This program computes exp(A) for a given matrix A.
c
c 2 algorithms are employed:
c    1. The Taylor expansion of order "ntaylor," denoted by T(ntaylor).
c    2. The Pade diagonal approximation of order
c       "npade" denoted by P(npade), is used instead, IFF  ntaylor = 0

c       The algorithm is applied twice to calculate
c       T(ntaylor) and T(ntaylor+10) [or, when ntaylor=0,
c       to calculate P(npade) and P(npade+10)].

c	An upper bound for the L2 norm of the Cauchy error
c       T(ntaylor+10)-T(ntaylor) [or, when ntaylor=0,
c       P(npade+10)-P(npade)] is computed.

c       The result exp(A) is returned via the first argument.
c
c	This version works with R (i.e., is written as a subroutine)
c
c	To use it, first do
c	% R SHLIB matrexp.f
c       to make the shared library matrexp.so
c
c	and then in R,
c	> dyn.load("matrexp.so")
c	> .C("matrexp_",as.double(runif(9)),
c            as.integer(3),as.integer(0),as.integer(8))
c
c       (This is all done automatically in the R package 'expm'.)
c
      subroutine matrexp(a, n, ntaylor, npade, accuracy)

      integer n, ntaylor, npade
      double precision a(n,n), accuracy


c	"n" is the order of the given matrix A
c
      double precision sum(n,n), sol10(n,n)

      double precision dsqrt, dl1norm, dlinfnorm
      integer log2

      integer npower, i,j

c FIXME: consider computing
c       sqrt(dl1norm(n,a) * dlinfnorm(n,a))
c in one function -- no need for the two separate l*norm() functions
c
      npower= log2(dsqrt(dl1norm(n,a)*dlinfnorm(n,a))) + 4
c
c Use the algorithm to compute T(ntaylor) or P(npade)
c
      if(ntaylor.gt.0)then
         call taylor(n,ntaylor,npower,a,sum)
      else
         call pade(n,npade,npower,a,sum)
      endif
c
c	computing the "solution" T(ntaylor+10) or P(npade+10)
c
      if(ntaylor.gt.0) then
        call taylor(n,ntaylor+10,npower,a,sol10)
      else
        call pade(n,npade+10,npower,a,sol10)
      endif

      call powMat(n,sum, npower)
c       copy the result back into a
c
      do i=1,n
         do j=1,n
            a(i,j) = sum(i,j)
         end do
      end do

      call powMat(n,sol10, npower)

c
c	compute the Cauchy error T(ntaylor+10)- T(ntaylor)
c                             or P(npade+10)  - P(npade)
c
cc- BLAS:  daxpy(m,  alpha, x, 1, y, 1)  :  y := y + alpha*x
c
       call subtract(n,sum,sol10, sum)
       accuracy = dsqrt(dl1norm(n,sum) * dlinfnorm(n,sum))

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine taylor(m,ntaylor,npower,a,sum)

c	Taylor series for exp(a/2**npower)

      integer m, ntaylor, npower
      double precision a(m,m),sum(m,m)

      double precision T(m,m)
      integer nscale,n

      nscale=2**npower
c	print*,'A is scaled by 2**',npower,'  =',nscale

      call identity(m, sum)

      do n=ntaylor,1,-1

C FIXME: use multScalAdd() instead of these three
C -----  add tests for this  taylor code  !! --> ../tests/ex.

        call multiplymatrix(m,sum,a,T)
        call multiplyscalar(m,T,1.d0/dble(n*nscale),sum)
        call addtodiag(m,sum,1.d0)

C   ---------  1) multiplymatrix(., A,. B,. C) ; 2) multiplyscalar(., C, s, D)
C              3) addtodiag (D, ., 1.0)
C --->  D := s (A B)  + Id() ---> do this directly via DGEMM :
C       D := Id() ;  D := s * A B + 1 * D
C__ almost:  call multScalAdd(m, 1.d0/dble(n*nscale), sum,a, T)
c  NOTE: result must be 'sum'
      enddo
      return
      end


      subroutine pade(m,npade,npower,a,approx)

c	Pade approximation for exp(a/2**npower)

      integer m, npade, npower
      double precision a(m,m), approx(m,m)
c Var
      double precision aminus(m,m),T(m,m)
      double precision padenom(m,m),padedenom(m,m)
      integer nscale,n,i

      nscale=2**npower
c	print*,'A is scaled by 2**',npower,'  =',nscale
      call identity(m,padenom)
      call identity(m,padedenom)

      do  n=npade,1,-1
         call multiplymatrix(m,padenom,a,T)
         call multiplyscalar(m,T,
     $        dble(npade-n+1)/dble(n*(2*npade-n+1)*nscale),padenom)
         call addtodiag(m,padenom,1.d0)
      enddo
      call minus(m,a,aminus)

      do n=npade,1,-1
        call multiplymatrix(m,padedenom,aminus,T)
        call multiplyscalar(m,T,
     $      dble(npade-n+1)/dble(n*(2*npade-n+1)*nscale),padedenom)
        call addtodiag(m,padedenom,1.d0)
      enddo
      do i=1,m
         call solve(m,padedenom,padenom(1,i),approx(1,i))
      end do
      return
      end

c MM: FIXME: Use BLAS for the following !!!
c     -----  (also R internal things --> use C, not Fortran !


c  initializing an identity matrix

      subroutine identity(m,x)

      integer m
      double precision x(m,m)

      integer i,j
      do i=1,m
         do j=1,m
            x(i,j)= 0d0
         enddo
         x(i,i)= 1d0
      enddo
      return
      end


C-- NOTA BENE: We also have
C   ---------  1) multiplymatrix(., A,. B,. C) ; 2) multiplyscalar(., C, s, D)
C              3) addtodiag (D, ., 1.0)
C --->  D := s (A B)  + Id() ---> do this directly via DGEMM :
C       D := Id() ;  D := s * A B + 1 * D

      subroutine multScalAdd(m,s, x,y,z)

      integer m
      double precision x(m,m),y(m,m),z(m,m), s

c  GEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA,  C, LDC )
c      C := alpha * AB  + beta * C

      call identity(m, z)
      call dgemm('N','N', m,m,m, s,  x, m,  y, m, 1.d0, z, m)

      return
      end

c  DSCAL ( N,  ALPHA, X,  1) :  x <- alpha * x  (x = x[1:n])

c FIXME: replace by using DSCAL :

      subroutine multiplyscalar(m,x,s,y)

c	multiplying a matrix x by a scalar s		Y :=  s X

      integer m
      double precision x(m,m),y(m,m), s

      integer i,j
      do i=1,m
        do j=1,m
           y(i,j)=x(i,j)*s
        enddo
      enddo
      return
      end

c  DGEMM (TRANSA, TRANSB,      M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
c      C := alpha A B + beta C

      subroutine multiplymatrix(m,x,y,z)

c	multiplying two   m x m  matrices		Z := X %*% Y

      integer m
      double precision x(m,m),y(m,m),z(m,m)

c  GEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA,  C, LDC )
c      C := a * AB + b * C

      call dgemm('N','N', m,m,m, 1.d0,  x, m,  y, m, 0.d0, z, m)

      return
      end


      subroutine powMat(m,x,ipower)

c Compute  x ^ (2^i) -- by simple squaring


      integer m, ipower
      double precision x(m,m)

      double precision xx(m,m)
      integer i, m2
      m2 = m*m

      call dcopy(m2, x,1, xx,1)
      do i=1,ipower
        call multiplymatrix(m,xx,xx,x)
        if(i .lt. ipower) call dcopy(m2, x,1, xx,1)
      enddo
      return
      end


cc-- BLAS:  daxpy(m,  alpha, x, 1, y, 1)  :  y := y + alpha*x

C this is now in mexp-common.f --  FIXME:  use daxpy() or ...
c
c      subroutine comb(m, x,a,y, z)
c                 ----
c	linear combination of 2 vectors		z[] := x[] + a* y[]



C  BLAS   DDOT ( N, X, INCX, Y, INCY )  --- inner product  X' Y


      subroutine solve(m,A,f,x)

c	CGS iteration

      integer m
      double precision A(m,m), f(m), x(m)
c
      double precision save(m),rcgs(m),r(m)
      double precision p(m),u(m)
      double precision rbar(m),v(m),q(m)

      external dnrm2, ddot
      double precision dnrm2, ddot

      double precision alpha, beta, thresh, eps, eta, omega0,omega1,
     +     omegainit, rho0,rho1, scalar,sigma, tau,vv
      integer l

      thresh=1.d-100
      eps=1.d-30
      call zero(m,x)
      call dcopy(m, f,1, r,   1)
      call dcopy(m, r,1, rcgs,1)
      call dcopy(m, r,1, p,   1)
      call dcopy(m, r,1, u,   1)
      omega0= dnrm2(m,rcgs, 1)
      omegainit=omega0
c	print*,'res0=',dabs(omegainit)
      tau=omega0
      vv=0.d0
      eta=0.d0
      call dcopy(m, r,1, rbar,1)
      rho0= ddot(m, rbar,1, r,1)
      if(dabs(rho0).le.thresh)then
c	  print*,'rho0=',rho0,'   MG iteration number=1'
        return
      endif
      do 10 l=1,m
        call multiplyvector(m,A,p,v)
        sigma=ddot(m, rbar,1, v,1)
        if(dabs(sigma).le.thresh)then
c	    print*,'sigma=',sigma,'  iteration number=',2*l+1
          return
        endif
        alpha=rho0/sigma
        if(dabs(alpha).le.thresh)then
c	    print*,'alpha=',alpha,'  iteration number=',2*l+1
          return
        endif
        scalar=-alpha
        call comb(m,u,scalar,v,q)
        call add(m,u,q,v)
        call multiplyvector(m,A,v,save)
        call comb(m,rcgs,scalar,save,rcgs)
        omega1=dnrm2(m,rcgs, 1)
        call comb(m,x,alpha,v,x)
c	  print*,'residual=',dabs(omega1),'  iteration number=',2*l+1
        if(dabs(omega1/omegainit).le.eps)then
c	    print*,'iteration number=',2*l+1
          return
        endif
        omega0=omega1

        rho1=ddot(m, rbar,1, rcgs,1)
        if(dabs(rho1).le.thresh)then
c	    print*,'rho1=',rho1,'   iteration number=',2*l+1
          return
        endif
        beta=rho1/rho0
        rho0=rho1
c       u[]    := rcgs[] + beta q[]
        call comb(m,rcgs,beta,q,u)
c       save[] := q[] + beta p[]
        call comb(m,q,beta,p,save)
c       p[] := u[] + beta save[] = (rcgs[] + beta q[]) + beta(q[] + beta p[])
        call comb(m,u,beta,save,p)
c	print*,'iteration number=',2*l+1
 10   continue
      return
      end
