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
c       (This is all done automatically in the R package "expm".)
c
c -- MM: This is *legacy* code - but we provide the "padeO" ... methods
c -- ===>  Fix the fortran code enough that it does not give "--as-cran" warnings:
c
      subroutine matrexpO(a, ndim, ntaylor, npade, accuracy)

      integer ndim, ntaylor, npade
      double precision a(ndim,ndim), accuracy


c	"ndim" is the order of the given matrix A
c
      double precision sum(ndim,ndim)
      double precision solution(ndim,ndim)
      double precision error(ndim,ndim)
      double precision dkeep(ndim,ndim)

      double precision dsqrt, dl1norm, dlinfnorm
      integer log2

c
c	use the algorithm to compute T(ntaylor) or P(npade)
c
      npower=log2(dsqrt(dl1norm(ndim,a)*dlinfnorm(ndim,a)))+4
      if(ntaylor.gt.0)then
        	call taylorO(ndim,ntaylor,npower,a,sum)
      else
        	call padeO(ndim,npade,npower,a,sum)
      endif
      call powermatrix(ndim,sum,npower,dkeep)
      call id(ndim,dkeep,sum)
c
c	computing the "solution" T(ntaylor+10) or P(npade+10)
c
      if(ntaylor.gt.0)then
        call taylorO(ndim,ntaylor+10,npower,a,solution)
      else
        call padeO(ndim,npade+10,npower,a,solution)
      endif
      call powermatrix(ndim,solution,npower,dkeep)
      call id(ndim,dkeep,solution)
c
c       copy the result back into a
c
      do i=1,ndim
         do j=1,ndim
            a(i,j) = sum(i,j)
         end do
      end do
c
c	compute the Cauchy error T(ntaylor+10)-T(ntaylor)
c                             or P(npade+10)-P(npade)
c
      call subtract(ndim,sum,solution,error)
      accuracy = dsqrt(dl1norm(ndim,error)*dlinfnorm(ndim,error))
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine taylorO(m,ntaylor,npower,a,sum)

c	Taylor series for exp(a/2**npower)

      implicit double precision (a-h,o-z)
      double precision a(m,m),sum(m,m),dkeep(m,m)
      nscale=2**npower
c	print*,'A is scaled by 2**',npower,'  =',nscale
      call initialize(m,sum,0.d0)
      call addtodiag(m,sum,1.d0)
      do n=ntaylor,1,-1
        call multiplymatrixO(m,sum,a,dkeep)
        call multiplyscalarO(m,dkeep,1.d0/dble(n*nscale),sum)
        call addtodiag(m,sum,1.d0)
      end do
      return
      end


      subroutine padeO(m,npade,npower,a,approx)

c	Pade approximation for exp(a/2**npower)

      integer m, npade, npower
      double precision a(m,m), approx(m,m)
c Var
      double precision aminus(m,m),dkeep(m,m)
      double precision padenom(m,m),padedenom(m,m)
      integer nscale,n,i

      nscale=2**npower
c	print*,'A is scaled by 2**',npower,'  =',nscale
      call initialize(m,padenom,0.d0)
      call addtodiag(m,padenom,1.d0)
      do n=npade,1,-1
        call multiplymatrixO(m,padenom,a,dkeep)
        call multiplyscalarO(m,dkeep,
     $      dble(npade-n+1)/dble(n*(2*npade-n+1)*nscale),padenom)
        call addtodiag(m,padenom,1.d0)
      end do
      call minus(m,a,aminus)
      call initialize(m,padedenom,0.d0)
      call addtodiag(m,padedenom,1.d0)
      do n=npade,1,-1
        call multiplymatrixO(m,padedenom,aminus,dkeep)
        call multiplyscalarO(m,dkeep,
     $      dble(npade-n+1)/dble(n*(2*npade-n+1)*nscale),padedenom)
        call addtodiag(m,padedenom,1.d0)
      end do
      do i=1,m
         call solveO(m,padedenom,padenom(1,i),approx(1,i))
      end do
      return
      end

c MM: FIXME: Use BLAS for the following !!!
c     -----  (also R internal things --> use C, not Fortran !

c	initializing a matrix to a scalar s

      subroutine initialize(m,x,s)

      integer m
      double precision x(m,m), s

      integer i,j
      do i=1,m
         do j=1,m
            x(i,j)=s
         end do
      end do
      return
      end

      subroutine multiplyscalarO(m,x,s,y)

c	multiplying a matrix x by a scalar s

      implicit double precision (a-h,o-z)
      double precision x(m,m),y(m,m)
      do i=1,m
        do j=1,m
            y(i,j)=x(i,j)*s
        end do
      end do
      return
      end

      subroutine multiplymatrixO(m,x,y,z)

c	multiplying 2 matrices

      implicit double precision (a-h,o-z)
      double precision x(m,m),y(m,m),z(m,m)
      do i=1,m
        do j=1,m
          z(i,j)=0.d0
          do k=1,m
              z(i,j)=z(i,j)+x(i,k)*y(k,j)
          end do
        end do
      end do
      return
      end


      subroutine id(m,x,y)

c	assign a matrix x to y

      implicit double precision (a-h,o-z)
      double precision x(m,m),y(m,m)
      do i=1,m
        do j=1,m
            y(i,j)=x(i,j)
        end do
      end do
      return
      end


      subroutine powermatrix(m,x,ipower,y)

c	computing the ith power of a matrix x

      implicit double precision (a-h,o-z)
      double precision x(m,m),y(m,m),dkeep(m,m)
      call id(m,x,y)
      call id(m,x,dkeep)
      do i=1,ipower
        call multiplymatrixO(m,dkeep,dkeep,y)
        call id(m,y,dkeep)
      end do
      return
      end


      subroutine iden(m,x,y)

c	assign a vector x to y

      implicit double precision (a-h,o-z)
      double precision x(m),y(m)
      do i=1,m
         y(i)=x(i)
      end do
      return
      end


      double precision function dip(m,u,v)

c	inner product of 2 vectors

      integer m
      double precision u(m),v(m)

      integer i
      dip=0.d0
      do i=1,m
         dip = dip+u(i)*v(i)
      end do
      return
      end


      double precision function dl2norm(m,u)

c	l2 norm of a vector

      implicit double precision (a-h,o-z)
      double precision u(m)
      dl2norm=dsqrt(dip(m,u,u))
      return
      end


      subroutine solveO(m,A,f,x)

c	CGS iteration

      integer m
      double precision A(m,m), f(m), x(m)
c
      double precision save(m),rcgs(m),r(m)
      double precision p(m),u(m)
      double precision rbar(m),v(m),q(m)

      external dl2norm, dip
      double precision dl2norm, dip

      double precision alpha, beta, thresh, eps, omega0,omega1,
     +          omegainit, rho0,rho1, scalar,sigma, tau,vv
      integer l


      thresh=1.d-100
      eps=1.d-30
      call zero(m,x)
      call iden(m,f,r)
      call iden(m,r,rcgs)
      call iden(m,r,p)
      call iden(m,r,u)
      omega0= dl2norm(m,rcgs)
      omegainit=omega0
c	print*,'res0=',dabs(omegainit)
      tau=omega0
      vv=0.d0
      eta=0.d0
      call iden(m,r,rbar)
      rho0= dip(m,rbar,r)
      if(dabs(rho0).le.thresh)then
c	  print*,'rho0=',rho0,'   MG iteration number=1'
        return
      endif
      do l=1,m
        call multiplyvector(m,A,p,v)
        sigma=dip(m,rbar,v)
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
        omega1=dl2norm(m,rcgs)
        call comb(m,x,alpha,v,x)
c	  print*,'residual=',dabs(omega1),'  iteration number=',2*l+1
        if(dabs(omega1/omegainit).le.eps)then
c	    print*,'iteration number=',2*l+1
          return
        endif
        omega0=omega1

        rho1=dip(m,rbar,rcgs)
        if(dabs(rho1).le.thresh)then
c	    print*,'rho1=',rho1,'   iteration number=',2*l+1
          return
        endif
        beta=rho1/rho0
        rho0=rho1
        call comb(m,rcgs,beta,q,u)
        call comb(m,q,beta,p,save)
        call comb(m,u,beta,save,p)
      end do
c	print*,'iteration number=',2*l+1
      return
      end
