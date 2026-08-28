! SPDX-License-Identifier: GPL-2.0-or-later
! Converted nleqslv numerical kernels collected in a module so all internal
! solver calls have explicit interfaces under modern Fortran compilers.
module nleqslv_legacy_kernels
   use nleqslv_blas_lapack
   use nleqslv_trace_state, only : nwsnot, nwckot, nwjerr, nwprot, nwlsot, nwdgot, nwpwot, nwmhot
   implicit none (type, external)
   private
   public :: liqsiz, nwnleq

   abstract interface
      subroutine legacy_fvec(x, f, n, flag)
         integer, intent(in) :: n, flag
         double precision, intent(in) :: x(*)
         double precision, intent(out) :: f(*)
      end subroutine legacy_fvec
      subroutine legacy_fjac(rjac, ldr, x, n)
         integer, intent(in) :: ldr, n
         double precision, intent(out) :: rjac(ldr,*)
         double precision, intent(in) :: x(*)
      end subroutine legacy_fjac
   end interface
contains

! ============================================================================
! From converted legacy source: lautil.f90
! ============================================================================

!-----------------------------------------------------------------------------

subroutine liqrfa(a, lda, n, tau, work, wsiz, info)
integer  lda, n, wsiz, info
double precision  a(lda,*), tau(*), work(*)

!-------------------------------------------------------------
!
!     QR decomposition of A(n,n)
!
!     Arguments
!
!      Inout A        Real(Lda,n)    Matrix to transform.
!      In    lda      Integer        Leading dimension of A
!      In    n        Integer        number of rows/cols A
!      Out   tau      Real(n)        Information for recovering
!      Out   work     Real(*)        workspace
!      In    wsiz     Integer        size of work()

!     Lapack blocked QR (much faster for larger n)
!
!-------------------------------------------------------------

call dgeqrf(n,n,A,lda,tau,work,wsiz,info)

return
end

!=============================================================

subroutine liqsiz(n,wrksiz)
integer n, wrksiz

!-------------------------------------------------------------------------
!     Query the size of the double precision work array required
!     for optimal operation of the Lapack QR routines
!-------------------------------------------------------------------------

double precision A(1), work(1)
integer lwork, info

lwork = -1
call dgeqrf(n,n,A,n,work,work,lwork,info)
if( info .ne. 0 ) then
wrksiz = -1
else
wrksiz = int(work(1))
endif

return
end

!=============================================================

subroutine liqrqt(a, lda, n, tau, qty, work, wsiz, info)
integer lda, n, wsiz, info
double precision a(lda,*), tau(*), qty(*), work(*)

!-------------------------------------------------------------
!      Arguments
!
!      In    A     Real(Lda, n)    QR decomposition
!      In    Lda   Integer         Leading dimension A
!      In    n     Integer         Number of rows/columns in A
!      In    tau   Integer         Householder constants from QR
!      Inout qty   Real(n)         On input, vector y
!                                  On output, trans(Q)*y
!      Out   work  Real(*)         workspace
!      In    wsiz  Integer         size of work()
!
!     Liqrqt calculates trans(Q)*y from the QR decomposition
!
!     Lapack blocked
!-------------------------------------------------------------

call dormqr('L','T',n,1,n,A,lda,tau,qty,n,work,wsiz,info)

return
end

!=============================================================

subroutine liqrqq(q,ldq,tau,n,work,wsiz,info)
integer n, ldq, wsiz, info
double precision  q(ldq,*),tau(*),work(*)

!     Arguments
!
!     Inout  Q     Real(ldq,*)     On input, QR decomposition
!                                    On output, the full Q
!     In     ldq   Integer         leading dimension of Q
!     In     tau   Real(n)         Householder constants of
!                                     the QR decomposition
!     In     n     Integer         number of rows/columns in Q
!     Out    work  Real(n)         workspace of length n
!     In     wsiz  Integer         size of work()
!
!     Generate Q from QR decomposition Liqrfa (dgeqr2)
!
!     Lapack blocked
!-------------------------------------------------------------

call dorgqr(n,n,n,q,ldq,tau,work,wsiz,info)

return
end

!-----------------------------------------------------------------------------

subroutine nuzero(n,x)
integer n
double precision x(*)

!     Parameters:
!
!     In    n        Integer           Number of elements.
!     In    x        Real(*)           Vector of reals.
!
!     Description:
!
!     Nuzero sets all elements of x to 0.
!     Does nothing when n <= 0

double precision Rzero
parameter(Rzero=0.0d0)

integer i

do i=1,n
x(i) = Rzero
enddo

return
end

! ----------------------------------------------------------------------

subroutine nuvgiv(x,y,c,s)
double precision x,y,c,s

!     Parameters
!
!     Inout   x     Real       x input / c*x+s*y on output
!     Inout   y     Real       y input / 0       on output
!     Out     c     Real       c of tranformation (cosine)
!     Out     s     Real       s of tranformation (  sine)
!
!     Description
!
!     Nuvgiv calculates the givens rotator
!
!             |  c   s |
!         G = |        |
!             | -s   c |
!
!     with  c*c+s*s=1
!
!     for which G * | x | = | t |
!                   | y |   | 0 |
!
!     resulting in
!
!            c * x + s * y = t
!           -s * x + c * y = 0   ==>  s/c = y/x or c/s = x/y
!
!     Use Lapack dlartg routine
!     return c and s and the modified x and y
!     This differs from dlartg which does not modify input arguments.
!     See http://www.netlib.org/lapack/explore-html/dd/d24/dlartg_8f_source.html
!     c * x + s * y may differ from t with machine precision

double precision t

double precision Rzero
parameter(Rzero=0.0d0)

call dlartg(x,y,c,s,t)
x = t
y = Rzero
return
end

! ============================================================================
! From converted legacy source: limhpar.f90
! ============================================================================
subroutine limhpar(R, ldr, n, sdiag, qtf, dn, dnlen, &
& glen, delta, mu, d, work)
integer ldr, n
double precision R(ldr,*), sdiag(*)
double precision qtf(*),dn(*), dnlen, glen, d(*),work(*)
double precision delta, mu

!-----------------------------------------------------------------------------
!
!     Arguments
!
!       Inout  R      Real(ldr,*)     upper triangular matrix R from QR (unaltered)
!                                     strict lower triangle contains
!                                        transposed strict upper triangle of the upper
!                                        triangular matrix S.
!
!       In     n      Integer         order of R.
!
!       In     ldr    Integer         leading dimension of the array R.
!
!       Out    sdiag  Real(*)         vector of size n, containing the
!                                     diagonal elements of the upper
!                                     triangular matrix S.
!
!       In     qtr    Real(*)         trans(Q)*f vector of size n
!       In     dn     Real(*)         Newton step
!       In     dnlen  Real            length Newton step
!       In     glen   Real            length gradient vector
!
!       Inout  mu     Real            Levenberg-Marquardt parameter
!       In     delta  Real            size of trust region (euclidian norm)
!
!       Out    d      Real(*)         vector with step with norm very close to delta
!       Out    work   Real(*)         workspace of size n.
!
!     Description
!
!     determine Levenberg-Marquardt parameter mu such that
!     norm[(R**T R + mu*I)**(-1) * qtf] - delta approximately 0
!     See description in liqrev.f for further details
!
!     Algorithm comes from More: The Levenberg-Marquardt algorithm, Implementation and Theory
!     Lecture Notes in Mathematics, 1978, no. 630.
!     uses liqrev (in file liqrev.f) which is based on Nocedal's method (see comments in file)
!-----------------------------------------------------------------------------

double precision phi, pnorm, qnorm, mulo, muhi,dmu, sqmu
integer iter
logical done


double precision Rone
parameter(Rone=1.0D0)

phi = dnlen - delta
muhi = glen/delta

call dcopy(n,dn,1,d,1)
call dscal(n, Rone/dnlen, d, 1)

!     solve R**T * x = dn
call dtrsv("U","T","N",n,R,ldr,d,1)
qnorm = dnrm2(n,d,1)
mulo = (phi/dnlen)/qnorm**2
mu = mulo

iter = 0
done = .false.
do while( .not. done )
iter = iter + 1
sqmu = sqrt(mu)
call liqrev(n, R, ldr, sqmu, qtf, d, sdiag, work)
pnorm = dnrm2(n,d,1)
call dcopy(n,d,1,work,1)
call dtrstt(R, ldr, n, sdiag, work)
done = abs(pnorm-delta) .le. .1d0*delta .or. iter .gt. 5
if( .not. done ) then
qnorm = dnrm2(n,work,1)
if( pnorm .gt. delta ) then
mulo = max(mulo,mu)
else if( pnorm .lt. delta ) then
muhi = min(muhi,mu)
endif
dmu = (pnorm-delta)/delta * (pnorm/qnorm)**2
mu = max(mulo, mu + dmu)
endif
enddo
return
end

! ============================================================================
! From converted legacy source: liqrev.f90
! ============================================================================

!-----------------------------------------------------------------------------

subroutine liqrev(n,r,ldr,diag,b,x,sdiag,wrk)
integer n,ldr
double precision  r(ldr,*),b(*),x(*),sdiag(*),wrk(*)
double precision diag

!-----------------------------------------------------------------------------
!
!     Arguments
!
!       In     n      Integer         order of R.
!       Inout  R      Real(ldr,*)     upper triangular matrix R from QR
!                                     unaltered
!                                     strict lower triangle contains
!                                        transposed strict upper triangle of the upper
!                                        triangular matrix S.
!
!       In     diag   Real            scalar for matrix D
!
!       In     ldr    Integer         leading dimension of the array R.
!       In     b      Real(*)         vector of size n
!
!       Out    x      Real(*)         vector of size n
!                                     on output contains least squares solution
!                                     of the system R*x = b, D*x = 0.
!
!       Out    sdiag  Real(*)         vector of size n, containing the
!                                     diagonal elements of the upper
!                                     triangular matrix S.
!
!       Out    wrk    Real(*)         workspace of size n.
!
!     Description
!
!     Given an n by n upper triangular matrix R, a diagonal matrix D with positive entries
!     and an n-vector b, determine an x which solves the system
!
!         |R*x| = |b|
!         |D*x| = |0|
!
!     in the least squares sense where D=diag*I.
!     This routine can be used for two different purposes.
!     The first is to provide a method of slightly modifying a singular or ill-conditioned matrix.
!     The second is for calculating a least squares solution to the above problem within
!     the context of e.g. a Levenberg-Marquardt algorithm combined with a More-Hebden algorithm
!     to determine a value of D (diagonal mu) such that x has a predetermined 2-norm.
!
!     The routine could also be used when the matrix R from the QR decomposition of a Jacobian
!     is ill-conditioned (or singular). Then it is difficult to calculate a Newton step
!     accurately (Dennis and Schnabel). D&S advise perturbing trans(J)*J with a positive
!     diagonal matrix.
!
!     The idea is  to solve (J^T * J + mu*I)x=b where mu is a small positive number.
!     Calculation of mu must be done in the calling routine.
!     Using a QR decomposition of J solving this system
!     is equivalent solving (R^T*R + mu*I)x=b, where R comes from the QR decomposition.
!     Solving this system is equivalent to solving the above least squares problem with the
!     elements of the matrix D set to sqrt(mu) which should be done in the calling routine.
!
!     On output the routine also provides an upper triangular matrix S such that
!     (see description of arguments above for the details)
!
!         (trans(R)*R + D*D) = trans(S)*S .
!
!     Method used here is described in
!     Nocedal and Wright, 2006, Numerical Optimization, Springer, ISBN 978-0-387-30303-1
!     page 258--261 (second edition)
!-----------------------------------------------------------------------------

integer j,k
double precision  bj,c,s,sum,temp

double precision Rzero
parameter(Rzero=0.0d0)

!     copy R and b to preserve input and initialise S.
!     Save the diagonal elements of R in wrk.
!     Beware: the algorithm operates on an upper triangular matrix,
!     which is stored in lower triangle of R.
!
do j=1,n
call dcopy(n-j+1,r(j,j),ldr,r(j,j),1)
wrk(j) = r(j,j)
enddo
call dcopy(n,b,1,x,1)

!     eliminate the diagonal matrix D using givens rotations.
!     Nocedal method: start at the bottom right
!     at end of loop R contains diagonal of S
!     save in sdiag and restore original diagonal of R

do j=n,1,-1

!        initialise the row of D to be eliminated

call nuzero(n-j+1,sdiag(j))
sdiag(j) = diag

!        the transformations to eliminate the row of D

bj = Rzero
do k=j,n

!           determine a givens rotation which eliminates the
!           appropriate element in the current row of D.
!           accumulate the transformation in the row of S.

!           eliminate the diagonal element in row j of D
!           this generates fill-in in columns [j+1 .. n] of row j of D
!           successively eliminate the fill-in with givens rotations
!           for R[j+1,j+1] and D[j,j+1].
!           rows of R have been copied into the columns of R initially (see above)
!           perform all operations on those columns to preserve the original R

if (sdiag(k) .ne. Rzero) then

call nuvgiv(r(k,k),sdiag(k),c,s)
if( k .lt. n ) then
call drot(n-k,r(k+1,k),1,sdiag(k+1),1,c,s)
endif

!              compute the modified element of (b,0).

temp =  c*x(k) + s*bj
bj   = -s*x(k) + c*bj
x(k) = temp

endif

enddo

enddo

!     retrieve diagonal of S from diagonal of R
!     restore original diagonal of R

do k=1,n
sdiag(k) = r(k,k)
r(k,k) = wrk(k)
enddo

!     x now contains modified b
!     solve trans(S)*x = x
!     still to be done: guard against division by 0 to be absolutely safe
!     call dblepr('liqrev sdiag', 12, sdiag, n)
x(n) = x(n) / sdiag(n)
do j=n-1,1,-1
sum  = ddot(n-j,r(j+1,j),1,x(j+1),1)
x(j) = (x(j) - sum)/sdiag(j)
enddo

return
end

! ----------------------------------------------------------------------

subroutine dtrstt(S,ldr,n,sdiag,x)
integer ldr, n
double precision S(ldr,*), sdiag(*), x(*)
integer j
double precision sum

!     solve S*x = x where x is the result from subroutine liqrev
!     S is a lower triangular matrix with diagonal entries in sdiag()
!     and here it is in the lower triangular part of R as returned by liqrev

x(1) = x(1) / sdiag(1)
do j=2,n
sum  = ddot(j-1,S(j,1),n,x,1)
x(j) = (x(j) - sum)/sdiag(j)
enddo

return
end

! ============================================================================
! From converted legacy source: liqrup.f90
! ============================================================================
subroutine liqrup(Q,ldq,n,R,ldr,u,v,wk)
integer ldq,n,ldr
double precision Q(ldq,*),R(ldr,*),u(*),v(*),wk(*)

!-----------------------------------------------------------------------------
!
!     Arguments
!
!     Inout  Q       Real(ldq,n)      orthogonal matrix from QR
!     In     ldq     Integer          leading dimension of Q
!     In     n       Integer          order of Q and R
!     Inout  R       Real(ldr,n)      upper triangular matrix R from QR
!     In     ldr     Integer          leading dimension of R
!     In     u       Real(*)          vector u of size n
!     In     v       Real(*)          vector v of size n
!     Out    wk      Real(*)          workspace of size n
!
!     on return
!
!        Q       Q is the matrix with orthonormal columns in a QR
!                decomposition of the matrix B = A + u*v'
!
!        R       R is the upper triangular matrix in a QR
!                decomposition of the matrix B = A + u*v'
!
!     Description
!
!     The matrices Q and R are a QR decomposition of a square matrix
!     A = Q*R.
!     Given Q and R, qrupdt computes a QR decomposition of the rank one
!     modification B = A + u*trans(v) of A. Here u and v are vectors.
!
!     Code based on Reichel and Cragg and simplified and modernized
!           Reichel, L. and W.B. Gragg (1990), Algorithm 686: FORTRAN subroutines for updating the QR decomposition,
!           ACM Trans. Math. Softw., 16, 4, pp. 369--377.
!
!-----------------------------------------------------------------------------

!     Local variables and functions

integer k,i
double precision c, s

!     calculate wk = trans(Q)*u

do i=1,n
wk(i) = ddot(n,Q(1,i),1,u,1)
enddo

!     nuvgiv uses Lapack dlartg and
!     sets its first argument x to value of the last argument returned by dlartg.

!     zero components wk(n),wk(n-1)...wk(2)
!     and apply rotators to R and Q.

do k=n-1,1,-1
call nuvgiv(wk(k),wk(k+1),c,s)
call drot(n-k+1,R(k,k),ldr,R(k+1,k),ldr,c,s)
call drot(n    ,Q(1,k),1  ,Q(1,k+1),1  ,c,s)
enddo

!     r(1,1:n) += wk(1)*v(1:n)
call daxpy(n,wk(1),v,1,R(1,1),ldr)

!     R is of upper hessenberg form. Triangularize R.

do k=1,n-1
call nuvgiv(R(k,k),R(k+1,k),c,s)
call drot(n-k,R(k,k+1),ldr,R(k+1,k+1),ldr,c,s)
call drot(n  ,Q(1,k)  ,1  ,Q(1,k+1)  ,1  ,c,s)
enddo

return
end

! ============================================================================
! From converted legacy source: lirslv.f90
! ============================================================================
!-----------------------------------------------------------------------

subroutine lirslv(R,ldr,n,cndtol, stepadj, &
& qtf,dn,ierr,rcond, rcdwrk,icdwrk)

integer ldr,n,ierr
double precision  cndtol,R(ldr,*)
double precision  dn(*),qtf(*)
double precision  rcdwrk(*)
integer           icdwrk(*)
double precision  rcond
logical           stepadj

!-----------------------------------------------------------------------
!
!     Solve R * dn = -qtf safely with quard against ill-conditioning
!
!     Arguments
!
!     Inout    R       Real(ldr,*)     upper triangular matrix from QR
!                                      if ill conditioned result from liqrev
!                                      Levenberg-Marquardt correction
!                                      Warning: lower triangular part of R destroyed
!     In       ldr     Integer         leading dimension of R
!     In       n       Integer         dimension of problem
!     In       cndtol  Real            tolerance of test for ill conditioning
!     In       stepadj Logical         allow adjusting step for singular/illconditioned jacobian
!     In       qtf     Real(*)         trans(Q)*f()
!     Out      dn      Real(*)         Newton direction
!     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
!                                      1 indicating Jacobian ill-conditioned
!                                      2 indicating Jacobian completely singular
!                                      3 indicating almost zero LM correction
!     Out      rcond   Real            inverse condition of upper triangular R of QR
!     Wk       rcdwrk  Real(*)         workspace
!     Wk       icdwrk  Integer(*)      workspace
!
!-----------------------------------------------------------------------

integer k

double precision Rone
parameter(Rone=1.0d0)
double precision mu

!     check for singularity or ill conditioning

call cndjac(n,R,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)

if( ierr .eq. 0 ) then
!         Normal Newton step
!         solve Jacobian*dn  =  -fn
!         ==> R*dn = - qtf

call dcopy(n,qtf,1,dn,1)
call dtrsv('U','N','N',n,r,ldr,dn,1)
call dscal(n, -Rone, dn, 1)

elseif( stepadj ) then
!         Adjusted Newton step
!         approximately from pseudoinverse(Jac+)
!         use mu to solve (trans(R)*R + mu*I*mu*I) * x = - trans(R) * fn
!         directly from the QR decomposition of R stacked with mu*I
!         a la Levenberg-Marquardt
call compmu(R,ldr,n,mu,rcdwrk,ierr)
if( ierr .eq. 0 ) then
call liqrev(n,R,ldr,mu,qtf,dn, &
& rcdwrk(1+n),rcdwrk(2*n+1))
call dscal(n, -Rone, dn, 1)

!            copy lower triangular R to upper triangular
do k=1,n
call dcopy (n-k+1,R(k,k),1,R(k,k),ldr)
R(k,k) = rcdwrk(1+n+k-1)
enddo
endif
endif

return
end

! ============================================================================
! From converted legacy source: nwbjac.f90
! ============================================================================
subroutine nwbjac(rjac,r,ldr,n,xc,fc,fq,fvec,fjac,epsm,jacflg, &
& wrk1,wrk2,wrk3, &
& xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn, &
& qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)

!-----------------------------------------------------------------------
!
!     Compute Jacobian matrix in xc, fc
!     scale it, compute gradient in xc and generate QR decomposition
!     calculate Newton step
!
!     Arguments
!
!     Out      rjac    Real(ldr,*)     jacobian (n columns) and used for storing full Q from Q
!     Out      r       Real(ldr,*)     used for storing R from QR factorization
!     In       ldr     Integer         leading dimension of rjac
!     In       n       Integer         dimensions of problem
!     In       xc      Real(*)         initial estimate of solution
!     Inout    fc      Real(*)         function values f(xc)
!     Wk       fq      Real(*)         workspace
!     In       fjac    Name            name of routine to calculate jacobian
!                                      (optional)
!     In       fvec    Name            name of routine to calculate f()
!     In       epsm    Real            machine precision
!     In       jacflg  Integer(*)      jacobian flag array
!                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
!                                                  3: user supplied banded
!                                      jacflg[2]: number of sub diagonals or -1 if not banded
!                                      jacflg[3]: number of super diagonals or -1 if not banded
!                                      jacflg[4]: 1 if adjusting step allowed when
!                                                   singular or illconditioned
!     Wk       wrk1    Real(*)         workspace
!     Wk       wrk2    Real(*)         workspace
!     Wk       wrk3    Real(*)         workspace
!     In       xscalm  Integer         x scaling method
!                                        1 from column norms of first jacobian
!                                          increased if needed after first iteration
!                                        0 scaling user supplied
!     Inout    scalex  Real(*)         scaling factors x(*)
!     Out      gp      Real(*)         gradient at xp()
!     In       cndtol  Real            tolerance of test for ill conditioning
!     Wk       rcdwrk  Real(*)         workspace
!     Wk       icdwrk  Integer(*)      workspace
!     Out      dn      Real(*)         Newton step
!     Out      qtf     Real(*)         workspace for nwnstp
!     Out      rcond   Real            estimated inverse condition of R from QR
!     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
!     In       qrwsiz  Integer         size of qrwork
!     Out      njcnt   Integer         number of jacobian evaluations
!     In       iter    Integer         iteration counter (used in scaling)
!     Inout    fstjac  logical         .true. if initial jacobian is available
!                                      on exit set to .false.
!     Out      ierr    Integer         error code
!                                        0 no error
!                                       >0 error in nwnstp (singular ...)
!
!-----------------------------------------------------------------------

integer ldr,n,iter, njcnt, ierr
integer jacflg(*),xscalm,qrwsiz
logical fstjac
double precision  epsm, cndtol, rcond
double precision  rjac(ldr,*),r(ldr,*)
double precision  xc(*),fc(*),dn(*)
double precision  wrk1(*),wrk2(*),wrk3(*)
double precision  qtf(*),gp(*),fq(*)
double precision  scalex(*)
double precision  rcdwrk(*),qrwork(*)
integer           icdwrk(*)
procedure(legacy_fjac) :: fjac
procedure(legacy_fvec) :: fvec

logical stepadj
double precision Rzero, Rone
parameter(Rzero=0.0d0, Rone=1.0d0)

!     evaluate the jacobian at the current iterate xc

if( .not. fstjac ) then
call nwfjac(xc,scalex,fc,fq,n,epsm,jacflg,fvec,fjac,rjac, &
& ldr,wrk1,wrk2,wrk3)
njcnt = njcnt + 1
else
fstjac = .false.
endif

!     if requested calculate x scale from jacobian column norms a la Minpack

if( xscalm .eq. 1 ) then
call vunsc(n,xc,scalex)
call nwcpsx(n,rjac,ldr,scalex,epsm,iter)
call vscal(n,xc,scalex)
endif

call nwscjac(n,rjac,ldr,scalex)

!     evaluate the gradient at the current iterate xc
!     gp = trans(Rjac) * fc
call dgemv('T',n,n,Rone,rjac,ldr,fc,1,Rzero,gp,1)

!     get broyden (newton) step
stepadj = jacflg(4) .eq. 1
call dcopy(n,fc,1,fq,1)
call brdstp(rjac,r,ldr,fq,n,cndtol, stepadj, &
& wrk1,dn,qtf,ierr,rcond, &
& rcdwrk,icdwrk,qrwork,qrwsiz)

!     save some data about jacobian for later output
call nwsnot(0,ierr,rcond)

return
end

!-----------------------------------------------------------------------

subroutine brdstp(rjac,r,ldr,fn,n,cndtol, stepadj, &
& qraux,dn,qtf,ierr,rcond, &
& rcdwrk,icdwrk,qrwork,qrwsiz)

integer ldr,n,ierr,qrwsiz
double precision  cndtol,rjac(ldr,*),r(ldr,*),qraux(*),fn(*)
double precision  dn(*),qtf(*)
double precision  rcdwrk(*),qrwork(*)
integer           icdwrk(*)
double precision  rcond
logical           stepadj

!-----------------------------------------------------------------------
!
!     Calculate the newton step
!
!     Arguments
!
!     Inout    rjac    Real(ldr,*)     jacobian matrix at current iterate; on return full Q
!     Inout    r       Real(ldr,*)     jacobian matrix at current iterate; on return R fom QR
!     In       ldr     Integer         leading dimension of rjac
!     In       fn      Real(*)         function values at current iterate
!     In       n       Integer         dimension of problem
!     In       cndtol  Real            tolerance of test for ill conditioning
!     In       stepadj Logical         allow adjusting step for singular/illconditioned jacobian
!     Inout    qraux   Real(*)         QR info from liqrfa (calling Lapack dgeqrf)
!     Out      dn      Real(*)         Newton direction
!     Out      qtf     Real(*)         trans(Q)*f()
!     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
!                                      1 indicating Jacobian ill-conditioned
!                                      2 indicating Jacobian completely singular
!                                      3 indicating almost zero LM correction
!     Out      rcond   Real            inverse condition of upper triangular R of QR
!     Wk       rcdwrk  Real(*)         workspace
!     Wk       icdwrk  Integer(*)      workspace
!     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
!     In       qrwsiz  Integer         size of qrwork
!
!-----------------------------------------------------------------------

integer info

double precision Rone
parameter(Rone=1.0d0)

!     perform a QR factorization of rjac (simple Lapack routine)
!     check for singularity or ill conditioning
!     form qtf = trans(Q) * fn

call liqrfa(rjac,ldr,n,qraux,qrwork,qrwsiz,ierr)

!     compute qtf = trans(Q)*fn

call dcopy(n,fn,1,qtf,1)
call liqrqt(rjac, ldr, n, qraux, qtf, qrwork, qrwsiz, info)

!     copy the upper triangular part of the QR decomposition
!     contained in Rjac into R[*, 1..n].
!     form Q from the QR decomposition (taur/qraux in wrk1)

call dlacpy('U',n,n,rjac,ldr,r,ldr)
call liqrqq(rjac,ldr,qraux,n,qrwork,qrwsiz,info)

!     now Rjac[* ,1..n] holds expanded Q
!     now R[* ,1..n] holds full upper triangle R

call lirslv(R,ldr,n,cndtol, stepadj, &
& qtf,dn,ierr,rcond, rcdwrk,icdwrk)

return
end

! ============================================================================
! From converted legacy source: nwbrdn.f90
! ============================================================================

subroutine brsolv(ldr,xc,n,scalex,maxit, &
& jacflg,xtol,ftol,btol,cndtol,global,xscalm, &
& stepmx,delta,sigma, &
& rjac,r,wrk1,wrk2,wrk3,wrk4,fc,fq,dn,d,qtf, &
& rcdwrk,icdwrk,qrwork,qrwsiz,epsm, &
& fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,iter, &
& termcd)

integer ldr,n,termcd,njcnt,nfcnt,iter
integer maxit,jacflg(*),global,xscalm,qrwsiz
integer outopt(*)
double precision  xtol,ftol,btol,cndtol
double precision  stepmx,delta,sigma,fpnorm,epsm
double precision  rjac(ldr,*),r(ldr,*)
double precision  xc(*),fc(*),xp(*),fp(*),dn(*),d(*)
double precision  wrk1(*),wrk2(*),wrk3(*),wrk4(*)
double precision  qtf(*),gp(*),fq(*)
double precision  scalex(*)
double precision  rcdwrk(*),qrwork(*)
integer           icdwrk(*)
procedure(legacy_fjac) :: fjac
procedure(legacy_fvec) :: fvec

!-----------------------------------------------------------------------
!
!     Solve system of nonlinear equations with Broyden and global strategy
!
!
!     Arguments
!
!     In       ldr     Integer         leading dimension of rjac
!     In       xc      Real(*)         initial estimate of solution
!     In       n       Integer         dimensions of problem
!     Inout    scalex  Real(*)         scaling factors x(*)
!     In       maxit   Integer         maximum number of allowable iterations
!     In       jacflg  Integer(*)      jacobian flag array
!                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
!                                                  3: user supplied banded
!                                      jacflg[2]: number of sub diagonals or -1 if not banded
!                                      jacflg[3]: number of super diagonals or -1 if not banded
!                                      jacflg[4]: 1 if adjusting step allowed when
!                                                   singular or illconditioned
!     In       xtol    Real            tolerance at which successive iterates x()
!                                      are considered close enough to
!                                      terminate algorithm
!     In       ftol    Real            tolerance at which function values f()
!                                      are considered close enough to zero
!     Inout    btol    Real            x tolerance for backtracking
!     Inout    cndtol  Real            tolerance of test for ill conditioning
!     In       global  Integer         global strategy to use
!                                        1 cubic linesearch
!                                        2 quadratic linesearch
!                                        3 geometric linesearch
!                                        4 double dogleg
!                                        5 powell dogleg
!                                        6 hookstep (More-Hebden Levenberg-Marquardt)
!     In       xscalm  Integer         x scaling method
!                                        1 from column norms of first jacobian
!                                          increased if needed after first iteration
!                                        0 scaling user supplied
!     In       stepmx  Real            maximum allowable step size
!     In       delta   Real            trust region radius
!     In       sigma   Real            reduction factor geometric linesearch
!     Inout    rjac    Real(ldr,*)     jacobian (n columns)(compact QR decomposition/Q matrix)
!     Inout    r       Real(ldr,*)     stored R from QR decomposition
!     Wk       wrk1    Real(*)         workspace
!     Wk       wrk2    Real(*)         workspace
!     Wk       wrk3    Real(*)         workspace
!     Wk       wrk4    Real(*)         workspace
!     Inout    fc      Real(*)         function values f(xc)
!     Wk       fq      Real(*)         workspace
!     Wk       dn      Real(*)         workspace
!     Wk       d       Real(*)         workspace
!     Wk       qtf     Real(*)         workspace
!     Wk       rcdwrk  Real(*)         workspace
!     Wk       icdwrk  Integer(*)      workspace
!     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
!     In       qrwsiz  Integer         size of qrwork
!     In       epsm    Real            machine precision
!     In       fjac    Name            name of routine to calculate jacobian
!                                      (optional)
!     In       fvec    Name            name of routine to calculate f()
!     In       outopt  Integer(*)      output options
!     Out      xp      Real(*)         final x()
!     Out      fp      Real(*)         final f(xp)
!     Out      gp      Real(*)         gradient at xp()
!     Out      njcnt   Integer         number of jacobian evaluations
!     Out      nfcnt   Integer         number of function evaluations
!     Out      iter    Integer         number of (outer) iterations
!     Out      termcd  Integer         termination code
!
!-----------------------------------------------------------------------

integer gcnt,retcd,ierr
double precision  dum(2),dlt0,fcnorm,rcond
logical fstjac
logical jacevl,jacupd
logical stepadj
integer priter



double precision Rone
parameter(Rone=1.0d0)

!     initialization

retcd = 0
iter  = 0
njcnt = 0
nfcnt = 0
ierr  = 0

dum(1) = 0
dlt0 = delta

if( outopt(1) .eq. 1 ) then
priter = 1
else
priter = -1
endif

!     evaluate function

call vscal(n,xc,scalex)
call nwfvec(xc,n,scalex,fvec,fc,fcnorm,wrk1)

!     evaluate user supplied or finite difference jacobian and check user supplied
!     jacobian, if requested

fstjac = .false.
if(mod(jacflg(1),2) .eq. 1) then

if( outopt(2) .eq. 1 ) then
fstjac = .true.
njcnt = njcnt + 1
call nwfjac(xc,scalex,fc,fq,n,epsm,jacflg,fvec,fjac,rjac, &
& ldr,wrk1,wrk2,wrk3)
call chkjac(rjac,ldr,xc,fc,n,epsm,jacflg,scalex, &
& fq,wrk1,wrk2,fvec,termcd)
if(termcd .lt. 0) then
!              copy initial values
call dcopy(n,xc,1,xp,1)
call dcopy(n,fc,1,fp,1)
call vunsc(n,xp,scalex)
fpnorm = fcnorm
return
endif
endif

endif

!     check stopping criteria for input xc

call nwtcvg(xc,fc,xc,xtol,retcd,ftol,iter,maxit,n,ierr,termcd)

if(termcd .gt. 0) then
call dcopy(n,xc,1,xp,1)
call dcopy(n,fc,1,fp,1)
fpnorm = fcnorm
if( outopt(3) .eq. 1 .and. .not. fstjac ) then
njcnt = njcnt + 1
call nwfjac(xp,scalex,fp,fq,n,epsm,jacflg,fvec,fjac,rjac, &
& ldr,wrk1,wrk2,wrk3)
endif
return
endif

if( priter .gt. 0 ) then

dum(1) = fcnorm
dum(2) = abs(fc(idamax(n,fc,1)))

if( global .eq. 0 ) then
call nwprot(iter, -1, dum)
elseif( global .le. 3 ) then
call nwlsot(iter,-1,dum)
elseif( global .eq. 4 ) then
call nwdgot(iter,-1,0,dum)
elseif( global .eq. 5 ) then
call nwpwot(iter,-1,0,dum)
elseif( global .eq. 6 ) then
call nwmhot(iter,-1,0,dum)
endif

endif

jacevl  = .true.
stepadj = jacflg(4) .eq. 1

do while( termcd .eq. 0 )
iter = iter+1

if( jacevl ) then

call nwbjac(rjac,r,ldr,n,xc,fc,fq,fvec,fjac,epsm,jacflg, &
& wrk1,wrk2,wrk3, &
& xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn, &
& qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)

else

!          - get broyden step
!          - calculate approximate gradient

call dcopy(n,fc,1,fq,1)
call brodir(rjac,ldr,r,fq,n,cndtol, stepadj, &
& dn,qtf,ierr,rcond,rcdwrk,icdwrk)

if( ierr .eq. 0 ) then
call dcopy(n,qtf,1,gp,1)
call dtrmv('U','T','N',n,r,ldr,gp,1)
endif
endif
!      - choose the next iterate xp by a global strategy

if( ierr .gt. 0 ) then
!           jacobian singular or too ill-conditioned
call nweset(n,xc,fc,fcnorm,xp,fp,fpnorm,gcnt,priter,iter)
elseif(global .eq. 0) then
call nwpure(n,xc,dn,stepmx,scalex, &
& fvec,xp,fp,fpnorm,wrk1,retcd,gcnt, &
& priter,iter)
elseif(global .eq. 1) then
call nwclsh(n,xc,fcnorm,dn,gp,stepmx,btol,scalex, &
& fvec,xp,fp,fpnorm,wrk1,retcd,gcnt, &
& priter,iter)
elseif(global .eq. 2) then
call nwqlsh(n,xc,fcnorm,dn,gp,stepmx,btol,scalex, &
& fvec,xp,fp,fpnorm,wrk1,retcd,gcnt, &
& priter,iter)
elseif(global .eq. 3) then
call nwglsh(n,xc,fcnorm,dn,gp,sigma,stepmx,btol,scalex, &
& fvec,xp,fp,fpnorm,wrk1,retcd,gcnt, &
& priter,iter)
elseif(global .eq. 4) then
call nwddlg(n,r,ldr,dn,gp,xc,fcnorm,stepmx, &
& btol,delta,qtf,scalex, &
& fvec,d,fq,wrk1,wrk2,wrk3,wrk4, &
& xp,fp,fpnorm,retcd,gcnt,priter,iter)
elseif(global .eq. 5) then
call nwpdlg(n,r,ldr,dn,gp,xc,fcnorm,stepmx, &
& btol,delta,qtf,scalex, &
& fvec,d,fq,wrk1,wrk2,wrk3,wrk4, &
& xp,fp,fpnorm,retcd,gcnt,priter,iter)
elseif(global .eq. 6) then
call nwmhlm(n,r,ldr,dn,gp,xc,fcnorm,stepmx, &
& btol,delta,qtf,scalex, &
& fvec,d,fq,wrk1,wrk2,wrk3,wrk4, &
& xp,fp,fpnorm,retcd,gcnt,priter,iter)
endif

nfcnt = nfcnt + gcnt

!      - check stopping criteria for the new iterate xp

call nwtcvg(xp,fp,xc,xtol,retcd,ftol,iter,maxit,n,ierr,termcd)

if( termcd .eq. 3 .and. .not. jacevl ) then
!           global strategy failed but jacobian is out of date
!           try again with proper jacobian
!           reset trust region radius

jacevl = .true.
jacupd = .false.
delta = dlt0
termcd = 0

elseif(termcd .gt. 0) then
jacupd = .false.
else
jacupd = .true.
jacevl = .false.
endif

if( jacupd ) then
!           perform Broyden update of current jacobian
!           update xc, fc, and fcnorm
call brupdt(n,rjac,r,ldr,xc,xp,fc,fp,epsm, &
& wrk1,wrk2,wrk3)
call dcopy(n,xp,1,xc,1)
call dcopy(n,fp,1,fc,1)
fcnorm = fpnorm
endif

enddo

if( outopt(3) .eq. 1 ) then
!        final update of jacobian
call brupdt(n,rjac,r,ldr,xc,xp,fc,fp,epsm, &
& wrk1,wrk2,wrk3)
!        reconstruct Broyden matrix
!        calculate Q * R where Q is overwritten by result
!        Q is in rjac and R is in r
call dtrmm('R','U','N','N',n,n,Rone,r,n,rjac,n)
!        unscale
call nwunscjac(n,rjac,ldr,scalex)
endif

call vunsc(n,xp,scalex)

return
end

!-----------------------------------------------------------------------

subroutine brupdt(n,q,r,ldr,xc,xp,fc,fp,epsm,dx,df,wa)
integer n,ldr
double precision  q(ldr,*),r(ldr,*)
double precision  xc(*),xp(*),fc(*),fp(*),dx(*),df(*),wa(*)
double precision  epsm

!-----------------------------------------------------------------------
!
!     Calculate new Q and R from rank-1 update with xp-xc and fp-fc
!     using Broyden method
!
!     Arguments
!
!     In       n       Integer         size of xc() etc.
!     Inout    Q       Real(ldr,n)     orthogonal matrix Q from QR
!                                       On output updated Q
!     Inout    R       Real(ldr,n)     upper triangular R from QR
!                                       On output updated R
!     In       ldr     Integer         leading dimension of Q and R
!     In       xc      Real(*)         current x() values
!     In       xp      Real(*)         new     x() values
!     In       fc      Real(*)         current f(xc)
!     In       fp      Real(*)         new     f(xp)
!     In       epsm    Real            machine precision
!     Wk       dx      Real(*)         workspace
!     Wk       df      Real(*)         workspace
!     Wk       wa      Real(*)         workspace
!
!-----------------------------------------------------------------------

integer i
double precision  eta,sts

logical doupdt

double precision Rzero, Rone, Rtwo, Rhund
parameter(Rzero=0.0d0, Rone=1.0d0, Rtwo=2.0d0, Rhund=100d0)

eta    = Rhund * Rtwo * epsm
doupdt = .false.

do i=1,n
dx(i) = xp(i) - xc(i)
df(i) = fp(i) - fc(i)
enddo

!     clear lower triangle

do i=1,n-1
call nuzero(n-i,r(i+1,i))
enddo

!     calculate df - B*dx = df - Q*R*dx
!     wa = R*dx
!     df = df - Q*(R*dx) (!not really needed if qrupdt were to be changed)
!     do not update with noise

call dcopy(n,dx,1,wa,1)
call dtrmv('U','N','N',n,r,ldr,wa,1)
call dgemv('N',n,n,-Rone,q,ldr,wa,1,Rone,df,1)

do i=1,n
if( abs(df(i)) .gt. eta*( abs(fp(i)) + abs(fc(i)) ) ) then
doupdt = .true.
else
df(i)  = Rzero
endif
enddo

if( doupdt ) then
!        equation 8.3.1 from Dennis and Schnabel (page 187)(Siam edition)
sts = dnrm2(n,dx,1)
call dscal(n,Rone/sts,dx,1)
call dscal(n,Rone/sts,df,1)
call liqrup(q,ldr,n,r,ldr,df,dx,wa)
endif

return
end

!-----------------------------------------------------------------------

subroutine brodir(q,ldr,r,fn,n,cndtol,stepadj,dn,qtf, &
& ierr,rcond,rcdwrk,icdwrk)

integer ldr,n,ierr
double precision  cndtol,q(ldr,*),r(ldr,*),fn(*)
double precision  dn(*),qtf(*)
double precision  rcdwrk(*)
integer           icdwrk(*)
double precision  rcond
logical           stepadj

!-----------------------------------------------------------------------
!
!     Calculate the approximate newton direction
!
!     Arguments
!
!     Inout    Q       Real(ldr,*)     Q part from QR at current iterate
!     In       ldr     Integer         leading dimension of Q and R
!     In       R       Real(ldr,*)     upper triangular R from QR decomposition
!     In       fn      Real(*)         function values at current iterate
!     In       n       Integer         dimension of problem
!     In       cndtol  Real            tolerance of test for ill conditioning
!     In       stepadj Logical         allow adjusting step for singular/illconditioned jacobian
!     Out      dn      Real(*)         Newton direction
!     Out      qtf     Real(*)         trans(Q)*f()
!     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
!                                      1 indicating Jacobian ill-conditioned
!                                      2 indicating Jacobian completely singular
!                                      3 indicating almost zero LM correction
!     Out      rcond   Real            inverse condition of matrix
!     Wk       rcdwrk  Real(*)         workspace
!     Wk       icdwrk  Integer(*)      workspace
!
!     QR decomposition with no pivoting.
!
!-----------------------------------------------------------------------

double precision Rzero, Rone
parameter(Rzero=0.0d0, Rone=1.0d0)

!     form qtf = trans(Q) * fn
call dgemv('T',n,n,Rone,q,ldr,fn,1,Rzero,qtf,1)

call lirslv(R,ldr,n,cndtol, stepadj, &
& qtf,dn,ierr,rcond, rcdwrk,icdwrk)
call nwsnot(1,ierr,rcond)

return
end

! ============================================================================
! From converted legacy source: nwclsh.f90
! ============================================================================

subroutine nwclsh(n,xc,fcnorm,d,g,stepmx,xtol,scalex,fvec, &
& xp,fp,fpnorm,xw,retcd,gcnt,priter,iter)

integer n,retcd,gcnt
double precision  stepmx,xtol,fcnorm,fpnorm
double precision  xc(*)
double precision  d(*),g(*),xp(*),fp(*),xw(*)
double precision  scalex(*)
procedure(legacy_fvec) :: fvec

integer priter,iter

!-------------------------------------------------------------------------
!
!     Find a next acceptable iterate using a safeguarded cubic line search
!     along the newton direction
!
!     Arguments
!
!     In       n       Integer          dimension of problem
!     In       xc      Real(*)          current iterate
!     In       fcnorm  Real             0.5 * || f(xc) ||**2
!     In       d       Real(*)          newton direction
!     In       g       Real(*)          gradient at current iterate
!     In       stepmx  Real             maximum stepsize
!     In       xtol    Real             relative step size at which
!                                       successive iterates are considered
!                                       close enough to terminate algorithm
!     In       scalex  Real(*)          diagonal scaling matrix for x()
!     In       fvec    Name             name of routine to calculate f()
!     In       xp      Real(*)          new x()
!     In       fp      Real(*)          new f(x)
!     In       fpnorm  Real             .5*||fp||**2
!     Out      xw      Real(*)           workspace for unscaling x(*)
!
!     Out      retcd   Integer          return code
!                                         0 new satisfactory x() found
!                                         1 no  satisfactory x() found
!                                           sufficiently distinct from xc()
!
!     Out      gcnt    Integer          number of steps taken
!     In       priter  Integer           >0 unit if intermediate steps to be printed
!                                        -1 if no printing
!
!-------------------------------------------------------------------------

integer i
double precision  alpha,slope,rsclen,oarg(4)
double precision  lambda,lamhi,lamlo,t
double precision ftarg
double precision  dlen
double precision a, b, disc, fpt, fpt0, fpnorm0, lambda0

logical firstback

parameter (alpha = 1.0d-4)

double precision Rhalf, Rone, Rtwo, Rthree, Rten, Rzero
parameter(Rzero=0.0d0)
parameter(Rhalf=0.5d0, Rone=1.0d0, Rtwo=2.0d0, Rten=10.0d0)
parameter(Rthree=3.0d0)
double precision t1,t2


!     silence warnings issued by ftncheck

lambda0 = Rzero
fpnorm0 = Rzero

!     safeguard initial step size

dlen = dnrm2(n,d,1)
if( dlen .gt. stepmx ) then
lamhi  = stepmx / dlen
else
lamhi  = Rone
endif

!     compute slope  =  g-trans * d

slope = ddot(n,g,1,d,1)

!     compute the smallest value allowable for the damping
!     parameter lambda ==> lamlo

rsclen = nudnrm(n,d,xc)
lamlo  = xtol / rsclen

!     initialization of retcd and lambda (linesearch length)

retcd  = 2
lambda = lamhi
gcnt   = 0
firstback = .true.

do while( retcd .eq. 2 )

!        compute next x

do i=1,n
xp(i) = xc(i) + lambda*d(i)
enddo

!        evaluate functions and the objective function at xp

call nwfvec(xp,n,scalex,fvec,fp,fpnorm,xw)
gcnt = gcnt + 1
ftarg = fcnorm + alpha * lambda * slope

if( priter .gt. 0) then
oarg(1) = lambda
oarg(2) = ftarg
oarg(3) = fpnorm
oarg(4) = abs(fp(idamax(n,fp,1)))
call nwlsot(iter,1,oarg)
endif

!        first is quadratic
!        test whether the standard step produces enough decrease
!        of the objective function.
!        If not update lambda and compute a new next iterate

if( fpnorm .le. ftarg ) then
retcd = 0
else
if( fpnorm .gt. lamlo**2 * sqrt(dlamch('O')) ) then
!               safety against overflow in what follows (use lamlo**2 for safety)
lambda = lambda/Rten
firstback = .true.
else
if( firstback ) then
t = ((-lambda**2)*slope/Rtwo)/ &
& (fpnorm-fcnorm-lambda*slope)
firstback = .false.
else
fpt  = fpnorm -fcnorm - lambda*slope
fpt0 = fpnorm0-fcnorm - lambda0*slope
a = fpt/lambda**2 - fpt0/lambda0**2
b = -lambda0*fpt/lambda**2 + lambda*fpt0/lambda0**2
a = a /(lambda - lambda0)
b = b /(lambda - lambda0)
if( abs(a) .le. dlamch('E') ) then
!                      not quadratic but linear
t = -slope/(2*b)
else
!                      use Higham procedure to compute roots acccurately
!                      Higham: Accuracy and Stability of Numerical Algorithms, second edition,2002, page 10.
!                      Actually solving 3*a*x^2+2*b*x+c=0 ==> (3/2)*a*x^2+b*x+(c/2)=0
!                      use max to prevent NaN in sqrt(disc)
disc = max(b**2 - Rthree * a * slope,Rzero)
t1 = -(b+sign(Rone,b)*sqrt(disc))/(Rthree*a)
t2 = slope/(Rthree*a)/t1
if(a .gt. Rzero ) then
!                          upward opening parabola ==> rightmost is solution
t = max(t1,t2)
else
!                          downward opening parabola ==> leftmost is solution
t = min(t1,t2)
endif
endif
t = min(t, Rhalf*lambda)
endif
lambda0 = lambda
fpnorm0 = fpnorm
lambda = max(t,lambda/Rten)
if(lambda .lt. lamlo) then
retcd = 1
endif
endif
endif
enddo

return
end

! ============================================================================
! From converted legacy source: nwddlg.f90
! ============================================================================

subroutine nwddlg(n,rjac,ldr,dn,g,xc,fcnorm,stepmx,xtol, &
& delta,qtf,scalex,fvec,d,xprev, &
& ssd,v,wa,fprev,xp,fp,fpnorm,retcd,gcnt, &
& priter,iter)

integer ldr, n, retcd, gcnt, priter, iter
double precision  fcnorm, stepmx, xtol, fpnorm, delta
double precision  rjac(ldr,*), dn(*), g(*), xc(*), qtf(*)
double precision  scalex(*), d(*)
double precision  xprev(*), xp(*), fp(*)
double precision  ssd(*), v(*), wa(*), fprev(*)
procedure(legacy_fvec) :: fvec

!-------------------------------------------------------------------------
!
!     Find a next iterate xp by the double dogleg method
!
!     Arguments
!
!     In       n       Integer         size of problem: dimension x, f
!     In       Rjac    Real(ldr,*)     R of QR-factored jacobian
!     In       ldr     Integer         leading dimension of Rjac
!     Inout    dn      Real(*)         newton direction
!     Inout    g       Real(*)         gradient at current point
!                                      trans(jac)*f()
!     In       xc      Real(*)         current iterate
!     In       fcnorm  Real            .5*||f(xc)||**2
!     In       stepmx  Real            maximum stepsize
!     In       xtol    Real            x-tolerance (stepsize)
!     Inout    delta   Real            on input: initial trust region radius
!                                                if -1 then set to something
!                                                reasonable
!                                      on output: final value
!                                      ! Do not modify between calls while
!                                        still iterating
!     In       qtf     Real(*)         trans(Q)*f(xc)
!     In       scalex  Real(*)         scaling factors for x()
!     In       fvec    Name            name of subroutine to evaluate f(x)
!                                      ! must be declared external in caller
!     Wk       d       Real(*)         work vector
!     Wk       xprev   Real(*)         work vector
!     Wk       ssd     Real(*)         work vector
!     Wk       v       Real(*)         work vector
!     Wk       wa      Real(*)         work vector
!     Wk       fprev   Real(*)         work vector
!     Out      xp      Real(*)         new x()
!     Out      fp      Real(*)         new f(xp)
!     Out      fpnorm  Real            new .5*||f(xp)||**2
!     Out      retcd   Integer         return code
!                                       0  new satisfactory x() found
!                                       1  no  satisfactory x() found
!     Out      gcnt    Integer         number of steps taken
!     In       priter  Integer         print flag
!                                       -1 no intermediate printing
!                                       >0 yes for print of intermediate results
!     In       iter    Integer         current iteration (only used for above)
!
!     All vectors at least size n
!
!-------------------------------------------------------------------------

integer i
double precision  dnlen,ssdlen,alpha,beta,lambda,fpred
double precision  sqalpha,eta,gamma,fpnsav,oarg(7)

logical nwtstep
integer dtype



double precision Rone, Rtwo, Rten, Rhalf, Rp2, Rp8
parameter(Rhalf=0.5d0)
parameter(Rone=1.0d0, Rtwo=2.0d0, Rten=10.0d0)
parameter(Rp2 = Rtwo/Rten, Rp8 = Rone - Rp2)
double precision Rzero
parameter(Rzero=0.0d0)

!     length newton direction

dnlen = dnrm2(n, dn, 1)

!     steepest descent direction and length

sqalpha = dnrm2(n,g,1)
alpha   = sqalpha**2

call dcopy(n, g, 1, d, 1)
call dtrmv('U','N','N',n,rjac,ldr,d,1)
beta = dnrm2(n,d,1)**2

call dcopy(n, g, 1, ssd, 1)
call dscal(n, -(alpha/beta), ssd, 1)

ssdlen = alpha*sqalpha/beta

!     set trust radius to ssdlen or dnlen if required

if( delta .eq. -Rone ) then
delta = min(ssdlen, stepmx)
elseif( delta .eq. -Rtwo ) then
delta = min(dnlen, stepmx)
endif

!     calculate double dogleg parameter

gamma = alpha*alpha/(-beta*ddot(n,g,1,dn,1))
!      call dgdbg(gamma, alpha*alpha, -beta*ddot(n,g,1,dn,1))
!     precautionary (just in case)
eta = max(Rzero, min(Rone,Rp2 + Rp8*gamma))

retcd = 4
gcnt  = 0

do while( retcd .gt. 1 )
!        find new step by double dogleg algorithm

call ddlgstp(n,dn,dnlen,delta,v, &
& ssd,ssdlen,eta,d,dtype,lambda)
nwtstep = dtype .eq. 4
!        compute the model prediction 0.5*||F + J*d||**2 (L2-norm)

call dcopy(n,d,1,wa,1)
call dtrmv('U','N','N',n,rjac,ldr,wa,1)
call daxpy(n, Rone, qtf,1,wa,1)
fpred = Rhalf * dnrm2(n,wa,1)**2

!        evaluate function at xp = xc + d

do i=1,n
xp(i) = xc(i) + d(i)
enddo

call nwfvec(xp,n,scalex,fvec,fp,fpnorm,wa)
gcnt = gcnt + 1

!        check whether the global step is acceptable

oarg(2) = delta
call nwtrup(n,fcnorm,g,d,nwtstep,stepmx,xtol,delta, &
& fpred,retcd,xprev,fpnsav,fprev,xp,fp,fpnorm)

if( priter .gt. 0 ) then
oarg(1) = lambda
oarg(3) = delta
oarg(4) = eta
oarg(5) = fpnorm
oarg(6) = abs(fp(idamax(n,fp,1)))
call nwdgot(iter,dtype,retcd,oarg)
endif

enddo

return
end

!-----------------------------------------------------------------------

subroutine ddlgstp(n,dn,dnlen,delta,v, &
& ssd,ssdlen,eta,d,dtype,lambda)
integer n
double precision  dn(*), ssd(*), v(*), d(*)
double precision  dnlen, delta, ssdlen, eta, lambda
integer dtype

!-------------------------------------------------------------------------
!
!     Find a new step by the double dogleg algorithm
!     Internal routine for nwddlg
!
!     Arguments
!
!     In       n       Integer         size of problem
!     In       dn      Real(*)         current newton step
!     Out      dnlen   Real            length dn()
!     In       delta   Real            current trust region radius
!     Out      v       Real(*)         (internal) eta * dn() - ssd()
!     In       ssd     Real(*)         (internal) steepest descent direction
!     In       ssdlen  Real            (internal) length ssd
!     In       eta     Real            (internal) double dogleg parameter
!     Out      d       Real(*)         new step for x()
!     Out      dtype   Integer         steptype
!                                       1 steepest descent
!                                       2 combination of dn and ssd
!                                       3 partial newton step
!                                       4 full newton direction
!     Out      lambda  Real            weight of eta*dn() in d()
!                                      closer to 1 ==> more of eta*dn()
!
!-----------------------------------------------------------------------

integer i
double precision vssd, vlen


if(dnlen .le. delta) then

!        Newton step smaller than trust radius ==> take it

call dcopy(n, dn, 1, d, 1)
delta = dnlen
dtype = 4

elseif(eta*dnlen .le. delta) then

!        take partial step in newton direction

call dcopy(n, dn, 1, d, 1)
call dscal(n, delta / dnlen, d, 1)
dtype = 3

elseif(ssdlen .ge. delta) then

!        take step in steepest descent direction

call dcopy(n, ssd, 1, d, 1)
call dscal(n, delta / ssdlen, d, 1)
dtype = 1

else

!        calculate convex combination of ssd and eta*dn with length delta

do i=1,n
v(i) = eta*dn(i) - ssd(i)
enddo

vssd = ddot(n,v,1,ssd,1)
vlen = dnrm2(n,v,1)**2

lambda =(-vssd+sqrt(vssd**2-vlen*(ssdlen**2-delta**2)))/vlen
call dcopy(n, ssd, 1, d, 1)
call daxpy(n, lambda, v, 1, d, 1)
dtype = 2

endif

return
end

! ============================================================================
! From converted legacy source: nwglsh.f90
! ============================================================================

subroutine nwglsh(n,xc,fcnorm,d,g,sigma,stepmx,xtol,scalex,fvec, &
& xp,fp,fpnorm,xw,retcd,gcnt,priter,iter)

integer n,retcd,gcnt
double precision  sigma,stepmx,xtol,fcnorm,fpnorm
double precision  xc(*)
double precision  d(*),g(*),xp(*),fp(*),xw(*)
double precision  scalex(*)
procedure(legacy_fvec) :: fvec

integer priter,iter

!-------------------------------------------------------------------------
!
!     Find a next acceptable iterate using geometric line search
!     along the newton direction
!
!     Arguments
!
!     In       n       Integer          dimension of problem
!     In       xc      Real(*)          current iterate
!     In       fcnorm  Real             0.5 * || f(xc) ||**2
!     In       d       Real(*)          newton direction
!     In       g       Real(*)          gradient at current iterate
!     In       sigma   Real             reduction factor for lambda
!     In       stepmx  Real             maximum stepsize
!     In       xtol    Real             relative step size at which
!                                       successive iterates are considered
!                                       close enough to terminate algorithm
!     In       scalex  Real(*)          diagonal scaling matrix for x()
!     In       fvec    Name             name of routine to calculate f()
!     In       xp      Real(*)          new x()
!     In       fp      Real(*)          new f(x)
!     In       fpnorm  Real             .5*||fp||**2
!     Out      xw      Real(*)          workspace for unscaling x
!
!     Out      retcd   Integer          return code
!                                         0 new satisfactory x() found
!                                         1 no  satisfactory x() found
!                                           sufficiently distinct from xc()
!
!     Out      gcnt    Integer          number of steps taken
!     In       priter  Integer           >0 if intermediate steps to be printed
!                                        -1 if no printing
!
!-------------------------------------------------------------------------

integer i
double precision  alpha,slope,rsclen,oarg(4)
double precision  lambda,lamhi,lamlo
double precision ftarg
double precision  dlen



parameter (alpha = 1.0d-4)

double precision Rone
parameter(Rone=1.0d0)

!     safeguard initial step size

dlen = dnrm2(n,d,1)
if( dlen .gt. stepmx ) then
lamhi  = stepmx / dlen
else
lamhi  = Rone
endif

!     compute slope  =  g-trans * d

slope = ddot(n,g,1,d,1)

!     compute the smallest value allowable for the damping
!     parameter lambda ==> lamlo

rsclen = nudnrm(n,d,xc)
lamlo  = xtol / rsclen

!     initialization of retcd and lambda (linesearch length)

retcd  = 2
lambda = lamhi
gcnt   = 0

do while( retcd .eq. 2 )

!        compute next x

do i=1,n
xp(i) = xc(i) + lambda*d(i)
enddo

!        evaluate functions and the objective function at xp

call nwfvec(xp,n,scalex,fvec,fp,fpnorm,xw)
gcnt = gcnt + 1
ftarg = fcnorm + alpha * lambda * slope

if( priter .gt. 0) then
oarg(1) = lambda
oarg(2) = ftarg
oarg(3) = fpnorm
oarg(4) = abs(fp(idamax(n,fp,1)))
call nwlsot(iter,1,oarg)
endif

!        test if the standard step produces enough decrease
!        of the objective function.
!        If not update lambda and compute a new next iterate

if( fpnorm .le. ftarg ) then
retcd = 0
else
lambda  = sigma * lambda
if(lambda .lt. lamlo) then
retcd = 1
endif
endif

enddo

return
end

! ============================================================================
! From converted legacy source: nwmhlm.f90
! ============================================================================

subroutine nwmhlm(n,rjac,ldr,dn,g,xc,fcnorm,stepmx,xtol, &
& delta,qtf,scalex,fvec,d,xprev, &
& ssd,v,wa,fprev,xp,fp,fpnorm,retcd,gcnt, &
& priter,iter)

integer ldr, n, retcd, gcnt, priter, iter
double precision  fcnorm, stepmx, xtol, fpnorm, delta
double precision  rjac(ldr,*), dn(*), g(*), xc(*), qtf(*)
double precision  scalex(*), d(*)
double precision  xprev(*), xp(*), fp(*)
double precision  ssd(*), v(*), wa(*), fprev(*)
procedure(legacy_fvec) :: fvec

!-------------------------------------------------------------------------
!
!     Find a next iterate xp by the More-Hebden-Levenberg-Marquardt method
!
!     Arguments
!
!     In       n       Integer         size of problem: dimension x, f
!     In       Rjac    Real(ldr,*)     R of QR-factored jacobian
!     In       ldr     Integer         leading dimension of Rjac
!     Inout    dn      Real(*)         newton direction
!     Inout    g       Real(*)         gradient at current point
!                                      trans(jac)*f()
!     In       xc      Real(*)         current iterate
!     In       fcnorm  Real            .5*||f(xc)||**2
!     In       stepmx  Real            maximum stepsize
!     In       xtol    Real            x-tolerance (stepsize)
!     Inout    delta   Real            on input: initial trust region radius
!                                                if -1 then set to something
!                                                reasonable
!                                      on output: final value
!                                      ! Do not modify between calls while
!                                        still iterating
!     In       qtf     Real(*)         trans(Q)*f(xc)
!     In       scalex  Real(*)         scaling factors for x()
!     In       fvec    Name            name of subroutine to evaluate f(x)
!                                      ! must be declared external in caller
!     Wk       d       Real(*)         work vector
!     Wk       xprev   Real(*)         work vector
!     Wk       ssd     Real(*)         work vector
!     Wk       v       Real(*)         work vector
!     Wk       wa      Real(*)         work vector
!     Wk       fprev   Real(*)         work vector
!     Out      xp      Real(*)         new x()
!     Out      fp      Real(*)         new f(xp)
!     Out      fpnorm  Real            new .5*||f(xp)||**2
!     Out      retcd   Integer         return code
!                                       0  new satisfactory x() found
!                                       1  no  satisfactory x() found
!     Out      gcnt    Integer         number of steps taken
!     In       priter  Integer         print flag
!                                       -1 no intermediate printing
!                                       >0 yes for print of intermediate results
!     In       iter    Integer         current iteration (only used for above)
!
!     All vectors at least size n
!
!-------------------------------------------------------------------------

integer i
double precision  dnlen,glen,ssdlen,alpha,beta,mu,fpred
double precision  fpnsav,oarg(6)

logical nwtstep
integer dtype



double precision Rone, Rtwo, Rhalf
parameter(Rhalf=0.5d0)
parameter(Rone=1.0d0, Rtwo=2.0d0)

!     length newton direction

dnlen = dnrm2(n, dn, 1)

!     gradient length and steepest descent direction and length

glen  = dnrm2(n,g,1)
alpha = glen**2

call dcopy(n, g, 1, d, 1)
call dtrmv('U','N','N',n,rjac,ldr,d,1)
beta = dnrm2(n,d,1)**2

call dcopy(n, g, 1, ssd, 1)
call dscal(n, -(alpha/beta), ssd, 1)

ssdlen = alpha*glen/beta

!     set trust radius to ssdlen or dnlen if required

if( delta .eq. -Rone ) then
delta = min(ssdlen, stepmx)
elseif( delta .eq. -Rtwo ) then
delta = min(dnlen, stepmx)
endif

retcd = 4
gcnt  = 0

do while( retcd .gt. 1 )
!        find new step by More Hebden LM  algorithm
!        reuse ssd as sdiag

call nwmhstep(Rjac,ldr,n,ssd,qtf,dn,dnlen,glen,delta,mu, &
& d, v, dtype)
nwtstep = dtype .eq. 2
!        compute the model prediction 0.5*||F + J*d||**2 (L2-norm)

call dcopy(n,d,1,wa,1)
call dtrmv('U','N','N',n,rjac,ldr,wa,1)
call daxpy(n, Rone, qtf,1,wa,1)
fpred = Rhalf * dnrm2(n,wa,1)**2

!        evaluate function at xp = xc + d

do i=1,n
xp(i) = xc(i) + d(i)
enddo

call nwfvec(xp,n,scalex,fvec,fp,fpnorm,wa)
gcnt = gcnt + 1

!        check whether the global step is acceptable

oarg(2) = delta
call nwtrup(n,fcnorm,g,d,nwtstep,stepmx,xtol,delta, &
& fpred,retcd,xprev,fpnsav,fprev,xp,fp,fpnorm)

if( priter .gt. 0 ) then
oarg(1) = mu
oarg(3) = delta
oarg(4) = dnrm2(n, d, 1)
oarg(5) = fpnorm
oarg(6) = abs(fp(idamax(n,fp,1)))
call nwmhot(iter,dtype,retcd,oarg)
endif

enddo

return
end

!-----------------------------------------------------------------------

subroutine nwmhstep(R,ldr,n,sdiag,qtf,dn,dnlen,glen,delta,mu, &
& d, work, dtype)
integer ldr, n
double precision R(ldr,*)
double precision sdiag(*), qtf(*), dn(*), d(*), work(*)
double precision  dnlen, glen, delta, mu
integer dtype

!-------------------------------------------------------------------------
!
!     Find a new step by the More Hebden Levemberg Marquardt algorithm
!     Internal routine for nwmhlm
!
!     Arguments
!
!     In       R       Real(ldr,*)     R of QR-factored jacobian
!     In       ldr     Integer         leading dimension of R
!     In       n       Integer         size of problem
!     Out      sdiag   Real(*)         diagonal of LM lower triangular modified R
!     In       qtf     Real(*)         trans(Q)*f(xc)
!     In       dn      Real(*)         current newton step
!     Out      dnlen   Real            length dn()
!     In       glen    Real            length gradient
!     In       delta   Real            current trust region radius
!     Inout    mu      Real            Levenberg-Marquardt parameter
!     Out      d       Real(*)         new step for x()
!     Work     work    Real(*)         work vector for limhpar
!     Out      dtype   Integer         steptype
!                                       1 LM step
!                                       2 full newton direction
!
!-----------------------------------------------------------------------

double precision Rone
parameter(Rone=1.0D0)

if(dnlen .le. delta) then

!        Newton step smaller than trust radius ==> take it

call dcopy(n, dn, 1, d, 1)
delta = dnlen
dtype = 2

else

!        calculate LM step
call limhpar(R, ldr, n, sdiag, qtf, dn, dnlen, glen, delta, &
& mu, d, work)
!        change sign of step d (limhpar solves for trans(R)*R+mu*I)=qtf instead of -qtf)
call dscal(n,-Rone,d,1)
dtype = 1
endif

return
end

! ============================================================================
! From converted legacy source: nwnjac.f90
! ============================================================================
subroutine nwnjac(rjac,ldr,n,xc,fc,fq,fvec,fjac,epsm,jacflg, &
& wrk1,wrk2,wrk3, &
& xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn, &
& qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)

!-----------------------------------------------------------------------
!
!     Compute Jacobian matrix in xc, fc
!     scale it, compute gradient in xc and generate QR decomposition
!     calculate Newton step
!
!     Arguments
!
!     Out      rjac    Real(ldr,*)     jacobian (n columns)
!     In       ldr     Integer         leading dimension of rjac
!     In       n       Integer         dimensions of problem
!     In       xc      Real(*)         initial estimate of solution
!     Inout    fc      Real(*)         function values f(xc)
!     Wk       fq      Real(*)         workspace
!     In       fjac    Name            name of routine to calculate jacobian
!                                      (optional)
!     In       fvec    Name            name of routine to calculate f()
!     In       epsm    Real            machine precision
!     In       jacflg  Integer(*)      jacobian flag array
!                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
!                                                  3: user supplied banded
!                                      jacflg[2]: number of sub diagonals or -1 if not banded
!                                      jacflg[3]: number of super diagonals or -1 if not banded
!                                      jacflg[4]: 1 if adjusting step allowed when
!                                                   singular or illconditioned
!     Wk       wrk1    Real(*)         workspace
!     Wk       wrk2    Real(*)         workspace
!     Wk       wrk3    Real(*)         workspace
!     In       xscalm  Integer         x scaling method
!                                        1 from column norms of first jacobian
!                                          increased if needed after first iteration
!                                        0 scaling user supplied
!     Inout    scalex  Real(*)         scaling factors x(*)
!     Out      gp      Real(*)         gradient at xp()
!     In       cndtol  Real            tolerance of test for ill conditioning
!     Wk       rcdwrk  Real(*)         workspace
!     Wk       icdwrk  Integer(*)      workspace
!     Out      dn      Real(*)         Newton step
!     Out      qtf     Real(*)         workspace for nwnstp
!     Out      rcond   Real            estimated inverse condition of R from QR
!     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
!     In       qrwsiz  Integer         size of qrwork
!     Out      njcnt   Integer         number of jacobian evaluations
!     In       iter    Integer         iteration counter (used in scaling)
!     Inout    fstjac  logical         .true. if initial jacobian is available
!                                      on exit set to .false.
!     Out      ierr    Integer         error code
!                                        0 no error
!                                       >0 error in nwnstp (singular ...)
!
!-----------------------------------------------------------------------

integer ldr,n,iter, njcnt, ierr
integer jacflg(*),xscalm,qrwsiz
logical fstjac
double precision  epsm, cndtol, rcond
double precision  rjac(ldr,*)
double precision  xc(*),fc(*),dn(*)
double precision  wrk1(*),wrk2(*),wrk3(*)
double precision  qtf(*),gp(*),fq(*)
double precision  scalex(*)
double precision  rcdwrk(*),qrwork(*)
integer           icdwrk(*)
procedure(legacy_fjac) :: fjac
procedure(legacy_fvec) :: fvec

logical stepadj
double precision Rzero, Rone
parameter(Rzero=0.0d0, Rone=1.0d0)

!     evaluate the jacobian at the current iterate xc

if( .not. fstjac ) then
call nwfjac(xc,scalex,fc,fq,n,epsm,jacflg,fvec,fjac,rjac, &
& ldr,wrk1,wrk2,wrk3)
njcnt = njcnt + 1
else
fstjac = .false.
endif

!     if requested calculate x scale from jacobian column norms a la Minpack

if( xscalm .eq. 1 ) then
call vunsc(n,xc,scalex)
call nwcpsx(n,rjac,ldr,scalex,epsm,iter)
call vscal(n,xc,scalex)
endif

call nwscjac(n,rjac,ldr,scalex)

!     evaluate the gradient at the current iterate xc
!     gp = trans(Rjac) * fc
call dgemv('T',n,n,Rone,rjac,ldr,fc,1,Rzero,gp,1)

!     get newton step
stepadj = jacflg(4) .eq. 1
call dcopy(n,fc,1,fq,1)
call nwnstp(rjac,ldr,fq,n,cndtol, stepadj, &
& wrk1,dn,qtf,ierr,rcond, &
& rcdwrk,icdwrk,qrwork,qrwsiz)

!     save some data about jacobian for later output
call nwsnot(0,ierr,rcond)

return
end

!-----------------------------------------------------------------------

subroutine nwnstp(rjac,ldr,fn,n,cndtol, stepadj, &
& qraux,dn,qtf,ierr,rcond, &
& rcdwrk,icdwrk,qrwork,qrwsiz)

integer ldr,n,ierr,qrwsiz
double precision  cndtol,rjac(ldr,*),qraux(*),fn(*)
double precision  dn(*),qtf(*)
double precision  rcdwrk(*),qrwork(*)
integer           icdwrk(*)
double precision  rcond
logical           stepadj

!-----------------------------------------------------------------------
!
!     Calculate the newton step
!
!     Arguments
!
!     Inout    rjac    Real(ldr,*)     jacobian matrix at current iterate
!                                      overwritten with QR decomposition
!     In       ldr     Integer         leading dimension of rjac
!     In       fn      Real(*)         function values at current iterate
!     In       n       Integer         dimension of problem
!     In       cndtol  Real            tolerance of test for ill conditioning
!     In       stepadj Logical         allow adjusting step for singular/illconditioned jacobian
!     Inout    qraux   Real(*)         QR info from liqrfa (calling Lapack dgeqrf)
!     Out      dn      Real(*)         Newton direction
!     Out      qtf     Real(*)         trans(Q)*f()
!     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
!                                      1 indicating Jacobian ill-conditioned
!                                      2 indicating Jacobian completely singular
!                                      3 indicating almost zero LM correction
!     Out      rcond   Real            inverse condition of upper triangular R of QR
!     Wk       rcdwrk  Real(*)         workspace
!     Wk       icdwrk  Integer(*)      workspace
!     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
!     In       qrwsiz  Integer         size of qrwork
!
!-----------------------------------------------------------------------

integer info

double precision Rone
parameter(Rone=1.0d0)

!     perform a QR factorization of rjac (simple Lapack routine)
!     check for singularity or ill conditioning
!     form qtf = trans(Q) * fn

call liqrfa(Rjac,ldr,n,qraux,qrwork,qrwsiz,ierr)

!     compute qtf = trans(Q)*fn

call dcopy(n,fn,1,qtf,1)
call liqrqt(Rjac, ldr, n, qraux, qtf, qrwork, qrwsiz, info)

call lirslv(Rjac,ldr,n,cndtol, stepadj, &
& qtf,dn,ierr,rcond, rcdwrk,icdwrk)

return
end

! ============================================================================
! From converted legacy source: nwnleq.f90
! ============================================================================

subroutine nwnleq(x0,n,scalex,maxit, &
& jacflg,xtol,ftol,btol,cndtol,method,global, &
& xscalm,stepmx,delta,sigma,rjac,ldr, &
& rwork,lrwork, &
& rcdwrk,icdwrk,qrwork,qrwsiz,fjac,fvec,outopt,xp, &
& fp,gp,njcnt,nfcnt,iter,termcd)

integer n,jacflg(*),maxit,njcnt,nfcnt,iter,termcd,method
integer global,xscalm,ldr,lrwork,qrwsiz
integer outopt(*)
double precision  xtol,ftol,btol,cndtol,stepmx,delta,sigma
double precision  xp(*),fp(*),gp(*),x0(*)
double precision  rjac(ldr,*),rwork(*),rcdwrk(*),qrwork(*)
double precision  scalex(*)
integer           icdwrk(*)
procedure(legacy_fjac) :: fjac
procedure(legacy_fvec) :: fvec

!-------------------------------------------------------------------------
!
!     Solves systems of nonlinear equations using the Newton / Broyden
!     method with a global strategy either linesearch or double dogleg
!
!     In       x0      Real(*)         starting vector for x
!     In       n       Integer         dimension of problem
!     Inout    scalex  Real(*)         scaling factors x()
!     Inout    maxit   Integer         maximum number iterations
!     Inout    jacflg  Integer(*)      jacobian flag array
!                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
!                                                  3: user supplied banded
!                                      jacflg[2]: number of sub diagonals or -1 if not banded
!                                      jacflg[3]: number of super diagonals or -1 if not banded
!                                      jacflg[4]: 1 if adjusting step allowed when
!                                                   singular or illconditioned
!     Inout    xtol    Real            x tolerance
!     Inout    ftol    Real            f tolerance
!     Inout    btol    Real            x tolerance for backtracking
!     Inout    cndtol  Real            tolerance of test for ill conditioning
!     Inout    method  Integer         method to use
!                                        0 Newton
!                                        1 Broyden
!     In       global  Integer         global strategy to use
!                                        1 cubic linesearch
!                                        2 quadratic linesearch
!                                        3 geometric linesearch
!                                        4 double dogleg
!                                        5 powell dogleg
!                                        6 hookstep (More-Hebden Levenberg-Marquardt)
!     In       xscalm  Integer         scaling method
!                                        0 scale fixed and supplied by user
!                                        1 for scale from jac. columns a la Minpack
!     Inout    stepmx  Real            maximum stepsize
!     Inout    delta   Real            trust region radius
!                                        > 0.0 or special value for initial value
!                                        -1.0  ==> use min(Cauchy length, stepmx)
!                                        -2.0  ==> use min(Newton length, stepmx)
!     Inout    sigma   Real            reduction factor geometric linesearch
!     Inout    rjac    Real(ldr,*)     workspace jacobian
!                                         2*n*n for Broyden and n*n for Newton
!     In       ldr     Integer         leading dimension rjac
!     Out      rwork   Real(*)         real workspace (9n)
!     In       lrwork  Integer         size real workspace
!     In       rcdwrk  Real(*)         workspace for Dtrcon (3n)
!     In       icdwrk  Integer(*)      workspace for Dtrcon (n)
!     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
!     In       qrwsiz  Integer         size of qrwork
!     In       fjac    Name            optional name of routine to calculate
!                                      user supplied jacobian
!     In       fvec    Name            name of routine to calculate f(x)
!     In       outopt  Integer(*)      output options
!                                       outopt(1)
!                                         0 no output
!                                         1 output an iteration report
!                                       outopt(2)
!                                         0 do not check any user supplied jacobian
!                                         1 check user supplied jacobian if supplied
!     Out      xp      Real(*)         final values for x()
!     Out      fp      Real(*)         final values for f(x)
!     Out      gp      Real(*)         gradient of f() at xp()
!     Out      njcnt   Integer         number of jacobian evaluations
!     Out      nfcnt   Integer         number of function evaluations
!     Out      iter    Integer         number of (outer) iterations
!     Out      termcd  Integer         termination code
!                                       > 0 process terminated
!                                             1  function criterion near zero
!                                             2  no better point found
!                                             3  x-values within tolerance
!                                             4  iteration limit exceeded
!                                             5  singular/ill-conditioned jacobian
!                                             6  totally singular jacobian
!                                                (when allowSingular=TRUE)
!
!                                       < 0 invalid input parameters
!                                            -1  n not positive
!                                            -2  insufficient workspace rwork
!                                            -3  cannot check user supplied jacobian (not supplied)
!
!    The subroutine fvec must be declared as
!
!!        subroutine fvec(x,f,n,flag)
!         double precision x(*), f(*)
!         integer  n, flag
!
!         x() are the x values for which to calculate the function values f(*)
!         The dimension of these vectors is n
!         The flag argument is set to
!            0  for calculation of function values
!           >0  indicating that jacobian column <flag> is being computed
!               so that fvec can abort.
!
!    The subroutine fjac must be declared as
!
!!        subroutine mkjac(rjac,ldr,x,n)
!         integer ldr
!         double precision rjac(ldr,*), x(*)
!         integer  n
!
!         The routine calculates the jacobian in point x(*) of the
!         function. If any illegal values are encountered during
!         calculation of the jacobian it is the responsibility of
!         the routine to quit.

!-------------------------------------------------------------------------

double precision epsm

!     check input parameters

call nwpchk(n,lrwork,xtol,ftol,btol,cndtol,maxit, &
& jacflg,method,global,stepmx,delta,sigma, &
& epsm,outopt,scalex,xscalm,termcd)
if(termcd .lt. 0) then
return
endif

!     first argument of nwsolv/brsolv is leading dimension of rjac in those routines
!     should be at least n

if( method .eq. 0 ) then

call nwsolv(ldr,x0,n,scalex,maxit,jacflg, &
& xtol,ftol,btol,cndtol,global,xscalm, &
& stepmx,delta,sigma, &
& rjac, &
& rwork(1    ),rwork(1+  n), &
& rwork(1+2*n),rwork(1+3*n), &
& rwork(1+4*n),rwork(1+5*n), &
& rwork(1+6*n),rwork(1+7*n), &
& rwork(1+8*n),rcdwrk,icdwrk,qrwork,qrwsiz,epsm, &
& fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,iter,termcd)

elseif( method .eq. 1 ) then

call brsolv(ldr,x0,n,scalex,maxit,jacflg, &
& xtol,ftol,btol,cndtol,global,xscalm, &
& stepmx,delta,sigma, &
& rjac,rjac(1,n+1), &
& rwork(1    ),rwork(1+  n), &
& rwork(1+2*n),rwork(1+3*n), &
& rwork(1+4*n),rwork(1+5*n), &
& rwork(1+6*n),rwork(1+7*n), &
& rwork(1+8*n),rcdwrk,icdwrk, &
& qrwork,qrwsiz,epsm, &
& fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,iter,termcd)

endif

return
end

!-----------------------------------------------------------------------

subroutine nwpchk(n,lrwk, &
& xtol,ftol,btol,cndtol,maxit,jacflg,method, &
& global,stepmx,delta,sigma,epsm,outopt, &
& scalex,xscalm,termcd)

integer n,lrwk,jacflg(*)
integer method,global,maxit,xscalm,termcd
integer outopt(*)
double precision  xtol,ftol,btol,cndtol,stepmx,delta,sigma,epsm
double precision  scalex(*)

!-------------------------------------------------------------------------
!
!     Check input arguments for consistency and modify if needed/harmless
!
!     Arguments
!
!     In       n       Integer         dimension of problem
!     In       lrwk    Integer         size real workspace
!     Inout    xtol    Real            x tolerance
!     Inout    ftol    Real            f tolerance
!     Inout    btol    Real            x tolerance for backtracking
!     Inout    cndtol  Real            tolerance of test for ill conditioning
!     Inout    maxit   Integer         maximum number iterations
!     Inout    jacflg  Integer(*)      jacobian flag
!     Inout    method  Integer         method to use (Newton/Broyden)
!     Inout    global  Integer         global strategy to use
!     Inout    stepmx  Real            maximum stepsize
!     Inout    delta     Real            trust region radius
!     Inout    sigma   Real            reduction factor geometric linesearch
!     Out      epsm                    machine precision
!     Inout    scalex  Real(*)         scaling factors x()
!     Inout    xscalm  integer         0 for fixed scaling, 1 for automatic scaling
!     Out      termcd  Integer         termination code (<0 on errors)
!
!-------------------------------------------------------------------------

integer i,len
double precision Rhuge

double precision Rzero, Rone, Rtwo, Rthree
parameter(Rzero=0.0d0, Rone=1.0d0, Rtwo=2.0d0, Rthree=3.0d0)

double precision Rhalf
parameter(Rhalf = 0.5d0)

!     check that parameters only take on acceptable values
!     if not, set them to default values

!     initialize termcd to all ok

termcd = 0

!     compute machine precision

epsm = epsmch()

!     get largest double precision number
Rhuge = dblhuge()

!     check dimensions of the problem

if(n .le. 0) then
termcd = -1
return
endif

!     check dimensions of workspace arrays

len = 9*n
!      +2*n*n
if(lrwk .lt. len) then
termcd = -2
return
endif

!     check jacflg, method, and global

if(jacflg(1) .gt. 3 .or. jacflg(1) .lt. 0) jacflg(1) = 0

if(method .lt. 0 .or. method .gt. 1) method = 0

if(global .lt. 0 .or. global .gt. 6) global = 4

!     set outopt to correct values

if(outopt(1) .ne. 0 ) then
outopt(1) = 1
endif

if(outopt(2) .ne. 0 ) then
outopt(2) = 1
endif

!     check scaling scale matrices

if(xscalm .eq. 0) then
do i = 1,n
if(scalex(i) .lt. Rzero) scalex(i) = -scalex(i)
if(scalex(i) .eq. Rzero) scalex(i) = Rone
enddo
else
xscalm = 1
do i = 1,n
scalex(i) = Rone
enddo
endif
!     check step and function tolerances

if(xtol .lt. Rzero) then
xtol = epsm**(Rtwo/Rthree)
endif

if(ftol .lt. Rzero) then
ftol = epsm**(Rtwo/Rthree)
endif

if( btol .lt. xtol ) btol = xtol

cndtol = max(cndtol, epsm)

!     check reduction in geometric linesearch

if( sigma .le. Rzero .or. sigma .ge. Rone ) then
sigma = Rhalf
endif

!     check iteration limit

if(maxit .le. 0) then
maxit = 150
endif

!     set stepmx

if(stepmx .le. Rzero) stepmx = Rhuge

!     check delta
if(delta .le. Rzero) then
if( delta .ne. -Rtwo ) then
delta = -Rone
endif
elseif(delta .gt. stepmx) then
delta = stepmx
endif

return
end

! ============================================================================
! From converted legacy source: nwnwtn.f90
! ============================================================================

subroutine nwsolv(ldr,xc,n,scalex,maxit, &
& jacflg,xtol,ftol,btol,cndtol,global,xscalm, &
& stepmx,delta,sigma, &
& rjac,wrk1,wrk2,wrk3,wrk4,fc,fq,dn,d,qtf, &
& rcdwrk,icdwrk,qrwork,qrwsiz,epsm, &
& fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,iter, &
& termcd)

integer ldr,n,termcd,njcnt,nfcnt,iter
integer maxit,jacflg(*),global,xscalm,qrwsiz
integer outopt(*)
double precision  xtol,ftol,btol,cndtol
double precision  stepmx,delta,sigma,fpnorm,epsm
double precision  rjac(ldr,*)
double precision  xc(*),fc(*),xp(*),fp(*),dn(*),d(*)
double precision  wrk1(*),wrk2(*),wrk3(*),wrk4(*)
double precision  qtf(*),gp(*),fq(*)
double precision  scalex(*)
double precision  rcdwrk(*),qrwork(*)
integer           icdwrk(*)
procedure(legacy_fjac) :: fjac
procedure(legacy_fvec) :: fvec

!-----------------------------------------------------------------------
!
!     Solve system of nonlinear equations with Newton and global strategy
!
!
!     Arguments
!
!     In       ldr     Integer         leading dimension of rjac
!     In       xc      Real(*)         initial estimate of solution
!     In       n       Integer         dimensions of problem
!     Inout    scalex  Real(*)         scaling factors x(*)
!     In       maxit   Integer         maximum number of allowable iterations
!     In       jacflg  Integer(*)      jacobian flag array
!                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
!                                                  3: user supplied banded
!                                      jacflg[2]: number of sub diagonals or -1 if not banded
!                                      jacflg[3]: number of super diagonals or -1 if not banded
!                                      jacflg[4]: 1 if adjusting step allowed when
!                                                   singular or illconditioned
!     In       xtol    Real            tolerance at which successive iterates x()
!                                      are considered close enough to
!                                      terminate algorithm
!     In       ftol    Real            tolerance at which function values f()
!                                      are considered close enough to zero
!     Inout    btol    Real            x tolerance for backtracking
!     Inout    cndtol  Real            tolerance of test for ill conditioning
!     In       global  Integer         global strategy to use
!                                        1 cubic linesearch
!                                        2 quadratic linesearch
!                                        3 geometric linesearch
!                                        4 double dogleg
!                                        5 powell dogleg
!                                        6 hookstep (More-Hebden Levenberg-Marquardt)
!     In       xscalm  Integer         x scaling method
!                                        1 from column norms of first jacobian
!                                          increased if needed after first iteration
!                                        0 scaling user supplied
!     In       stepmx  Real            maximum allowable step size
!     In       delta   Real            trust region radius
!     In       sigma   Real            reduction factor geometric linesearch
!     Inout    rjac    Real(ldr,*)     jacobian (n columns)
!     Wk       wrk1    Real(*)         workspace
!     Wk       wrk2    Real(*)         workspace
!     Wk       wrk3    Real(*)         workspace
!     Wk       wrk4    Real(*)         workspace
!     Inout    fc      Real(*)         function values f(xc)
!     Wk       fq      Real(*)         workspace
!     Wk       dn      Real(*)         workspace
!     Wk       d       Real(*)         workspace
!     Wk       qtf     Real(*)         workspace
!     Wk       rcdwrk  Real(*)         workspace
!     Wk       icdwrk  Integer(*)      workspace
!     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
!     In       qrwsiz  Integer         size of qrwork
!     In       epsm    Real            machine precision
!     In       fjac    Name            name of routine to calculate jacobian
!                                      (optional)
!     In       fvec    Name            name of routine to calculate f()
!     In       outopt  Integer(*)      output options
!     Out      xp      Real(*)         final x()
!     Out      fp      Real(*)         final f(xp)
!     Out      gp      Real(*)         gradient at xp()
!     Out      njcnt   Integer         number of jacobian evaluations
!     Out      nfcnt   Integer         number of function evaluations
!     Out      iter    Integer         number of (outer) iterations
!     Out      termcd  Integer         termination code
!
!-----------------------------------------------------------------------

integer gcnt,retcd,ierr
double precision  dum(2),fcnorm,rcond
logical fstjac
integer priter



!     initialization

retcd = 0
iter  = 0
njcnt = 0
nfcnt = 0
ierr  = 0

dum(1) = 0

if( outopt(1) .eq. 1 ) then
priter = 1
else
priter = -1
endif

!     evaluate function

call vscal(n,xc,scalex)
call nwfvec(xc,n,scalex,fvec,fc,fcnorm,wrk1)

!     evaluate user supplied or finite difference jacobian and check user supplied
!     jacobian, if requested

fstjac = .false.
if(mod(jacflg(1),2) .eq. 1) then

if( outopt(2) .eq. 1 ) then
fstjac = .true.
njcnt = njcnt + 1
call nwfjac(xc,scalex,fc,fq,n,epsm,jacflg,fvec,fjac,rjac, &
& ldr,wrk1,wrk2,wrk3)
call chkjac(rjac,ldr,xc,fc,n,epsm,jacflg,scalex, &
& fq,wrk1,wrk2,fvec,termcd)
if(termcd .lt. 0) then
!              copy initial values
call dcopy(n,xc,1,xp,1)
call dcopy(n,fc,1,fp,1)
call vunsc(n,xp,scalex)
fpnorm = fcnorm
return
endif
endif

endif

!     check stopping criteria for input xc

call nwtcvg(xc,fc,xc,xtol,retcd,ftol,iter,maxit,n,ierr,termcd)

if(termcd .gt. 0) then
call dcopy(n,xc,1,xp,1)
call dcopy(n,fc,1,fp,1)
fpnorm = fcnorm
if( outopt(3) .eq. 1 .and. .not. fstjac ) then
njcnt = njcnt + 1
call nwfjac(xp,scalex,fp,fq,n,epsm,jacflg,fvec,fjac,rjac, &
& ldr,wrk1,wrk2,wrk3)
endif
return
endif

if( priter .gt. 0 ) then

dum(1) = fcnorm
dum(2) = abs(fc(idamax(n,fc,1)))

if( global .eq. 0 ) then
call nwprot(iter, -1, dum)
elseif( global .le. 3 ) then
call nwlsot(iter,-1,dum)
elseif( global .eq. 4 ) then
call nwdgot(iter,-1,0,dum)
elseif( global .eq. 5 ) then
call nwpwot(iter,-1,0,dum)
elseif( global .eq. 6 ) then
call nwmhot(iter,-1,0,dum)
endif

endif

do while( termcd .eq. 0 )
iter = iter + 1

call nwnjac(rjac,ldr,n,xc,fc,fq,fvec,fjac,epsm,jacflg,wrk1, &
& wrk2,wrk3, &
& xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn, &
& qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)
!        - choose the next iterate xp by a global strategy

if( ierr .gt. 0 ) then
!           jacobian singular or too ill-conditioned
call nweset(n,xc,fc,fcnorm,xp,fp,fpnorm,gcnt,priter,iter)
elseif(global .eq. 0) then
call nwpure(n,xc,dn,stepmx,scalex, &
& fvec,xp,fp,fpnorm,wrk1,retcd,gcnt, &
& priter,iter)
elseif(global .eq. 1) then
call nwclsh(n,xc,fcnorm,dn,gp,stepmx,btol,scalex, &
& fvec,xp,fp,fpnorm,wrk1,retcd,gcnt, &
& priter,iter)
elseif(global .eq. 2) then
call nwqlsh(n,xc,fcnorm,dn,gp,stepmx,btol,scalex, &
& fvec,xp,fp,fpnorm,wrk1,retcd,gcnt, &
& priter,iter)
elseif(global .eq. 3) then
call nwglsh(n,xc,fcnorm,dn,gp,sigma,stepmx,btol,scalex, &
& fvec,xp,fp,fpnorm,wrk1,retcd,gcnt, &
& priter,iter)
elseif(global .eq. 4) then
call nwddlg(n,rjac,ldr,dn,gp,xc,fcnorm,stepmx, &
& btol,delta,qtf,scalex, &
& fvec,d,fq,wrk1,wrk2,wrk3,wrk4, &
& xp,fp,fpnorm,retcd,gcnt,priter,iter)
elseif(global .eq. 5) then
call nwpdlg(n,rjac,ldr,dn,gp,xc,fcnorm,stepmx, &
& btol,delta,qtf,scalex, &
& fvec,d,fq,wrk1,wrk2,wrk3,wrk4, &
& xp,fp,fpnorm,retcd,gcnt,priter,iter)
elseif(global .eq. 6) then
call nwmhlm(n,rjac,ldr,dn,gp,xc,fcnorm,stepmx, &
& btol,delta,qtf,scalex, &
& fvec,d,fq,wrk1,wrk2,wrk3,wrk4, &
& xp,fp,fpnorm,retcd,gcnt,priter,iter)
endif

nfcnt = nfcnt + gcnt

!        - check stopping criteria for the new iterate xp

call nwtcvg(xp,fp,xc,xtol,retcd,ftol,iter,maxit,n,ierr,termcd)

if(termcd .eq. 0) then
!           update xc, fc, and fcnorm
call dcopy(n,xp,1,xc,1)
call dcopy(n,fp,1,fc,1)
fcnorm = fpnorm
endif

enddo

if( outopt(3) .eq. 1 ) then
call nwfjac(xp,scalex,fp,fq,n,epsm,jacflg,fvec,fjac,rjac, &
& ldr,wrk1,wrk2,wrk3)
endif

call vunsc(n,xp,scalex)

return
end

! ============================================================================
! From converted legacy source: nwpdlg.f90
! ============================================================================

subroutine nwpdlg(n,rjac,ldr,dn,g,xc,fcnorm,stepmx,xtol, &
& delta,qtf,scalex,fvec,d,xprev, &
& ssd,v,wa,fprev,xp,fp,fpnorm,retcd,gcnt, &
& priter,iter)

integer ldr, n, retcd, gcnt, priter, iter
double precision  fcnorm, stepmx, xtol, fpnorm, delta
double precision  rjac(ldr,*), dn(*), g(*), xc(*), qtf(*)
double precision  scalex(*), d(*)
double precision  xprev(*), xp(*), fp(*)
double precision  ssd(*), v(*), wa(*), fprev(*)
procedure(legacy_fvec) :: fvec

!-------------------------------------------------------------------------
!
!     Find a next iterate xp by the Powell dogleg method
!
!     Arguments
!
!     In       n       Integer         size of problem: dimension x, f
!     In       Rjac    Real(ldr,*)     R of QR-factored jacobian
!     In       ldr     Integer         leading dimension of Rjac
!     Inout    dn      Real(*)         newton direction
!     Inout    g       Real(*)         gradient at current point
!                                      trans(jac)*f()
!     In       xc      Real(*)         current iterate
!     In       fcnorm  Real            .5*||f(xc)||**2
!     In       stepmx  Real            maximum stepsize
!     In       xtol    Real            x-tolerance (stepsize)
!     Inout    delta     Real            on input: initial trust region radius
!                                                if -1 then set to something
!                                                reasonable
!                                      on output: final value
!                                      ! Do not modify between calls while
!                                        still iterating
!     In       qtf     Real(*)         trans(Q)*f(xc)
!     In       scalex  Real(*)         scaling factors for x()
!     In       fvec    Name            name of subroutine to evaluate f(x)
!                                      ! must be declared external in caller
!     Wk       d       Real(*)         work vector
!     Wk       xprev   Real(*)         work vector
!     Wk       ssd     Real(*)         work vector
!     Wk       v       Real(*)         work vector
!     Wk       wa      Real(*)         work vector
!     Wk       fprev   Real(*)         work vector
!     Out      xp      Real(*)         new x()
!     Out      fp      Real(*)         new f(xp)
!     Out      fpnorm  Real            new .5*||f(xp)||**2
!     Out      retcd   Integer         return code
!                                       0  new satisfactory x() found
!                                       1  no  satisfactory x() found
!     Out      gcnt    Integer         number of steps taken
!     In       priter  Integer         print flag
!                                       -1 no intermediate printing
!                                       >0 yes for print of intermediate results
!     In       iter    Integer         current iteration (only used for above)
!
!     All vectors at least size n
!
!-------------------------------------------------------------------------

integer i
double precision  dnlen,ssdlen,alpha,beta,lambda,fpred
double precision  sqalpha,fpnsav,oarg(5)

logical nwtstep
integer dtype



double precision Rone, Rtwo, Rhalf
parameter(Rhalf=0.5d0)
parameter(Rone=1.0d0, Rtwo=2.0d0)

!     length newton direction

dnlen = dnrm2(n, dn, 1)

!     steepest descent direction and length

sqalpha = dnrm2(n,g,1)
alpha   = sqalpha**2

call dcopy(n, g, 1, d, 1)
call dtrmv('U','N','N',n,rjac,ldr,d,1)
beta = dnrm2(n,d,1)**2

call dcopy(n, g, 1, ssd, 1)
call dscal(n, -(alpha/beta), ssd, 1)

ssdlen = alpha*sqalpha/beta

!     set trust radius to ssdlen or dnlen if required

if( delta .eq. -Rone ) then
delta = min(ssdlen, stepmx)
elseif( delta .eq. -Rtwo ) then
delta = min(dnlen, stepmx)
endif

retcd = 4
gcnt  = 0

do while( retcd .gt. 1 )

!        find new step by single dogleg algorithm

call pwlstp(n,dn,dnlen,delta,v, &
& ssd,ssdlen,d,dtype,lambda)
nwtstep = dtype .eq. 3
!        compute the model prediction 0.5*||F + J*d||**2 (L2-norm)

call dcopy(n,d,1,wa,1)
call dtrmv('U','N','N',n,rjac,ldr,wa,1)
call daxpy(n, Rone, qtf,1,wa,1)
fpred = Rhalf * dnrm2(n,wa,1)**2

!        evaluate function at xp = xc + d

do i=1,n
xp(i) = xc(i) + d(i)
enddo

call nwfvec(xp,n,scalex,fvec,fp,fpnorm,wa)
gcnt = gcnt + 1

!        check whether the global step is acceptable

oarg(2) = delta
call nwtrup(n,fcnorm,g,d,nwtstep,stepmx,xtol,delta, &
& fpred,retcd,xprev,fpnsav,fprev,xp,fp,fpnorm)

if( priter .gt. 0 ) then
oarg(1) = lambda
oarg(3) = delta
oarg(4) = fpnorm
oarg(5) = abs(fp(idamax(n,fp,1)))
call nwpwot(iter,dtype,retcd,oarg)
endif

enddo

return
end

!-----------------------------------------------------------------------

subroutine pwlstp(n,dn,dnlen,delta,v, &
& ssd,ssdlen,d,dtype,lambda)
integer n
double precision  dn(*), ssd(*), v(*), d(*)
double precision  dnlen, delta, ssdlen, lambda
integer dtype

!-------------------------------------------------------------------------
!
!     Find a new step by the Powell dogleg algorithm
!     Internal routine for nwpdlg
!
!     Arguments
!
!     In       n       Integer         size of problem
!     In       dn      Real(*)         current newton step
!     Out      dnlen   Real            length dn()
!     In       delta   Real            current trust region radius
!     Out      v       Real(*)         (internal) dn() - ssd()
!     In       ssd     Real(*)         (internal) steepest descent direction
!     In       ssdlen  Real            (internal) length ssd
!     Out      d       Real(*)         new step for x()
!     Out      dtype   Integer         steptype
!                                       1 steepest descent
!                                       2 combination of dn and ssd
!                                       3 full newton direction
!     Out      lambda  Real            weight of dn() in d()
!                                      closer to 1 ==> more of dn()
!
!-----------------------------------------------------------------------

integer i
double precision vssd, vlen


if(dnlen .le. delta) then

!        Newton step smaller than trust radius ==> take it

call dcopy(n, dn, 1, d, 1)
delta = dnlen
dtype = 3

elseif(ssdlen .ge. delta) then

!        take step in steepest descent direction

call dcopy(n, ssd, 1, d, 1)
call dscal(n, delta / ssdlen, d, 1)
dtype = 1

else

!        calculate convex combination of ssd and dn with length delta

do i=1,n
v(i) = dn(i) - ssd(i)
enddo

vssd = ddot(n,v,1,ssd,1)
vlen = dnrm2(n,v,1)**2

lambda =(-vssd+sqrt(vssd**2-vlen*(ssdlen**2-delta**2)))/vlen
call dcopy(n, ssd, 1, d, 1)
call daxpy(n, lambda, v, 1, d, 1)
dtype = 2

endif

return
end

! ============================================================================
! From converted legacy source: nwpure.f90
! ============================================================================

subroutine nwpure(n,xc,d,stepmx,scalex,fvec, &
& xp,fp,fpnorm,xw,retcd,gcnt,priter,iter)

integer n,retcd,gcnt
double precision  stepmx,fpnorm
double precision  xc(*)
double precision  d(*),xp(*),fp(*),xw(*)
double precision  scalex(*)
procedure(legacy_fvec) :: fvec

integer priter,iter

!-------------------------------------------------------------------------
!
!     Find a next iterate using geometric line search
!     along the newton direction
!
!     Arguments
!
!     In       n       Integer          dimension of problem
!     In       xc      Real(*)          current iterate
!     In       d       Real(*)          newton direction
!     In       stepmx  Real             maximum stepsize
!     In       scalex  Real(*)          diagonal scaling matrix for x()
!     In       fvec    Name             name of routine to calculate f()
!     In       xp      Real(*)          new x()
!     In       fp      Real(*)          new f(x)
!     In       fpnorm  Real             .5*||fp||**2
!     Out      xw      Real(*)          workspace for unscaling x
!
!     Out      retcd   Integer          return code
!                                         0 new satisfactory x() found (!always)
!
!     Out      gcnt    Integer          number of steps taken
!     In       priter  Integer           >0 if intermediate steps to be printed
!                                        -1 if no printing
!
!-------------------------------------------------------------------------

integer i
double precision  oarg(3)
double precision  lambda

double precision  dlen



double precision Rone
parameter(Rone=1.0d0)

!     safeguard initial step size

dlen = dnrm2(n,d,1)
if( dlen .gt. stepmx ) then
lambda = stepmx / dlen
else
lambda = Rone
endif

retcd  = 0
gcnt   = 1

!     compute the next iterate xp

do i=1,n
xp(i) = xc(i) + lambda*d(i)
enddo

!     evaluate functions and the objective function at xp

call nwfvec(xp,n,scalex,fvec,fp,fpnorm,xw)

if( priter .gt. 0) then
oarg(1) = lambda
oarg(2) = fpnorm
oarg(3) = abs(fp(idamax(n,fp,1)))
call nwprot(iter,1,oarg)
endif

return
end

! ============================================================================
! From converted legacy source: nwqlsh.f90
! ============================================================================

subroutine nwqlsh(n,xc,fcnorm,d,g,stepmx,xtol,scalex,fvec, &
& xp,fp,fpnorm,xw,retcd,gcnt,priter,iter)

integer n,retcd,gcnt
double precision  stepmx,xtol,fcnorm,fpnorm
double precision  xc(*)
double precision  d(*),g(*),xp(*),fp(*),xw(*)
double precision  scalex(*)
procedure(legacy_fvec) :: fvec

integer priter,iter

!-------------------------------------------------------------------------
!
!     Find a next acceptable iterate using a safeguarded quadratic line search
!     along the newton direction
!
!     Arguments
!
!     In       n       Integer          dimension of problem
!     In       xc      Real(*)          current iterate
!     In       fcnorm  Real             0.5 * || f(xc) ||**2
!     In       d       Real(*)          newton direction
!     In       g       Real(*)          gradient at current iterate
!     In       stepmx  Real             maximum stepsize
!     In       xtol    Real             relative step size at which
!                                       successive iterates are considered
!                                       close enough to terminate algorithm
!     In       scalex  Real(*)          diagonal scaling matrix for x()
!     In       fvec    Name             name of routine to calculate f()
!     In       xp      Real(*)          new x()
!     In       fp      Real(*)          new f(x)
!     In       fpnorm  Real             .5*||fp||**2
!     Out      xw      Real(*)           workspace for unscaling x(*)
!
!     Out      retcd   Integer          return code
!                                         0 new satisfactory x() found
!                                         1 no  satisfactory x() found
!                                           sufficiently distinct from xc()
!
!     Out      gcnt    Integer          number of steps taken
!     In       priter  Integer           >0 unit if intermediate steps to be printed
!                                        -1 if no printing
!
!-------------------------------------------------------------------------

integer i
double precision  alpha,slope,rsclen,oarg(4)
double precision  lambda,lamhi,lamlo,t
double precision ftarg
double precision  dlen



parameter (alpha = 1.0d-4)

double precision Rone, Rtwo, Rten
parameter(Rone=1.0d0, Rtwo=2.0d0, Rten=10.0d0)

!     safeguard initial step size

dlen = dnrm2(n,d,1)
if( dlen .gt. stepmx ) then
lamhi  = stepmx / dlen
else
lamhi  = Rone
endif

!     compute slope  =  g-trans * d

slope = ddot(n,g,1,d,1)

!     compute the smallest value allowable for the damping
!     parameter lambda ==> lamlo

rsclen = nudnrm(n,d,xc)
lamlo  = xtol / rsclen

!     initialization of retcd and lambda (linesearch length)

retcd  = 2
lambda = lamhi
gcnt   = 0

do while( retcd .eq. 2 )

!        compute next x

do i=1,n
xp(i) = xc(i) + lambda*d(i)
enddo

!        evaluate functions and the objective function at xp

call nwfvec(xp,n,scalex,fvec,fp,fpnorm,xw)
gcnt = gcnt + 1
ftarg = fcnorm + alpha * lambda * slope

if( priter .gt. 0) then
oarg(1) = lambda
oarg(2) = ftarg
oarg(3) = fpnorm
oarg(4) = abs(fp(idamax(n,fp,1)))
call nwlsot(iter,1,oarg)
endif

!        test whether the standard step produces enough decrease
!        of the objective function.
!        If not update lambda and compute a new next iterate

if( fpnorm .le. ftarg ) then
retcd = 0
else
t = ((-lambda**2)*slope/Rtwo)/(fpnorm-fcnorm-lambda*slope)
lambda  = max(lambda / Rten , t)
if(lambda .lt. lamlo) then
retcd = 1
endif
endif

enddo

return
end

! ============================================================================
! From converted legacy source: nwtrup.f90
! ============================================================================

subroutine nwtrup(n,fcnorm,g,sc,nwtstep,stepmx,xtol,delta, &
& fpred,retcd,xprev,fpnsav,fprev,xp,fp, &
& fpnorm)

integer n,retcd
double precision  fcnorm,stepmx,xtol,delta,fpred,fpnsav,fpnorm
double precision  xp(*),g(*)
double precision  sc(*),xprev(*),fprev(*),fp(*)
logical nwtstep

!-------------------------------------------------------------------------
!
!     Decide whether to accept xp=xc+sc as the next iterate
!     and updates the trust region delta
!
!     Arguments
!
!     In       n       Integer         size of xc()
!     In       fcnorm  Real            .5*||f(xc)||**2
!     In       g       Real(*)         gradient at xc()
!     In       sc      Real(*)         current step
!     In       nwtstep Logical         true if sc is newton direction
!     In       stepmx  Real            maximum step size
!     In       xtol    Real            minimum step tolerance
!     Inout    delta   Real            trust region radius
!     In       fpred   Real            predicted value of .5*||f()||**2
!
!     Inout    retcd   Integer         return code
!                                       0 xp accepted as next iterate;
!                                         delta trust region for next iteration.
!
!                                       1 xp unsatisfactory but
!                                         accepted as next iterate because
!                                         xp-xc .lt. smallest allowable
!                                         step length.
!
!                                       2 f(xp) too large.
!                                         continue current iteration with
!                                         new reduced delta.
!
!                                       3 f(xp) sufficiently small, but
!                                         quadratic model predicts f(xp)
!                                         sufficiently well to continue current
!                                         iteration with new doubled delta.
!
!                                      On first entry, retcd must be 4
!
!     Wk       xprev   Real(*)         (internal) work
!     Wk       fpnsav  Real            (internal) work
!     Wk       fprev   Real(*)         (internal) work
!     Inout    xp      Real(*)         new iterate x()
!     Inout    fp      Real(*)         new f(xp)
!     Inout    fpnorm  Real            new .5*||f(xp)||**2
!
!-------------------------------------------------------------------------

double precision  ared,pred,slope,sclen,rln,dltmp

logical ret3ok

double precision Rone, Rtwo, Rthree, Rfour, Rten
double precision Rhalf, Rpten
parameter(Rpten = 0.1d0)
parameter(Rhalf=0.5d0)
parameter(Rone=1.0d0, Rtwo=2.0d0, Rthree=3.0d0, Rfour=4.0d0)
parameter(Rten=10.0d0)

double precision Rp99,Rp4, Rp75
parameter(Rp99=Rone-Rten**(-2), Rp4=Rten**(-4), Rp75=Rthree/Rfour)

double precision alpha
parameter(alpha = Rp4)

!     pred measures the predicted reduction in the function value

ared  = fpnorm - fcnorm
pred  = fpred  - fcnorm
slope = ddot(n,g,1,sc,1)

if(retcd .ne. 3) then
ret3ok = .false.
else
ret3ok = fpnorm .ge. fpnsav .or. ared .gt. alpha * slope
endif

if(retcd .eq. 3 .and. ret3ok) then

!        reset xp,fp,fpnorm to saved values and terminate global step

retcd = 0
call dcopy(n,xprev,1,xp,1)
call dcopy(n,fprev,1,fp,1)
fpnorm = fpnsav
!        reset delta to initial value
!        but beware
!           if the trial step was a Newton step then delta is reset to
!           .5 * length(Newton step) which will be smaller
!           because delta is set to length(Newton step) elsewhere
!           see nwddlg.f and nwpdlg.f
delta  = Rhalf*delta

elseif(ared .gt. alpha * slope) then

!        fpnorm too large (decrease not sufficient)

rln = nudnrm(n,sc,xp)
if(rln .lt. xtol) then

!           cannot find satisfactory xp sufficiently distinct from xc

retcd = 1

else

!           reduce trust region and continue current global step

retcd = 2
sclen = dnrm2(n,sc,1)
dltmp = -slope*sclen/(Rtwo*(ared-slope))

if(dltmp .lt. Rpten*delta) then
delta = Rpten*delta
else
delta = min(Rhalf*delta, dltmp)
endif

endif

elseif(retcd .ne. 2 .and. (abs(pred-ared) .le. Rpten*abs(ared)) &
& .and. (.not. nwtstep) .and. (delta .le. Rp99*stepmx)) then

!        pred predicts ared very well
!        attempt a doubling of the trust region and continue global step
!        when not taking a newton step and trust region not at maximum

call dcopy(n,xp,1,xprev,1)
call dcopy(n,fp,1,fprev,1)
fpnsav = fpnorm
delta  = min(Rtwo*delta,stepmx)
retcd  = 3

else

!        fpnorm sufficiently small to accept xp as next iterate.
!        Choose new trust region.

retcd = 0
if(ared .ge. Rpten*pred) then

!           Not good enough. Decrease trust region for next iteration

delta = Rhalf*delta
elseif( ared .le. Rp75*pred ) then

!           Wonderful. Increase trust region for next iteration

delta = min(Rtwo*delta,stepmx)
endif

endif

return
end

! ============================================================================
! From converted legacy source: nwutil.f90
! ============================================================================

subroutine nwtcvg(xplus,fplus,xc,xtol,retcd,ftol,iter, &
& maxit,n,ierr,termcd)

integer n,iter,maxit,ierr,termcd,retcd
double precision  xtol,ftol
double precision  xplus(*),fplus(*),xc(*)

!-------------------------------------------------------------------------
!
!     Decide whether to terminate the nonlinear algorithm
!
!     Arguments
!
!     In       xplus   Real(*)         new x values
!     In       fplus   Real(*)         new f values
!     In       xc      Real(*)         current x values
!     In       xtol    Real            stepsize tolerance
!     In       retcd   Integer         return code from global search routines
!     In       ftol    Real            function tolerance
!     In       iter    Integer         iteration number
!     In       maxit   Integer         maximum number of iterations allowed
!     In       n       Integer         size of x and f
!     In       ierr    Integer         return code of cndjac (condition estimation)
!
!     Out      termcd  Integer         termination code
!                                        0 no termination criterion satisfied
!                                          ==> continue iterating
!                                        1 norm of scaled function value is
!                                          less than ftol
!                                        2 scaled distance between last
!                                          two steps less than xtol
!                                        3 unsuccessful global strategy
!                                          ==> cannot find a better point
!                                        4 iteration limit exceeded
!                                        5 Jacobian too ill-conditioned
!                                        6 Jacobian singular
!                                        7 Jacobian not usable (all zero entries)
!
!-------------------------------------------------------------------------

double precision fmax, rsx


!     check whether function values are within tolerance

termcd = 0

if( ierr .ne. 0 ) then
termcd = 4 + ierr
return
endif

fmax = abs(fplus(idamax(n,fplus,1)))
if( fmax .le. ftol) then
termcd = 1
return
endif

!     initial check at start so there is no xplus
!     so only a check of function values is useful
if(iter .eq. 0) return

if(retcd .eq. 1) then
termcd = 3
return
endif

!     check whether relative step length is within tolerance
!     Dennis Schnabel Algorithm A7.2.3

rsx = nuxnrm(n, xplus, xc)
if(rsx .le. xtol) then
termcd = 2
return
endif

!     check iteration limit

if(iter .ge. maxit) then
termcd = 4
endif

return
end

!-----------------------------------------------------------------------

subroutine nweset(n,xc,fc,fcnorm,xp,fp,fpnorm,gcnt,priter,iter)
double precision xc(*),fc(*),fcnorm,xp(*),fp(*),fpnorm
integer n, gcnt, priter, iter

!-------------------------------------------------------------------------
!
!     calling routine got an error in decomposition/update of Jacobian/Broyden
!     jacobian an singular or too ill-conditioned
!     prepare return arguments
!
!     Arguments
!
!     In       n       Integer         size of x
!     In       xc      Real(*)         current (starting) x values
!     In       fc      Real(*)         function values f(xc)
!     In       fcnorm  Real            norm fc
!     Out      xp      Real(*)         final x values
!     Out      fp      Real(*)         function values f(xp)
!     Out      fpnorm  Real            final norm fp
!     Out      gcnt    Integer         # of backtracking steps (here set to 0)
!     In       priter  Integer         flag for type of output
!     In       iter    Integer         iteration counter
!
!-------------------------------------------------------------------------

call dcopy(n,xc,1,xp,1)
call dcopy(n,fc,1,fp,1)
fpnorm = fcnorm
gcnt   = 0
if( priter .gt. 0 ) then
call nwjerr(iter)
endif

return
end

!-----------------------------------------------------------------------

subroutine chkjac1(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,fvec,termcd)

integer lda,n,termcd
double precision  A(lda,*),xc(*),fc(*)
double precision  epsm,scalex(*)
double precision  fz(*),wa(*),xw(*)
procedure(legacy_fvec) :: fvec

!-------------------------------------------------------------------------
!
!     Check the user supplied jacobian against its finite difference approximation
!
!     Arguments
!
!     In       A       Real(lda,*)     user supplied jacobian
!     In       lda     Integer         leading dimension of ajanal
!     In       xc      Real(*)         vector of x values
!     In       fc      Real(*)         function values f(xc)
!     In       n       Integer         size of x
!     In       epsm    Real            machine precision
!     In       scalex  Real(*)         scaling vector for x()
!     Wk       fz      Real(*)         workspace
!     Wk       wa      Real(*)         workspace
!     Wk       xw      Real(*)         workspace
!     In       fvec    Name            name of routine to evaluate f(x)
!     Out      termcd  Integer         return code
!                                        0  user supplied jacobian ok
!                                      -10  user supplied jacobian NOT ok
!
!-------------------------------------------------------------------------

integer i,j,errcnt
double precision  ndigit,p,h,xcj,dinf
double precision  tol



integer MAXERR
parameter(MAXERR=10)

double precision Rquart, Rten
parameter(Rquart=0.25d0, Rten=10.0d0)

termcd = 0

!     compute the finite difference jacobian and check it against
!     the analytic one

ndigit = -log10(epsm)
p = sqrt(max(Rten**(-ndigit),epsm))
tol    = epsm**Rquart

errcnt = 0
call dcopy(n,xc,1,xw,1)
call vunsc(n,xw,scalex)

do j=1,n
h = p + p * abs(xw(j))
xcj   = xw(j)
xw(j) = xcj + h

!        avoid (small) rounding errors
!        h = xc(j) - xcj but not here to avoid clever optimizers

h = rnudif(xw(j), xcj)

call fvec(xw,fz,n,j)
xw(j) = xcj

do i=1,n
wa(i) = (fz(i)-fc(i))/h
enddo

dinf = abs(wa(idamax(n,wa,1)))

do i=1,n
if(abs(A(i,j)-wa(i)).gt.tol*dinf) then
errcnt = errcnt + 1
if( errcnt .gt. MAXERR ) then
termcd = -10
return
endif
call nwckot(i,j,A(i,j),wa(i))
endif
enddo
enddo

!      call vscal(n,xc,scalex)

if( errcnt .gt. 0 ) then
termcd = -10
endif
return
end

!-----------------------------------------------------------------------

subroutine chkjac2(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,fvec,termcd, &
& dsub,dsuper)

integer lda,n,termcd,dsub,dsuper
double precision  A(lda,*),xc(*),fc(*)
double precision  epsm,scalex(*)
double precision  fz(*),wa(*),xw(*)
procedure(legacy_fvec) :: fvec

!-------------------------------------------------------------------------
!
!     Check the user supplied jacobian against its finite difference approximation
!
!     Arguments
!
!     In       A       Real(lda,*)     user supplied jacobian
!     In       lda     Integer         leading dimension of ajanal
!     In       xc      Real(*)         vector of x values
!     In       fc      Real(*)         function values f(xc)
!     In       n       Integer         size of x
!     In       epsm    Real            machine precision
!     In       scalex  Real(*)         scaling vector for x()
!     Wk       fz      Real(*)         workspace
!     Wk       wa      Real(*)         workspace
!     Wk       xw      Real(*)         workspace
!     In       fvec    Name            name of routine to evaluate f(x)
!     Out      termcd  Integer         return code
!                                        0  user supplied jacobian ok
!                                      -10  user supplied jacobian NOT ok
!
!-------------------------------------------------------------------------

integer i,j,k,dsum,errcnt
double precision  ndigit,p,h,dinf
double precision  tol
double precision w(n),xstep(n)

integer MAXERR
parameter(MAXERR=10)

double precision Rquart, Rten, Rzero
parameter(Rquart=0.25d0, Rten=10.0d0, Rzero=0.0d0)

dsum = dsub + dsuper + 1

termcd = 0

!     compute the finite difference jacobian and check it against
!     the user supplied one

ndigit = -log10(epsm)
p = sqrt(max(Rten**(-ndigit),epsm))
tol    = epsm**Rquart

errcnt = 0
call dcopy(n,xc,1,xw,1)
call vunsc(n,xw,scalex)

do j=1,n
xstep(j) = p + p * abs(xw(j))
w(j) = xw(j)
enddo

do k=1,dsum
do j=k,n,dsum
xw(j) = xw(j) + xstep(j)
enddo

!        for non finite values error message will be wrong
call fvec(xw,fz,n,n+k)

do j=k,n,dsum
h = xstep(j)
xw(j) = w(j)
dinf = Rzero
do i=max(j-dsuper,1),min(j+dsub,n)
wa(i) = (fz(i)-fc(i)) / h
if(abs(wa(i)).gt.dinf) dinf = abs(wa(i))
enddo

do i=max(j-dsuper,1),min(j+dsub,n)
if(abs(A(i,j)-wa(i)).gt.tol*dinf) then
errcnt = errcnt + 1
if( errcnt .gt. MAXERR ) then
termcd = -10
return
endif
call nwckot(i,j,A(i,j),wa(i))
endif
enddo
enddo
enddo

!      call vscal(n,xc,scalex)

if( errcnt .gt. 0 ) then
termcd = -10
endif
return
end

!-----------------------------------------------------------------------

subroutine chkjac(A,lda,xc,fc,n,epsm,jacflg,scalex,fz,wa,xw,fvec, &
& termcd)

integer lda,n,termcd,jacflg(*)
double precision  A(lda,*),xc(*),fc(*)
double precision  epsm,scalex(*)
double precision  fz(*),wa(*),xw(*)
procedure(legacy_fvec) :: fvec

!-------------------------------------------------------------------------
!
!     Check the user supplied jacobian against its finite difference approximation
!
!     Arguments
!
!     In       A       Real(lda,*)     user supplied jacobian
!     In       lda     Integer         leading dimension of ajanal
!     In       xc      Real(*)         vector of x values
!     In       fc      Real(*)         function values f(xc)
!     In       n       Integer         size of x
!     In       epsm    Real            machine precision
!     In       jacflg  Integer(*)      indicates how to compute jacobian
!                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
!                                                  3: user supplied banded
!                                      jacflg[2]: number of sub diagonals or -1 if not banded
!                                      jacflg[3]: number of super diagonals or -1 if not banded
!                                      jacflg[4]: 1 if adjusting jacobian allowed when
!                                                   singular or illconditioned
!     In       scalex  Real(*)         scaling vector for x()
!     Wk       fz      Real(*)         workspace
!     Wk       wa      Real(*)         workspace
!     Wk       xw      Real(*)         workspace
!     In       fvec    Name            name of routine to evaluate f(x)
!     Out      termcd  Integer         return code
!                                        0  user supplied jacobian ok
!                                      -10  user supplied jacobian NOT ok
!
!-------------------------------------------------------------------------

if(jacflg(1) .eq. 3) then
!        user supplied and banded
call chkjac2(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,fvec,termcd, &
& jacflg(2),jacflg(3))
else
call chkjac1(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,fvec,termcd)
endif

return
end

!-----------------------------------------------------------------------

subroutine fdjac0(xc,fc,n,epsm,fvec,fz,rjac,ldr)

integer ldr,n
double precision  epsm
double precision  rjac(ldr,*),fz(*),xc(*),fc(*)
procedure(legacy_fvec) :: fvec

!-------------------------------------------------------------------------
!
!     Compute the finite difference jacobian at the current point xc
!
!     Arguments
!
!     In       xc      Real(*)         current point
!     In       fc      Real(*)         function values at current point
!     In       n       Integer         size of x and f
!     In       epsm    Real            machine precision
!     In       fvec    Name            name of routine to evaluate f(x)
!     Wk       fz      Real(*)         workspace
!     Out      rjac    Real(ldr,*)     jacobian matrix at x
!                                        entry [i,j] is derivative of
!                                        f(i) wrt to x(j)
!     In       ldr     Integer         leading dimension of rjac
!
!-------------------------------------------------------------------------

integer i,j
double precision  ndigit,p,h,xcj


double precision Rten
parameter(Rten=10d0)

ndigit = -log10(epsm)
p = sqrt(max(Rten**(-ndigit),epsm))

do j=1,n
h = p + p * abs(xc(j))

!        or as alternative h  = p * max(Rone, abs(xc(j)))

xcj   = xc(j)
xc(j) = xcj + h

!        avoid (small) rounding errors
!        h = xc(j) - xcj  but not here to avoid clever optimizers

h = rnudif(xc(j), xcj)
call fvec(xc,fz,n,j)
xc(j) = xcj
do i=1,n
rjac(i,j) = (fz(i)-fc(i)) / h
enddo
enddo

return
end

!-----------------------------------------------------------------------

subroutine fdjac2(xc,fc,n,epsm,fvec,fz,rjac,ldr,dsub,dsuper, &
& w,xstep)

integer ldr,n,dsub,dsuper
double precision  epsm
double precision  rjac(ldr,*),fz(*),xc(*),fc(*)
double precision  w(*), xstep(*)
procedure(legacy_fvec) :: fvec

!-------------------------------------------------------------------------
!
!     Compute a banded finite difference jacobian at the current point xc
!
!     Arguments
!
!     In       xc      Real(*)         current point
!     In       fc      Real(*)         function values at current point
!     In       n       Integer         size of x and f
!     In       epsm    Real            machine precision
!     In       fvec    Name            name of routine to evaluate f(x)
!     Wk       fz      Real(*)         workspace
!     Out      rjac    Real(ldr,*)     jacobian matrix at x
!                                        entry [i,j] is derivative of
!                                        f(i) wrt to x(j)
!     In       ldr     Integer         leading dimension of rjac
!     In       dsub    Integer         number of subdiagonals
!     In       dsuper  Integer         number of superdiagonals
!     Internal w       Real(*)         for temporary saving of xc
!     Internal xstep   Real(*)         stepsizes
!
!-------------------------------------------------------------------------

integer i,j,k
double precision  ndigit,p,h


double precision Rten
parameter(Rten=10d0)

integer dsum

dsum = dsub + dsuper + 1

ndigit = -log10(epsm)
p = sqrt(max(Rten**(-ndigit),epsm))

do k=1,n
xstep(k) = p + p * abs(xc(k))
enddo

do k=1,dsum
do j=k,n,dsum
w(j) = xc(j)
xc(j) = xc(j) + xstep(j)
enddo

call fvec(xc,fz,n,n+k)
do j=k,n,dsum
call nuzero(n,rjac(1,j))
!            fdjac0 for why
!            doing this ensures that results for fdjac2 and fdjac0 will be identical
h = rnudif(xc(j),w(j))
xc(j) = w(j)
do i=max(j-dsuper,1),min(j+dsub,n)
rjac(i,j) = (fz(i)-fc(i)) / h
enddo
enddo
enddo

return
end

!-----------------------------------------------------------------------

double precision function nudnrm(n, d, x)
integer n
double precision  d(*), x(*)


!-------------------------------------------------------------------------
!
!     calculate  max( abs(d[*]) / max(x[*],1) )
!
!     Arguments
!
!     In   n        Integer       number of elements in d() and x()
!     In   d        Real(*)       vector d
!     In   x        Real(*)       vector x
!
!-------------------------------------------------------------------------

integer i
double precision  t

double precision Rzero, Rone
parameter(Rzero=0.0d0, Rone=1.0d0)

t = Rzero
do i=1,n
t = max(t, abs(d(i)) / max(abs(x(i)),Rone))
enddo
nudnrm = t

return
end

!-----------------------------------------------------------------------

double precision function nuxnrm(n, xn, xc)
integer n
double precision  xn(*), xc(*)


!-------------------------------------------------------------------------
!
!     calculate  max( abs(xn[*]-xc[*]) / max(xn[*],1) )
!
!     Arguments
!
!     In   n        Integer       number of elements in xn() and xc()
!     In   xn       Real(*)       vector xn
!     In   xc       Real(*)       vector xc
!
!-------------------------------------------------------------------------

integer i
double precision  t

double precision Rzero, Rone
parameter(Rzero=0.0d0, Rone=1.0d0)

t = Rzero
do i=1,n
t = max(t, abs(xn(i)-xc(i)) / max(abs(xn(i)),Rone))
enddo
nuxnrm = t

return
end

!-----------------------------------------------------------------------

double precision function rnudif(x, y)
double precision x, y


!-------------------------------------------------------------------------
!
!     Return difference of x and y (x - y)
!
!     Arguments
!
!     In   x  Real      argument 1
!     In   y  Real      argument 2
!
!-------------------------------------------------------------------------

rnudif = x - y
return
end

!-----------------------------------------------------------------------

subroutine compmu(r,ldr,n,mu,y,ierr)

integer ldr,n,ierr
double precision r(ldr,*),mu,y(*)

!-------------------------------------------------------------------------
!
!     Compute a small perturbation mu for the (almost) singular matrix R.
!     mu is used in the computation of the Levenberg-Marquardt step.
!
!     Arguments
!
!     In       R       Real(ldr,*)     upper triangular matrix from QR
!     In       ldr     Integer         leading dimension of R
!     In       n       Integer         column dimension of R
!     Out      mu      Real            sqrt(l1 norm of R * infinity norm of R
!                                      * n * epsm * 100) designed to make
!                                        trans(R)*R + mu * I not singular
!     Wk       y       Real(*)         workspace for dlange
!     Out      ierr    Integer         0 indicating mu ok
!                                      3 indicating mu much too small
!
!-------------------------------------------------------------------------

double precision aifnrm,al1nrm,epsm



double precision Rhund
parameter(Rhund=100d0)

!     get the infinity norm of R
!     get the l1 norm of R
ierr = 0
aifnrm = dlantr('I','U','N',n,n,r,ldr,y)
al1nrm = dlantr('1','U','N',n,n,r,ldr,y)
epsm = epsmch()
mu = sqrt(n*epsm*Rhund)*aifnrm*al1nrm
!     matrix consists of zero's or near zero's
!     LM correction in liqrev will not work
if( mu .le. Rhund*epsm ) then
ierr = 3
endif
return
end

!-----------------------------------------------------------------------

subroutine cndjac(n,r,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)
integer n,ldr,icdwrk(*),ierr
double precision cndtol,rcond,r(ldr,*),rcdwrk(*)

!---------------------------------------------------------------------
!
!     Check r for singularity and/or ill conditioning
!
!     Arguments
!
!     In       n       Integer         dimension of problem
!     In       r       Real(ldr,*)     upper triangular R from QR decomposition
!     In       ldr     Integer         leading dimension of rjac
!     In       cndtol  Real            tolerance of test for ill conditioning
!                                       when rcond <= cndtol then ierr is set to 1
!                                       cndtol should be >= machine precision
!     Out      rcond   Real            inverse condition  of r
!     Wk       rcdwrk  Real(*)         workspace (for dtrcon)
!     Wk       icdwrk  Integer(*)      workspace (fordtrcon)
!     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
!                                      1 indicating Jacobian too ill-conditioned
!                                      2 indicating Jacobian completely singular
!
!---------------------------------------------------------------------

integer i,info
logical rsing
double precision Rzero
parameter(Rzero=0.0d0)

ierr = 0

rsing = .false.
do i=1,n
if( r(i,i) .eq. Rzero ) then
rsing = .true.
endif
enddo

if( rsing ) then
ierr = 2
rcond = Rzero
else
call dtrcon('1','U','N',n,r,ldr,rcond,rcdwrk,icdwrk,info)
if( rcond .eq. Rzero ) then
ierr = 2
elseif( rcond .le. cndtol ) then
ierr = 1
endif
endif

return
end

!-----------------------------------------------------------------------

subroutine nwfjac(x,scalex,f,fq,n,epsm,jacflg,fvec,mkjac,rjac, &
& ldr,xw,w,xstep)

integer ldr,n,jacflg(*)
double precision  epsm
double precision  x(*),f(*),scalex(*),xw(*),w(*),xstep(*)
double precision  rjac(ldr,*),fq(*)
procedure(legacy_fvec) :: fvec
procedure(legacy_fjac) :: mkjac

!-------------------------------------------------------------------------
!
!     Calculate the jacobian  matrix
!
!     Arguments
!
!     In       x       Real(*)         (scaled) current x values
!     In       scalex  Real(*)         scaling factors x
!     In       f       Real(*)         function values f(x)
!     Wk       fq      Real(*)         (internal) workspace
!     In       n       Integer         size of x and f
!     In       epsm    Real            machine precision
!     In       jacflg  Integer(*)      indicates how to compute jacobian
!                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
!                                                  3: user supplied banded
!                                      jacflg[2]: number of sub diagonals or -1 if not banded
!                                      jacflg[3]: number of super diagonals or -1 if not banded
!                                      jacflg[4]: 1 if adjusting jacobian allowed when
!                                                   singular or illconditioned
!     In       fvec    Name            name of routine to evaluate f()
!     In       mkjac   Name            name of routine to evaluate jacobian
!     Out      rjac    Real(ldr,*)     jacobian matrix (unscaled)
!     In       ldr     Integer         leading dimension of rjac
!     Internal xw      Real(*)         used for storing unscaled x
!     Internal w       Real(*)         workspace for banded jacobian
!     Internal xstep   Real(*)         workspace for banded jacobian
!
!-------------------------------------------------------------------------

!     compute the finite difference or analytic jacobian at x

call dcopy(n,x,1,xw,1)
call vunsc(n,xw,scalex)
if(jacflg(1) .eq. 0) then
call fdjac0(xw,f,n,epsm,fvec,fq,rjac,ldr)
elseif(jacflg(1) .eq. 2) then
call fdjac2(xw,f,n,epsm,fvec,fq,rjac,ldr,jacflg(2),jacflg(3), &
& w,xstep)
else
call mkjac(rjac,ldr,xw,n)
endif

return
end

!-----------------------------------------------------------------------

subroutine nwscjac(n,rjac,ldr,scalex)
integer n, ldr
double precision rjac(ldr,*), scalex(*)

!-------------------------------------------------------------------------
!
!     Scale jacobian
!
!     Arguments
!
!     In       n       Integer         size of x and f
!     Inout    rjac    Real(ldr,*)     jacobian matrix
!     In       ldr     Integer         leading dimension of rjac
!     In       scalex  Real(*)         scaling factors for x
!
!-------------------------------------------------------------------------

integer j
double precision t, Rone
parameter(Rone=1.0d0)

do j = 1,n
t = Rone/scalex(j)
call dscal(n,t,rjac(1,j),1)
enddo

return
end

!-----------------------------------------------------------------------

subroutine nwunscjac(n,rjac,ldr,scalex)
integer n, ldr
double precision rjac(ldr,*), scalex(*)

!-------------------------------------------------------------------------
!
!     Unscale jacobian
!
!     Arguments
!
!     In       n       Integer         size of x and f
!     Inout    rjac    Real(ldr,*)     jacobian matrix
!     In       ldr     Integer         leading dimension of rjac
!     In       scalex  Real(*)         scaling factors for x
!
!-------------------------------------------------------------------------

integer j
double precision t

do j = 1,n
t = scalex(j)
call dscal(n,t,rjac(1,j),1)
enddo

return
end

!-----------------------------------------------------------------------

subroutine nwcpsx(n,rjac,ldr,scalex,epsm, mode)

integer ldr,n,mode
double precision  epsm
double precision  scalex(*)
double precision  rjac(ldr,*)

!-------------------------------------------------------------------------
!
!     Calculate scaling factors from the jacobian  matrix
!
!     Arguments
!
!     In       n       Integer         size of x and f
!     In       rjac    Real(ldr,*)     jacobian matrix
!     In       ldr     Integer         leading dimension of rjac
!     Out      scalex  Real(*)         scaling factors for x
!     In       epsm    Real            machine precision
!     In       mode    Integer         1: initialise, >1: adjust
!-------------------------------------------------------------------------

integer k


if( mode .eq. 1 ) then
do k=1,n
scalex(k) = dnrm2(n,rjac(1,k),1)
if( scalex(k) .le. epsm ) scalex(k) = 1
enddo
else if( mode .gt. 1 ) then
do k=1,n
scalex(k) = max(scalex(k),dnrm2(n,rjac(1,k),1))
enddo
endif
return
end

!-----------------------------------------------------------------------

subroutine nwcpmt(n, x, scalex, factor, wrk, stepsiz)
double precision x(*), scalex(*), wrk(*)
double precision factor, stepsiz
integer n

!-------------------------------------------------------------------------
!
!     Calculate maximum stepsize
!
!     Arguments
!
!     In       n       Integer     size of x
!     In       x       Real(*)     x-values
!     In       scalex  Real(*)     scaling factors for x
!     In       factor  Real        multiplier
!     Inout    wrk     Real(*)     workspace
!     Out      stepsiz Real        stepsize
!
!     Currently not used
!     Minpack uses this to calculate initial trust region size
!     Not (yet) used in this code because it doesn't seem to help
!     Manually setting an initial trust region size works better
!
!-------------------------------------------------------------------------

double precision Rzero
parameter(Rzero=0.0d0)



call dcopy(n,x,1,wrk,1)
call vscal(n,wrk,scalex)
stepsiz = factor * dnrm2(n,wrk,1)
if( stepsiz .eq. Rzero ) stepsiz = factor
return
end

!-----------------------------------------------------------------------

subroutine vscal(n,x,sx)

integer n
double precision  x(*),sx(*)

!-------------------------------------------------------------------------
!
!     Scale a vector x
!
!     Arguments
!
!     In       n       Integer         size of x
!     Inout    x       Real(*)         vector to scale
!     In       sx      Real(*)         scaling vector
!
!-------------------------------------------------------------------------

integer i

do i = 1,n
x(i) = sx(i) * x(i)
enddo

return
end

!-----------------------------------------------------------------------

subroutine vunsc(n,x,sx)

integer n
double precision  x(*),sx(*)

!-------------------------------------------------------------------------
!
!     Unscale a vector x
!
!     Arguments
!
!     In       n       Integer         size of x
!     Inout    x       Real(*)         vector to unscale
!     In       sx      Real(*)         scaling vector
!
!-------------------------------------------------------------------------

integer i

do i = 1,n
x(i) = x(i) / sx(i)
enddo

return
end

!-----------------------------------------------------------------------

subroutine nwfvec(x,n,scalex,fvec,f,fnorm,xw)

integer n
double precision  x(*),xw(*),scalex(*),f(*),fnorm
procedure(legacy_fvec) :: fvec

!-------------------------------------------------------------------------
!
!     Evaluate the function at current iterate x and scale its value
!
!     Arguments
!
!     In       x       Real(*)         x
!     In       n       Integer         size of x
!     In       scalex  Real(*)         scaling vector for x
!     In       fvec    Name            name of routine to calculate f(x)
!     Out      f       Real(*)         f(x)
!     Out      fnorm   Real            .5*||f(x)||**2
!     Internal xw      Real(*)         used for storing unscaled xc
!
!-------------------------------------------------------------------------



double precision Rhalf
parameter(Rhalf=0.5d0)

call dcopy(n,x,1,xw,1)
call vunsc(n,xw,scalex)
call fvec(xw,f,n,0)

fnorm = Rhalf * dnrm2(n,f,1)**2

return
end

!-----------------------------------------------------------------------

double precision function epsmch()

!     Return machine precision
!     Use Lapack routine





!     dlamch('e') returns negeps (1-eps)
!     dlamch('p') returns 1+eps

epsmch = dlamch('p')

return
end

!-----------------------------------------------------------------------

double precision function dblhuge()

!     Return largest double precision number
!     Use Lapack routine





!     dlamch('o') returns max double precision

dblhuge = dlamch('o')

return
end

end module nleqslv_legacy_kernels
