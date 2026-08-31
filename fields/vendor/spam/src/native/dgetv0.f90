!\BeginDoc
!
!\Name: dgetv0
!
!\Description:
!  Generate a random initial residual vector for the Arnoldi process.
!  Force the residual vector to be in the range of the operator OP.
!
!\Usage:
!  call dgetv0
!     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM,
!       IPNTR, WORKD, IERR )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first
!          call to dgetv0.
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO = 99: done
!          -------------------------------------------------------------
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B in the (generalized)
!          eigenvalue problem A*x = lambda*B*x.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  ITRY    Integer.  (INPUT)
!          ITRY counts the number of times that dgetv0 is called.
!          It should be set to 1 on the initial call to dgetv0.
!
!  INITV   Logical variable.  (INPUT)
!          .TRUE.  => the initial residual vector is given in RESID.
!          .FALSE. => generate a random initial residual vector.
!
!  N       Integer.  (INPUT)
!          Dimension of the problem.
!
!  J       Integer.  (INPUT)
!          Index of the residual vector to be generated, with respect to
!          the Arnoldi process.  J > 1 in case of a "restart".
!
!  V       Double precision N by J array.  (INPUT)
!          The first J-1 columns of V contain the current Arnoldi basis
!          if this is a "restart".
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  RESID   Double precision array of length N.  (INPUT/OUTPUT)
!          Initial residual vector to be generated.  If RESID is
!          provided, force RESID into the range of the operator OP.
!
!  RNORM   Double precision scalar.  (OUTPUT)
!          B-norm of the generated residual.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!
!  WORKD   Double precision work array of length 2*N.  (REVERSE COMMUNICATION).
!          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
!
!  IERR    Integer.  (OUTPUT)
!          =  0: Normal exit.
!          = -1: Cannot generate a nontrivial restarted residual vector
!                in the range of the operator OP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!
!\Routines called:
!     second  ARPACK utility routine for timing.
!     dvout   ARPACK utility routine for vector output.
!     dlarnv  LAPACK routine for generating a random vector.
!     dgemv   Level 2 BLAS routine for matrix vector multiplication.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     ddot    Level 1 BLAS that computes the scalar product of two vectors.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: getv0.F   SID: 2.7   DATE OF SID: 04/07/99   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dgetv0 &
& ( ido, bmat, initv, n, j, v, ldv, resid, rnorm, &
& ipntr, workd, ierr )
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%

implicit none
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
character  bmat*1
logical    initv
integer    ido, ierr,  j, ldv, n
Double precision &
& rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
integer    ipntr(3)
Double precision &
& resid(n), v(ldv,j), workd(2*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
Double precision &
& one, zero
parameter (one = 1.0D+0, zero = 0.0D+0)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
logical    first, inits, orth
integer    idist, iseed(4), iter,  jj
Double     precision &
& rnorm0
save       first, iseed, inits, iter, orth, rnorm0
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
external   dlarnv, dcopy, dgemv
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
Double precision &
& ddot, dnrm2
external   ddot, dnrm2
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
intrinsic    abs, sqrt
!
!     %-----------------%
!     | Data Statements |
!     %-----------------%
!
data       inits /.true./
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!
!     %-----------------------------------%
!     | Initialize the seed of the LAPACK |
!     | random number generator           |
!     %-----------------------------------%
!
if (inits) then
iseed(1) = 1
iseed(2) = 3
iseed(3) = 5
iseed(4) = 7
inits = .true.
end if
!
if (ido .eq.  0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
ierr   = 0
iter   = 0
first  = .FALSE.
orth   = .FALSE.
!
!        %-----------------------------------------------------%
!        | Possibly generate a random starting vector in RESID |
!        | Use a LAPACK random number generator used by the    |
!        | matrix generation routines.                         |
!        |    idist = 1: uniform (0,1)  distribution;          |
!        |    idist = 2: uniform (-1,1) distribution;          |
!        |    idist = 3: normal  (0,1)  distribution;          |
!        %-----------------------------------------------------%
!
if (.not.initv) then
idist = 2
call dlarnv (idist, iseed, n, resid)
end if
!
!        %----------------------------------------------------------%
!        | Force the starting vector into the range of OP to handle |
!        | the generalized problem when B is possibly (singular).   |
!        %----------------------------------------------------------%
!
if (bmat .eq. 'G') then
ipntr(1) = 1
ipntr(2) = n + 1
call dcopy (n, resid, 1, workd, 1)
ido = -1
go to 9000
end if
end if
!
!     %-----------------------------------------%
!     | Back from computing OP*(initial-vector) |
!     %-----------------------------------------%
!
if (first) go to 20
!
!     %-----------------------------------------------%
!     | Back from computing B*(orthogonalized-vector) |
!     %-----------------------------------------------%
!
if (orth)  go to 40
!
!     %------------------------------------------------------%
!     | Starting vector is now in the range of OP; r = OP*r; |
!     | Compute B-norm of starting vector.                   |
!     %------------------------------------------------------%
!
first = .TRUE.
if (bmat .eq. 'G') then
call dcopy (n, workd(n+1), 1, resid, 1)
ipntr(1) = n + 1
ipntr(2) = 1
ido = 2
go to 9000
else if (bmat .eq. 'I') then
call dcopy (n, resid, 1, workd(:), 1)
end if
!
20 continue
!
first = .FALSE.
if (bmat .eq. 'G') then
rnorm0 = ddot (n, resid, 1, workd, 1)
rnorm0 = sqrt(abs(rnorm0))
else if (bmat .eq. 'I') then
rnorm0 = dnrm2(n, resid, 1)
end if
rnorm  = rnorm0
!
!     %---------------------------------------------%
!     | Exit if this is the very first Arnoldi step |
!     %---------------------------------------------%
!
if (j .eq. 1) go to 50
!
!     %----------------------------------------------------------------
!     | Otherwise need to B-orthogonalize the starting vector against |
!     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
!     | This is the case where an invariant subspace is encountered   |
!     | in the middle of the Arnoldi factorization.                   |
!     |                                                               |
!     |       s = V^{T}*B*r;   r = r - V*s;                           |
!     |                                                               |
!     | Stopping criteria used for iter. ref. is discussed in         |
!     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
!     %---------------------------------------------------------------%
!
orth = .TRUE.
30 continue
!
call dgemv ('T', n, j-1, one, v, ldv, workd(:), 1, &
& zero, workd(n+1), 1)
call dgemv ('N', n, j-1, -one, v, ldv, workd(n+1), 1, &
& one, resid(n), 1)
!
!     %----------------------------------------------------------%
!     | Compute the B-norm of the orthogonalized starting vector |
!     %----------------------------------------------------------%
!
if (bmat .eq. 'G') then
call dcopy (n, resid, 1, workd(n+1), 1)
ipntr(1) = n + 1
ipntr(2) = 1
ido = 2
go to 9000
else if (bmat .eq. 'I') then
call dcopy (n, resid, 1, workd, 1)
end if
!
40 continue
!
if (bmat .eq. 'G') then
rnorm = ddot (n, resid, 1, workd, 1)
rnorm = sqrt(abs(rnorm))
else if (bmat .eq. 'I') then
rnorm = dnrm2(n, resid, 1)
end if
!
!     %--------------------------------------%
!     | Check for further orthogonalization. |
!     %--------------------------------------%
!
if (rnorm .gt. 0.717*rnorm0) go to 50
!
iter = iter + 1
if (iter .le. 5) then
!
!        %-----------------------------------%
!        | Perform iterative refinement step |
!        %-----------------------------------%
!
rnorm0 = rnorm
go to 30
else
!
!        %------------------------------------%
!        | Iterative refinement step "failed" |
!        %------------------------------------%
!
do 45 jj = 1, n
resid(jj) = zero
45 continue
rnorm = zero
ierr = -1
end if
!
50 continue
!
ido = 99
!
9000 continue
return
!
!     %---------------%
!     | End of dgetv0 |
!     %---------------%
!
end
