!     Modifications:
!     2023-10-17: eliminated all kind=4 occurences.
!
!      
! dneupd.f
! dnaupd.f
! dnaup2.f
! dnapps.f
! dnconv.f
! dsortc.f
! dneigh.f
! dlaqrb.f
! dngets.f
! dgetv0.f
! dnaitr.f
!


!\BeginDoc
!
!\Name: dneupd
!
!\Description:
!
!  This subroutine returns the converged approximations to eigenvalues
!  of A*z = lambda*B*z and (optionally):
!
!      (1) The corresponding approximate eigenvectors;
!
!      (2) An orthonormal basis for the associated approximate
!          invariant subspace;
!
!      (3) Both.
!
!  There is negligible additional cost to obtain eigenvectors.  An orthonormal
!  basis is always computed.  There is an additional storage cost of n*nev
!  if both are requested (in this case a separate array Z must be supplied).
!
!  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
!  are derived from approximate eigenvalues and eigenvectors of
!  of the linear operator OP prescribed by the MODE selection in the
!  call to DNAUPD .  DNAUPD  must be called before this routine is called.
!  These approximate eigenvalues and vectors are commonly called Ritz
!  values and Ritz vectors respectively.  They are referred to as such
!  in the comments that follow.  The computed orthonormal basis for the
!  invariant subspace corresponding to these Ritz values is referred to as a
!  Schur basis.
!
!  See documentation in the header of the subroutine DNAUPD  for
!  definition of OP as well as other terms and the relation of computed
!  Ritz values and Ritz vectors of OP with respect to the given problem
!  A*z = lambda*B*z.  For a brief description, see definitions of
!  IPARAM(7), MODE and WHICH in the documentation of DNAUPD .
!
!\Usage:
!  call dneupd
!     ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, WORKEV, BMAT,
!       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL,
!       LWORKL, INFO )
!
!\Arguments:
!  RVEC    LOGICAL  (INPUT)
!          Specifies whether a basis for the invariant subspace corresponding
!          to the converged Ritz value approximations for the eigenproblem
!          A*z = lambda*B*z is computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
!                                See Remarks below.
!
!  HOWMNY  Character*1  (INPUT)
!          Specifies the form of the basis for the invariant subspace
!          corresponding to the converged Ritz values that is to be computed.
!
!          = 'A': Compute NEV Ritz vectors;
!          = 'P': Compute NEV Schur vectors;
!          = 'S': compute some of the Ritz vectors, specified
!                 by the logical array SELECT.
!
!  SELECT  Logical array of dimension NCV.  (INPUT)
!          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
!          computed. To select the Ritz vector corresponding to a
!          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE..
!          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
!
!  DR      Double precision  array of dimension NEV+1.  (OUTPUT)
!          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains
!          the real part of the Ritz  approximations to the eigenvalues of
!          A*z = lambda*B*z.
!          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit:
!          DR contains the real part of the Ritz values of OP computed by
!          DNAUPD . A further computation must be performed by the user
!          to transform the Ritz values computed for OP by DNAUPD  to those
!          of the original system A*z = lambda*B*z. See remark 3 below.
!
!  DI      Double precision  array of dimension NEV+1.  (OUTPUT)
!          On exit, DI contains the imaginary part of the Ritz value
!          approximations to the eigenvalues of A*z = lambda*B*z associated
!          with DR.
!
!          NOTE: When Ritz values are complex, they will come in complex
!                conjugate pairs.  If eigenvectors are requested, the
!                corresponding Ritz vectors will also come in conjugate
!                pairs and the real and imaginary parts of these are
!                represented in two consecutive columns of the array Z
!                (see below).
!
!  Z       Double precision  N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
!          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
!          Z represent approximate eigenvectors (Ritz vectors) corresponding
!          to the NCONV=IPARAM(5) Ritz values for eigensystem
!          A*z = lambda*B*z.
!
!          The complex Ritz vector associated with the Ritz value
!          with positive imaginary part is stored in two consecutive
!          columns.  The first column holds the real part of the Ritz
!          vector and the second column holds the imaginary part.  The
!          Ritz vector associated with the Ritz value with negative
!          imaginary part is simply the complex conjugate of the Ritz vector
!          associated with the positive imaginary part.
!
!          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
!
!          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
!          the array Z may be set equal to first NEV+1 columns of the Arnoldi
!          basis array V computed by DNAUPD .  In this case the Arnoldi basis
!          will be destroyed and overwritten with the eigenvector basis.
!
!  LDZ     Integer.  (INPUT)
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
!
!  SIGMAR  Double precision   (INPUT)
!          If IPARAM(7) = 3 or 4, represents the real part of the shift.
!          Not referenced if IPARAM(7) = 1 or 2.
!
!  SIGMAI  Double precision   (INPUT)
!          If IPARAM(7) = 3 or 4, represents the imaginary part of the shift.
!          Not referenced if IPARAM(7) = 1 or 2. See remark 3 below.
!
!  WORKEV  Double precision  work array of dimension 3*NCV.  (WORKSPACE)
!
!  **** The remaining arguments MUST be the same as for the   ****
!  **** call to DNAUPD  that was just completed.               ****
!
!  NOTE: The remaining arguments
!
!           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
!           WORKD, WORKL, LWORKL, INFO
!
!         must be passed directly to DNEUPD  following the last call
!         to DNAUPD .  These arguments MUST NOT BE MODIFIED between
!         the the last call to DNAUPD  and the call to DNEUPD .
!
!  Three of these parameters (V, WORKL, INFO) are also output parameters:
!
!  V       Double precision  N by NCV array.  (INPUT/OUTPUT)
!
!          Upon INPUT: the NCV columns of V contain the Arnoldi basis
!                      vectors for OP as constructed by DNAUPD  .
!
!          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
!                       contain approximate Schur vectors that span the
!                       desired invariant subspace.  See Remark 2 below.
!
!          NOTE: If the array Z has been set equal to first NEV+1 columns
!          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
!          Arnoldi basis held by V has been overwritten by the desired
!          Ritz vectors.  If a separate array Z has been passed then
!          the first NCONV=IPARAM(5) columns of V will contain approximate
!          Schur vectors that span the desired invariant subspace.
!
!  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          WORKL(1:ncv*ncv+3*ncv) contains information obtained in
!          dnaupd .  They are not changed by dneupd .
!          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the
!          real and imaginary part of the untransformed Ritz values,
!          the upper quasi-triangular matrix for H, and the
!          associated matrix representation of the invariant subspace for H.
!
!          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
!          of the above information computed by dneupd .
!          -------------------------------------------------------------
!          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
!                     original system.
!          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
!                     the original system.
!          IPNTR(11): pointer to the NCV corresponding error bounds.
!          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
!                     Schur matrix for H.
!          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
!                     of the upper Hessenberg matrix H. Only referenced by
!                     dneupd  if RVEC = .TRUE. See Remark 2 below.
!          -------------------------------------------------------------
!
!  INFO    Integer.  (OUTPUT)
!          Error flag on output.
!
!          =  0: Normal exit.
!
!          =  1: The Schur form computed by LAPACK routine dlahqr
!                could not be reordered by LAPACK routine dtrsen .
!                Re-enter subroutine dneupd  with IPARAM(5)=NCV and
!                increase the size of the arrays DR and DI to have
!                dimension at least dimension NCV and allocate at least NCV
!                columns for Z. NOTE: Not necessary if Z and V share
!                the same space. Please notify the authors if this error
!                occurs.
!
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work WORKL array is not sufficient.
!          = -8: Error return from calculation of a real Schur form.
!                Informational error from LAPACK routine dlahqr .
!          = -9: Error return from calculation of eigenvectors.
!                Informational error from LAPACK routine dtrevc .
!          = -10: IPARAM(7) must be 1,2,3,4.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: HOWMNY = 'S' not yet implemented
!          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
!          = -14: DNAUPD  did not find any eigenvalues to sufficient
!                 accuracy.
!          = -15: DNEUPD got a different count of the number of converged
!                 Ritz values than DNAUPD got.  This indicates the user
!                 probably made an error in passing data from DNAUPD to
!                 DNEUPD or that the data was modified before entering
!                 DNEUPD
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
!     Real Matrices", Linear Algebra and its Applications, vol 88/89,
!     pp 575-595, (1987).
!
!\Routines called:
!     ivout   ARPACK utility routine that prints integers.
!     dmout    ARPACK utility routine that prints matrices
!     dvout    ARPACK utility routine that prints vectors.
!     dgeqr2   LAPACK routine that computes the QR factorization of
!             a matrix.
!     dlacpy   LAPACK matrix copy routine.
!     dlahqr   LAPACK routine to compute the real Schur form of an
!             upper Hessenberg matrix.
!     dlamch   LAPACK routine that determines machine constants.
!     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dlaset   LAPACK matrix initialization routine.
!     dorm2r   LAPACK routine that applies an orthogonal matrix in
!             factored form.
!     dtrevc   LAPACK routine to compute the eigenvectors of a matrix
!             in upper quasi-triangular form.
!     dtrsen   LAPACK routine that re-orders the Schur form.
!     dtrmm    Level 3 BLAS matrix times an upper triangular matrix.
!     dger     Level 2 BLAS rank one update to a matrix.
!     dcopy    Level 1 BLAS that copies one vector to another .
!     ddot     Level 1 BLAS that computes the scalar product of two vectors.
!     dnrm2    Level 1 BLAS that computes the norm of a vector.
!     dscal    Level 1 BLAS that scales a vector.
!
!\Remarks
!
!  1. Currently only HOWMNY = 'A' and 'P' are implemented.
!
!     Let trans(X) denote the transpose of X.
!
!  2. Schur vectors are an orthogonal representation for the basis of
!     Ritz vectors. Thus, their numerical properties are often superior.
!     If RVEC = .TRUE. then the relationship
!             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
!     trans(V(:,1:IPARAM(5))) * V(:,1:IPARAM(5)) = I are approximately
!     satisfied. Here T is the leading submatrix of order IPARAM(5) of the
!     real upper quasi-triangular matrix stored workl(ipntr(12)). That is,
!     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
!     each 2-by-2 diagonal block has its diagonal elements equal and its
!     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
!     diagonal block is a complex conjugate pair of Ritz values. The real
!     Ritz values are stored on the diagonal of T.
!
!  3. If IPARAM(7) = 3 or 4 and SIGMAI is not equal zero, then the user must
!     form the IPARAM(5) Rayleigh quotients in order to transform the Ritz
!     values computed by DNAUPD  for OP to those of A*z = lambda*B*z.
!     Set RVEC = .true. and HOWMNY = 'A', and
!     compute
!           trans(Z(:,I)) * A * Z(:,I) if DI(I) = 0.
!     If DI(I) is not equal to zero and DI(I+1) = - D(I),
!     then the desired real and imaginary parts of the Ritz value are
!           trans(Z(:,I)) * A * Z(:,I) +  trans(Z(:,I+1)) * A * Z(:,I+1),
!           trans(Z(:,I)) * A * Z(:,I+1) -  trans(Z(:,I+1)) * A * Z(:,I),
!     respectively.
!     Another possibility is to set RVEC = .true. and HOWMNY = 'P' and
!     compute trans(V(:,1:IPARAM(5))) * A * V(:,1:IPARAM(5)) and then an upper
!     quasi-triangular matrix of order IPARAM(5) is computed. See remark
!     2 above.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Chao Yang                    Houston, Texas
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: neupd.F   SID: 2.7   DATE OF SID: 09/20/00   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
subroutine dneupd (rvec , howmny, select, dr    , di, &
& z    , ldz   , sigmar, sigmai, workev, &
& bmat , n     , which , nev   , tol, &
& resid, ncv   , v     , ldv   , iparam, &
& ipntr, workd , workl , lworkl, info)
!
implicit none
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
character  bmat, howmny, which*2
logical    rvec
integer    info, ldz, ldv, lworkl, n, ncv, nev
Double precision &
& sigmar, sigmai, tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
integer    iparam(8), ipntr(14)
logical    select(ncv)
Double precision &
& dr(nev+1)    , di(nev+1), resid(n)  , &
& v(ldv,ncv)   , z(ldz,*) , workd(3*n), &
& workl(lworkl), workev(3*ncv)
!
!     %------------%
!     | Parameters |
!     %------------%
!
Double precision &
& one, zero
parameter (one = 1.0D+0 , zero = 0.0D-5 )
!
integer(4) &
& ione
parameter (ione = 1)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
character  type*6
integer    bounds, ierr  , ih    , ihbds   , &
& iheigr, iheigi, iconj , nconv   , &
& invsub, iuptri, iwev  , iwork(1), &
& j     , k     , ldh   , ldq     , &
& mode  , outncv, ritzr   , &
& ritzi , wri   , wrr   , irr     , &
& iri   , ibd   , ishift, numcnv  , &
& np    , jj
logical    reord
Double precision &
& conds  , rnorm, sep  , temp, &
& vl(1,1), temp1, eps23
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
external   dcopy  , dger   , dgeqr2 , dlacpy , &
& dlahqr , dlaset , dorm2r , &
& dtrevc , dtrmm  , dtrsen , dscal
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
Double precision &
& dlapy2 , dnrm2 , dlamch , ddot
external   dlapy2 , dnrm2 , dlamch , ddot
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
intrinsic    abs, min, sqrt
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %------------------------%
!     | Set default parameters |
!     %------------------------%
!
mode = iparam(7)
nconv = iparam(5)
info = 0
!
!     %---------------------------------%
!     | Get machine dependent constant. |
!     %---------------------------------%
!
eps23 = dlamch ('Epsilon-Machine')
eps23 = eps23**(2.0D+0  / 3.0D+0 )
!
!     %--------------%
!     | Quick return |
!     %--------------%
!
ierr = 0
!
if (nconv .le. 0) then
ierr = -14
else if (n .le. 0) then
ierr = -1
else if (nev .le. 0) then
ierr = -2
else if (ncv .le. nev+1 .or.  ncv .gt. n) then
ierr = -3
else if (which .ne. 'LM' .and. &
& which .ne. 'SM' .and. &
& which .ne. 'LR' .and. &
& which .ne. 'SR' .and. &
& which .ne. 'LI' .and. &
& which .ne. 'SI') then
ierr = -5
else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
ierr = -6
else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
ierr = -7
else if ( (howmny .ne. 'A' .and. &
& howmny .ne. 'P' .and. &
& howmny .ne. 'S') .and. rvec ) then
ierr = -13
else if (howmny .eq. 'S' ) then
ierr = -12
end if
!
if (mode .eq. 1 .or. mode .eq. 2) then
type = 'REGULR'
else if (mode .eq. 3 .and. (abs(sigmai) .le. zero)) then
type = 'SHIFTI'
else if (mode .eq. 3 ) then
type = 'REALPT'
else if (mode .eq. 4 ) then
type = 'IMAGPT'
else
ierr = -10
end if
if (mode .eq. 1 .and. bmat .eq. 'G')    ierr = -11
!
!     %------------%
!     | Error Exit |
!     %------------%
!
if (ierr .ne. 0) then
info = ierr
go to 9000
end if
!
!     %--------------------------------------------------------%
!     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q   |
!     | etc... and the remaining workspace.                    |
!     | Also update pointer to be used on output.              |
!     | Memory is laid out as follows:                         |
!     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
!     | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary   |
!     |                                   parts of ritz values |
!     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds   |
!     %--------------------------------------------------------%
!
!     %-----------------------------------------------------------%
!     | The following is used and set by DNEUPD .                  |
!     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
!     |                             real part of the Ritz values. |
!     | workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed |
!     |                        imaginary part of the Ritz values. |
!     | workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed |
!     |                           error bounds of the Ritz values |
!     | workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper |
!     |                             quasi-triangular matrix for H |
!     | workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the    |
!     |       associated matrix representation of the invariant   |
!     |       subspace for H.                                     |
!     | GRAND total of NCV * ( 3 * NCV + 6 ) locations.           |
!     %-----------------------------------------------------------%
!
ih     = ipntr(5)
ritzr  = ipntr(6)
ritzi  = ipntr(7)
bounds = ipntr(8)
ldh    = ncv
ldq    = ncv
iheigr = bounds + ldh
iheigi = iheigr + ldh
ihbds  = iheigi + ldh
iuptri = ihbds  + ldh
invsub = iuptri + ldh*ncv
ipntr(9)  = iheigr
ipntr(10) = iheigi
ipntr(11) = ihbds
ipntr(12) = iuptri
ipntr(13) = invsub
wrr = 1
wri = ncv + 1
iwev = wri + ncv
!
!     %-----------------------------------------%
!     | irr points to the REAL part of the Ritz |
!     |     values computed by _neigh before    |
!     |     exiting _naup2.                     |
!     | iri points to the IMAGINARY part of the |
!     |     Ritz values computed by _neigh      |
!     |     before exiting _naup2.              |
!     | ibd points to the Ritz estimates        |
!     |     computed by _neigh before exiting   |
!     |     _naup2.                             |
!     %-----------------------------------------%
!
irr = ipntr(14)+ncv*ncv
iri = irr+ncv
ibd = iri+ncv
!
!     %------------------------------------%
!     | RNORM is B-norm of the RESID(1:N). |
!     %------------------------------------%
!
rnorm = workl(ih+2)
workl(ih+2) = zero
!
!
if (rvec) then
!
reord = .false.
!
!        %---------------------------------------------------%
!        | Use the temporary bounds array to store indices   |
!        | These will be used to mark the select array later |
!        %---------------------------------------------------%
!
do 10 j = 1,ncv
workl(bounds+j-1) = j
select(j) = .false.
10 continue
!
!        %-------------------------------------%
!        | Select the wanted Ritz values.      |
!        | Sort the Ritz values so that the    |
!        | wanted ones appear at the tailing   |
!        | NEV positions of workl(irr) and     |
!        | workl(iri).  Move the corresponding |
!        | error estimates in workl(bound)     |
!        | accordingly.                        |
!        %-------------------------------------%
!
np     = ncv - nev
ishift = 0
call dngets (ishift       , which     , nev       , &
& np           , workl(irr), workl(iri), &
& workl(bounds))
!
!
!        %-----------------------------------------------------%
!        | Record indices of the converged wanted Ritz values  |
!        | Mark the select array for possible reordering       |
!        %-----------------------------------------------------%
!
numcnv = 0
do 11 j = 1,ncv
temp1 = max(eps23, &
& dlapy2 ( workl(irr+ncv-j), workl(iri+ncv-j) ))
jj = INT(workl(bounds + ncv - j))
if (numcnv .lt. nconv .and. &
& workl(ibd+jj-1) .le. tol*temp1) then
select(jj) = .true.
numcnv = numcnv + 1
if (jj .gt. nev) reord = .true.
endif
11 continue
!
!        %-----------------------------------------------------------%
!        | Check the count (numcnv) of converged Ritz values with    |
!        | the number (nconv) reported by dnaupd.  If these two      |
!        | are different then there has probably been an error       |
!        | caused by incorrect passing of the dnaupd data.           |
!        %-----------------------------------------------------------%
!
!
if (numcnv .ne. nconv) then
info = -15
go to 9000
end if
!
!        %-----------------------------------------------------------%
!        | Call LAPACK routine dlahqr  to compute the real Schur form |
!        | of the upper Hessenberg matrix returned by DNAUPD .        |
!        | Make a copy of the upper Hessenberg matrix.               |
!        | Initialize the Schur vector matrix Q to the identity.     |
!        %-----------------------------------------------------------%
!
call dcopy (ldh*ncv, workl(ih), ione, &
& workl(iuptri), ione)
call dlaset ('All', ncv, ncv, &
& zero , one, workl(invsub), &
& ldq)
call dlahqr (.true., .true.       , ncv, &
& ione, ncv, workl(iuptri), &
& ldh   , workl(iheigr), workl(iheigi), &
& ione, ncv, workl(invsub), &
& ldq   , ierr)
call dcopy (ncv, workl(invsub+ncv-1), ldq, &
& workl(ihbds), 1)
!
if (ierr .ne. 0) then
info = -8
go to 9000
end if
!
if (reord) then
!
!           %-----------------------------------------------------%
!           | Reorder the computed upper quasi-triangular matrix. |
!           %-----------------------------------------------------%
!
call dtrsen ('None'       , 'V'          , &
& select       , ncv          , &
& workl(iuptri), ldh          , &
& workl(invsub), ldq          , &
& workl(iheigr), workl(iheigi), &
& nconv        , conds        , &
& sep          , workl(ihbds) , &
& ncv          , iwork        , &
& 1            , ierr)
!
if (ierr .eq. 1) then
info = 1
go to 9000
end if
!
end if
!
!        %---------------------------------------%
!        | Copy the last row of the Schur vector |
!        | into workl(ihbds).  This will be used |
!        | to compute the Ritz estimates of      |
!        | converged Ritz values.                |
!        %---------------------------------------%
!
call dcopy (ncv, workl(invsub+ncv-1), ldq, &
& workl(ihbds), ione)
!
!        %----------------------------------------------------%
!        | Place the computed eigenvalues of H into DR and DI |
!        | if a spectral transformation was not used.         |
!        %----------------------------------------------------%
!
if (type .eq. 'REGULR') then
call dcopy (nconv, workl(iheigr), 1, dr(ione), 1)
call dcopy (nconv, workl(iheigi), 1, di(ione), 1)
end if
!
!        %----------------------------------------------------------%
!        | Compute the QR factorization of the matrix representing  |
!        | the wanted invariant subspace located in the first NCONV |
!        | columns of workl(invsub,ldq).                            |
!        %----------------------------------------------------------%
!
call dgeqr2 (ncv, nconv , workl(invsub), &
& ldq, workev, workev(ncv+1), &
& ierr)
!
!        %---------------------------------------------------------%
!        | * Postmultiply V by Q using dorm2r .                     |
!        | * Copy the first NCONV columns of VQ into Z.            |
!        | * Postmultiply Z by R.                                  |
!        | The N by NCONV matrix Z is now a matrix representation  |
!        | of the approximate invariant subspace associated with   |
!        | the Ritz values in workl(iheigr) and workl(iheigi)      |
!        | The first NCONV columns of V are now approximate Schur  |
!        | vectors associated with the real upper quasi-triangular |
!        | matrix of order NCONV in workl(iuptri)                  |
!        %---------------------------------------------------------%
!
call dorm2r ('Right', 'Notranspose', n            , &
& ncv   , nconv        , workl(invsub), &
& ldq   , workev       , v            , &
& ldv   , workd(n+1)   , ierr)
call dlacpy ('All', n, nconv, v(1,1), ldv, z, ldz)
!
do 20 j=1, nconv
!
!           %---------------------------------------------------%
!           | Perform both a column and row scaling if the      |
!           | diagonal element of workl(invsub,ldq) is negative |
!           | I'm lazy and don't take advantage of the upper    |
!           | quasi-triangular form of workl(iuptri,ldq)        |
!           | Note that since Q is orthogonal, R is a diagonal  |
!           | matrix consisting of plus or minus ones           |
!           %---------------------------------------------------%
!
if (workl(invsub+(j-1)*ldq+j-1) .lt. zero) then
call dscal (nconv, -one, workl(iuptri+j-1), &
& ldq)
call dscal (nconv, -one, workl(iuptri+(j-1)*ldq), &
& ione)
end if
!
20 continue
!
if (howmny .eq. 'A') then
!
!           %--------------------------------------------%
!           | Compute the NCONV wanted eigenvectors of T |
!           | located in workl(iuptri,ldq).              |
!           %--------------------------------------------%
!
do 30 j=1, ncv
if (j .le. nconv) then
select(j) = .true.
else
select(j) = .false.
end if
30 continue
!
call dtrevc ('Right', 'Select'     , select       , &
& ncv    , workl(iuptri), ldq          , &
& vl, ione, workl(invsub), &
& ldq    , ncv          , outncv       , &
& workev(1) , ierr)
!
if (ierr .ne. 0) then
info = -9
go to 9000
end if
!
!           %------------------------------------------------%
!           | Scale the returning eigenvectors so that their |
!           | Euclidean norms are all one. LAPACK subroutine |
!           | dtrevc  returns each eigenvector normalized so  |
!           | that the element of largest magnitude has      |
!           | magnitude 1;                                   |
!           %------------------------------------------------%
!
iconj = 0
do 40 j=1, nconv
!
if ( abs( workl(iheigi+j-1) ) .le. zero ) then
!
!                 %----------------------%
!                 | real eigenvalue case |
!                 %----------------------%
!
temp = dnrm2 ( ncv, workl(invsub+(j-1)*ldq), &
& ione)
call dscal ( ncv, one / temp, &
& workl(invsub+(j-1)*ldq), ione )
!
else
!
!                 %-------------------------------------------%
!                 | Complex conjugate pair case. Note that    |
!                 | since the real and imaginary part of      |
!                 | the eigenvector are stored in consecutive |
!                 | columns, we further normalize by the      |
!                 | square root of two.                       |
!                 %-------------------------------------------%
!
if (iconj .eq. 0) then
temp = dlapy2 (dnrm2 (ncv, &
& workl(invsub+(j-1)*ldq), &
& 1), &
& dnrm2 (ncv, &
& workl(invsub+j*ldq), &
& 1))
call dscal (ncv, one/temp, &
& workl(invsub+(j-1)*ldq), &
& ione)
call dscal (ncv, one/temp, &
& workl(invsub+j*ldq), &
& ione )
iconj = 1
else
iconj = 0
end if
!
end if
!
40 continue
!
call dgemv ('T', ncv, nconv, one, workl(invsub), &
& ldq, workl(ihbds), ione, &
& zero, workev(ione), 1)
!
iconj = 0
do 45 j=1, nconv
if (workl(iheigi+j-1) .le. zero) then
!
!                 %-------------------------------------------%
!                 | Complex conjugate pair case. Note that    |
!                 | since the real and imaginary part of      |
!                 | the eigenvector are stored in consecutive |
!                 %-------------------------------------------%
!
if (iconj .eq. 0) then
workev(j) = dlapy2 (workev(j), workev(j+1))
workev(j+1) = workev(j)
iconj = 1
else
iconj = 0
end if
end if
45 continue
!
!           %---------------------------------------%
!           | Copy Ritz estimates into workl(ihbds) |
!           %---------------------------------------%
!
call dcopy (nconv, workev, 1, workl(ihbds), 1)
!
!           %---------------------------------------------------------%
!           | Compute the QR factorization of the eigenvector matrix  |
!           | associated with leading portion of T in the first NCONV |
!           | columns of workl(invsub,ldq).                           |
!           %---------------------------------------------------------%
!
call dgeqr2 (ncv, nconv , workl(invsub), &
& ldq, workev, workev(ncv+1), &
& ierr)
!
!           %----------------------------------------------%
!           | * Postmultiply Z by Q.                       |
!           | * Postmultiply Z by R.                       |
!           | The N by NCONV matrix Z is now contains the  |
!           | Ritz vectors associated with the Ritz values |
!           | in workl(iheigr) and workl(iheigi).          |
!           %----------------------------------------------%
!
call dorm2r ('Right', 'Notranspose', n            , &
& ncv  , nconv        , workl(invsub), &
& ldq  , workev       , z            , &
& ldz  , workd(n+1)   , ierr)
!
call dtrmm ('Right'   , 'Upper'       , 'No transpose', &
& 'Non-unit', n            , nconv         , &
& one       , workl(invsub), ldq           , &
& z         , ldz)
!
end if
!
else
!
!        %------------------------------------------------------%
!        | An approximate invariant subspace is not needed.     |
!        | Place the Ritz values computed DNAUPD  into DR and DI |
!        %------------------------------------------------------%
!
call dcopy (nconv, workl(ritzr), 1, dr(ione), 1)
call dcopy (nconv, workl(ritzi), 1, di(ione), 1)
call dcopy (nconv, workl(ritzr), 1, workl(iheigr), 1)
call dcopy (nconv, workl(ritzi), 1, workl(iheigi), 1)
call dcopy (nconv, workl(bounds), 1, workl(ihbds), 1)
end if
!
!     %------------------------------------------------%
!     | Transform the Ritz values and possibly vectors |
!     | and corresponding error bounds of OP to those  |
!     | of A*x = lambda*B*x.                           |
!     %------------------------------------------------%
!
if (type .eq. 'REGULR') then
!
if (rvec) &
& call dscal (ncv, rnorm, workl(ihbds), ione)
!
else
!
!        %---------------------------------------%
!        |   A spectral transformation was used. |
!        | * Determine the Ritz estimates of the |
!        |   Ritz values in the original system. |
!        %---------------------------------------%
!
if (type .eq. 'SHIFTI') then
!
if (rvec) &
& call dscal (ncv, rnorm, workl(ihbds), ione)
!
do 50 k=1, ncv
temp = dlapy2 ( workl(iheigr+k-1), &
& workl(iheigi+k-1) )
workl(ihbds+k-1) = abs( workl(ihbds+k-1) ) &
& / temp / temp
50 continue
!
else if (type .eq. 'REALPT') then
!
do 60 k=1, ncv
60 continue
!
else if (type .eq. 'IMAGPT') then
!
do 70 k=1, ncv
70 continue
!
end if
!
!        %-----------------------------------------------------------%
!        | *  Transform the Ritz values back to the original system. |
!        |    For TYPE = 'SHIFTI' the transformation is              |
!        |             lambda = 1/theta + sigma                      |
!        |    For TYPE = 'REALPT' or 'IMAGPT' the user must from     |
!        |    Rayleigh quotients or a projection. See remark 3 above.|
!        | NOTES:                                                    |
!        | *The Ritz vectors are not affected by the transformation. |
!        %-----------------------------------------------------------%
!
if (type .eq. 'SHIFTI') then
!
do 80 k=1, ncv
temp = dlapy2 ( workl(iheigr+k-1), &
& workl(iheigi+k-1) )
workl(iheigr+k-1) = workl(iheigr+k-1)/temp/temp &
& + sigmar
workl(iheigi+k-1) = -workl(iheigi+k-1)/temp/temp &
& + sigmai
80 continue
!
call dcopy (nconv, workl(iheigr), 1, dr(ione), 1)
call dcopy (nconv, workl(iheigi), 1, di(ione), 1)
!
else if (type .eq. 'REALPT' .or. type .eq. 'IMAGPT') then
!
call dcopy (nconv, workl(iheigr), 1, dr(ione), 1)
call dcopy (nconv, workl(iheigi), 1, di(ione), 1)
!
end if
!
end if
!
!     %-------------------------------------------------%
!     | Eigenvector Purification step. Formally perform |
!     | one of inverse subspace iteration. Only used    |
!     | for MODE = 2.                                   |
!     %-------------------------------------------------%
!
if (rvec .and. howmny .eq. 'A' .and. type .eq. 'SHIFTI') then
!
!        %------------------------------------------------%
!        | Purify the computed Ritz vectors by adding a   |
!        | little bit of the residual vector:             |
!        |                      T                         |
!        |          resid(:)*( e    s ) / theta           |
!        |                      NCV                       |
!        | where H s = s theta. Remember that when theta  |
!        | has nonzero imaginary part, the corresponding  |
!        | Ritz vector is stored across two columns of Z. |
!        %------------------------------------------------%
!
iconj = 0
do 110 j=1, nconv
!            if (workl(iheigi+j-1) .eq. zero) then
if (abs(workl(iheigi+j-1)) .le. zero) then
workev(j) =  workl(invsub+(j-1)*ldq+ncv-1) / &
& workl(iheigr+j-1)
else if (iconj .eq. 0) then
temp = dlapy2 ( workl(iheigr+j-1), workl(iheigi+j-1) )
workev(j) = ( workl(invsub+(j-1)*ldq+ncv-1) * &
& workl(iheigr+j-1) + &
& workl(invsub+j*ldq+ncv-1) * &
& workl(iheigi+j-1) ) / temp / temp
workev(j+1) = ( workl(invsub+j*ldq+ncv-1) * &
& workl(iheigr+j-1) - &
& workl(invsub+(j-1)*ldq+ncv-1) * &
& workl(iheigi+j-1) ) / temp / temp
iconj = 1
else
iconj = 0
end if
110 continue
!
!        %---------------------------------------%
!        | Perform a rank one update to Z and    |
!        | purify all the Ritz vectors together. |
!        %---------------------------------------%
!
call dger (n, nconv, one, resid, 1, workev, 1, z, ldz)
!
end if
!
9000 continue
!
return
!
!     %---------------%
!     | End of DNEUPD  |
!     %---------------%
!
end
!\BeginDoc
!
!\Name: dnaupd
!
!\Description:
!  Reverse communication interface for the Implicitly Restarted Arnoldi
!  iteration. This subroutine computes approximations to a few eigenpairs
!  of a linear operator "OP" with respect to a semi-inner product defined by
!  a symmetric positive semi-definite real matrix B. B may be the identity
!  matrix. NOTE: If the linear operator "OP" is real and symmetric
!  with respect to the real positive semi-definite symmetric matrix B,
!  i.e. B*OP = (OP`)*B, then subroutine dsaupd should be used instead.
!
!  The computed approximate eigenvalues are called Ritz values and
!  the corresponding approximate eigenvectors are called Ritz vectors.
!
!  dnaupd is usually called iteratively to solve one of the
!  following problems:
!
!  Mode 1:  A*x = lambda*x.
!           ===> OP = A  and  B = I.
!
!  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
!           ===> OP = inv[M]*A  and  B = M.
!           ===> (If M can be factored see remark 3 below)
!
!  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
!           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M.
!           ===> shift-and-invert mode (in real arithmetic)
!           If OP*x = amu*x, then
!           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
!           Note: If sigma is real, i.e. imaginary part of sigma is zero;
!                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M
!                 amu == 1/(lambda-sigma).
!
!  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
!           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M.
!           ===> shift-and-invert mode (in real arithmetic)
!           If OP*x = amu*x, then
!           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
!
!  Both mode 3 and 4 give the same enhancement to eigenvalues close to
!  the (complex) shift sigma.  However, as lambda goes to infinity,
!  the operator OP in mode 4 dampens the eigenvalues more strongly than
!  does OP defined in mode 3.
!
!  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
!        should be accomplished either by a direct method
!        using a sparse matrix factorization and solving
!
!           [A - sigma*M]*w = v  or M*w = v,
!
!        or through an iterative method for solving these
!        systems.  If an iterative method is used, the
!        convergence test must be more stringent than
!        the accuracy requirements for the eigenvalue
!        approximations.
!
!\Usage:
!  call dnaupd
!     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
!       IPNTR, WORKD, WORKL, LWORKL, INFO )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first
!          call to dnaupd.  IDO will be set internally to
!          indicate the type of operation to be performed.  Control is
!          then given back to the calling routine which has the
!          responsibility to carry out the requested operation and call
!          dnaupd with the result.  The operand is given in
!          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    In mode 3 and 4, the vector B * X is already
!                    available in WORKD(ipntr(3)).  It does not
!                    need to be recomputed in forming OP * X.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO =  3: compute the IPARAM(8) real and imaginary parts
!                    of the shifts where INPTR(14) is the pointer
!                    into WORKL for placing the shifts. See Remark
!                    5 below.
!          IDO = 99: done
!          -------------------------------------------------------------
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.
!          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
!          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  WHICH   Character*2.  (INPUT)
!          'LM' -> want the NEV eigenvalues of largest magnitude.
!          'SM' -> want the NEV eigenvalues of smallest magnitude.
!          'LR' -> want the NEV eigenvalues of largest real part.
!          'SR' -> want the NEV eigenvalues of smallest real part.
!          'LI' -> want the NEV eigenvalues of largest imaginary part.
!          'SI' -> want the NEV eigenvalues of smallest imaginary part.
!
!  NEV     Integer.  (INPUT/OUTPUT)
!          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
!
!  TOL     Double precision scalar.  (INPUT)
!          Stopping criterion: the relative accuracy of the Ritz value
!          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
!          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
!          DEFAULT = DLAMCH('EPS')  (machine precision as computed
!                    by the LAPACK auxiliary subroutine DLAMCH).
!
!  RESID   Double precision array of length N.  (INPUT/OUTPUT)
!          On INPUT:
!          If INFO .EQ. 0, a random initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          On OUTPUT:
!          RESID contains the final residual vector.
!
!  NCV     Integer.  (INPUT)
!          Number of columns of the matrix V. NCV must satisfy the two
!          inequalities 2 <= NCV-NEV and NCV <= N.
!          This will indicate how many Arnoldi vectors are generated
!          at each iteration.  After the startup phase in which NEV
!          Arnoldi vectors are generated, the algorithm generates
!          approximately NCV-NEV Arnoldi vectors at each subsequent update
!          iteration. Most of the cost in generating each Arnoldi vector is
!          in the matrix-vector operation OP*x.
!          NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of Ritz
!          values are kept together. (See remark 4 below)
!
!  V       Double precision array N by NCV.  (OUTPUT)
!          Contains the final set of Arnoldi basis vectors.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling program.
!
!  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
!          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
!          The shifts selected at each iteration are used to restart
!          the Arnoldi iteration in an implicit fashion.
!          -------------------------------------------------------------
!          ISHIFT = 0: the shifts are provided by the user via
!                      reverse communication.  The real and imaginary
!                      parts of the NCV eigenvalues of the Hessenberg
!                      matrix H are returned in the part of the WORKL
!                      array corresponding to RITZR and RITZI. See remark
!                      5 below.
!          ISHIFT = 1: exact shifts with respect to the current
!                      Hessenberg matrix H.  This is equivalent to
!                      restarting the iteration with a starting vector
!                      that is a linear combination of approximate Schur
!                      vectors associated with the "wanted" Ritz values.
!          -------------------------------------------------------------
!
!          IPARAM(2) = No longer referenced.
!          IPARAM(3) = MXITER
!          On INPUT:  maximum number of Arnoldi update iterations allowed.
!          On OUTPUT: actual number of Arnoldi update iterations taken.
!
!          IPARAM(4) = NB: blocksize to be used in the recurrence.
!          The code currently works only for NB = 1.
!
!          IPARAM(5) = NCONV: number of "converged" Ritz values.
!          This represents the number of Ritz values that satisfy
!          the convergence criterion.
!
!          IPARAM(6) = IUPD
!          No longer referenced. Implicit restarting is ALWAYS used.
!
!          IPARAM(7) = MODE
!          On INPUT determines what type of eigenproblem is being solved.
!          Must be 1,2,3,4; See under \Description of dnaupd for the
!          four modes available.
!
!          IPARAM(8) = NP
!          When ido = 3 and the user provides shifts through reverse
!          communication (IPARAM(1)=0), dnaupd returns NP, the number
!          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
!          5 below.
!
!rc          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
!rc          OUTPUT: NUMOP  = total number of OP*x operations,
!rc                  NUMOPB = total number of B*x operations if BMAT='G',
!rc                  NUMREO = total number of steps of re-orthogonalization.
!
!  IPNTR   Integer array of length 14.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD and WORKL
!          arrays for matrices/vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X in WORKD.
!          IPNTR(2): pointer to the current result vector Y in WORKD.
!          IPNTR(3): pointer to the vector B * X in WORKD when used in
!                    the shift-and-invert mode.
!          IPNTR(4): pointer to the next available location in WORKL
!                    that is untouched by the program.
!          IPNTR(5): pointer to the NCV by NCV upper Hessenberg matrix
!                    H in WORKL.
!          IPNTR(6): pointer to the real part of the ritz value array
!                    RITZR in WORKL.
!          IPNTR(7): pointer to the imaginary part of the ritz value array
!                    RITZI in WORKL.
!          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
!                    with the Ritz values located in RITZR and RITZI in WORKL.
!
!          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
!
!          Note: IPNTR(9:13) is only referenced by dneupd. See Remark 2 below.
!
!          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
!                     original system.
!          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
!                     the original system.
!          IPNTR(11): pointer to the NCV corresponding error bounds.
!          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
!                     Schur matrix for H.
!          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
!                     of the upper Hessenberg matrix H. Only referenced by
!                     dneupd if RVEC = .TRUE. See Remark 2 below.
!          -------------------------------------------------------------
!
!  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The user should not use WORKD
!          as temporary workspace during the iteration. Upon termination
!          WORKD(1:N) contains B*RESID(1:N). If an invariant subspace
!          associated with the converged Ritz values is desired, see remark
!          2 below, subroutine dneupd uses this output.
!          See Data Distribution Note below.
!
!  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  See Data Distribution Note below.
!
!  LWORKL  Integer.  (INPUT)
!          LWORKL must be at least 3*NCV**2 + 6*NCV.
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =  0: Normal exit.
!          =  1: Maximum number of iterations taken.
!                All possible eigenvalues of OP has been found. IPARAM(5)
!                returns the number of wanted converged Ritz values.
!          =  2: No longer an informational error. Deprecated starting
!                with release 2 of ARPACK.
!          =  3: No shifts could be applied during a cycle of the
!                Implicitly restarted Arnoldi iteration. One possibility
!                is to increase the size of NCV relative to NEV.
!                See remark 4 below.
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -4: The maximum number of Arnoldi update iteration
!                must be greater than zero.
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work array is not sufficient.
!          = -8: Error return from LAPACK eigenvalue calculation;
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
!          = -12: IPARAM(1) must be equal to 0 or 1.
!          = -9999: Could not build an Arnoldi factorization.
!                   IPARAM(5) returns the size of the current Arnoldi
!                   factorization.
!
!\Remarks
!  1. The computed Ritz values are approximate eigenvalues of OP. The
!     selection of WHICH should be made with this in mind when
!     Mode = 3 and 4.  After convergence, approximate eigenvalues of the
!     original problem may be obtained with the ARPACK subroutine dneupd.
!
!  2. If a basis for the invariant subspace corresponding to the converged Ritz
!     values is needed, the user must call dneupd immediately following
!     completion of dnaupd. This is new starting with release 2 of ARPACK.
!
!  3. If M can be factored into a Cholesky factorization M = LL`
!     then Mode = 2 should not be selected.  Instead one should use
!     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular
!     linear systems should be solved with L and L` rather
!     than computing inverses.  After convergence, an approximate
!     eigenvector z of the original problem is recovered by solving
!     L`z = x  where x is a Ritz vector of OP.
!
!  4. At present there is no a-priori analysis to guide the selection
!     of NCV relative to NEV.  The only formal requrement is that NCV > NEV + 2.
!     However, it is recommended that NCV .ge. 2*NEV+1.  If many problems of
!     the same type are to be solved, one should experiment with increasing
!     NCV while keeping NEV fixed for a given test problem.  This will
!     usually decrease the required number of OP*x operations but it
!     also increases the work and storage required to maintain the orthogonal
!     basis vectors.  The optimal "cross-over" with respect to CPU time
!     is problem dependent and must be determined empirically.
!     See Chapter 8 of Reference 2 for further information.
!
!  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
!     NP = IPARAM(8) real and imaginary parts of the shifts in locations
!         real part                  imaginary part
!         -----------------------    --------------
!     1   WORKL(IPNTR(14))           WORKL(IPNTR(14)+NP)
!     2   WORKL(IPNTR(14)+1)         WORKL(IPNTR(14)+NP+1)
!                        .                          .
!                        .                          .
!                        .                          .
!     NP  WORKL(IPNTR(14)+NP-1)      WORKL(IPNTR(14)+2*NP-1).
!
!     Only complex conjugate pairs of shifts may be applied and the pairs
!     must be placed in consecutive locations. The real part of the
!     eigenvalues of the current upper Hessenberg matrix are located in
!     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the imaginary part
!     in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). They are ordered
!     according to the order defined by WHICH. The complex conjugate
!     pairs are kept together and the associated Ritz estimates are located in
!     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
!
!-----------------------------------------------------------------------
!
!\Data Distribution Note:
!
!  Fortran-D syntax:
!  ================
!  Double precision resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
!  decompose  d1(n), d2(n,ncv)
!  align      resid(i) with d1(i)
!  align      v(i,j)   with d2(i,j)
!  align      workd(i) with d1(i)     range (1:n)
!  align      workd(i) with d1(i-n)   range (n+1:2*n)
!  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
!  distribute d1(block), d2(block,:)
!  replicated workl(lworkl)
!
!  Cray MPP syntax:
!  ===============
!  Double precision  resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
!  shared     resid(block), v(block,:), workd(block,:)
!  replicated workl(lworkl)
!
!  CM2/CM5 syntax:
!  ==============
!
!-----------------------------------------------------------------------
!
!     include   'ex-nonsym.doc'
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
!  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
!     Real Matrices", Linear Algebra and its Applications, vol 88/89,
!     pp 575-595, (1987).
!
!\Routines called:
!     dnaup2  ARPACK routine that implements the Implicitly Restarted
!             Arnoldi Iteration.
!     ivout   ARPACK utility routine that prints integers.
!     second  ARPACK utility routine for timing.
!     dvout   ARPACK utility routine that prints vectors.
!     dlamch  LAPACK routine that determines machine constants.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     12/16/93: Version '1.1'
!
!\SCCS Information: @(#)
! FILE: naupd.F   SID: 2.10   DATE OF SID: 08/23/02   RELEASE: 2
!
!\Remarks
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dnaupd &
& ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, &
& ipntr, workd, workl, lworkl, info )
!
implicit none
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
character  bmat*1, which*2
integer    ido, info, ldv, lworkl, n, ncv, nev
Double precision &
& tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
integer    iparam(8), ipntr(14)
Double precision &
& resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
!
!     %------------%
!     | Parameters |
!     %------------%
!
Double precision &
& zero
parameter (zero = 0.0D+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
integer    bounds, ierr, ih, iq, ishift, iupd, iw, &
& ldh, ldq, levec, mode, mxiter, nb, &
& nev0, next, np, ritzi, ritzr, j
save       bounds, ih, iq, ishift, iupd, iw, ldh, ldq, &
& levec, mode,  mxiter, nb, nev0, next, &
& np, ritzi, ritzr
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
external   dnaup2
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
Double precision &
& dlamch
external   dlamch
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
if (levec .gt. 10000000) then
goto 9000
end if
!
if (ido .eq. 0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
!        %----------------%
!        | Error checking |
!        %----------------%
!
ierr   = 0
ishift = iparam(1)
!         levec  = iparam(2)
mxiter = iparam(3)
!         nb     = iparam(4)
nb     = 1
!
!        %--------------------------------------------%
!        | Revision 2 performs only implicit restart. |
!        %--------------------------------------------%
!
iupd   = 1
mode   = iparam(7)
!
if (n .le. 0) then
ierr = -1
else if (nev .le. 0) then
ierr = -2
else if (ncv .le. nev+1 .or.  ncv .gt. n) then
ierr = -3
else if (mxiter .le.          0) then
ierr = 4
else if (which .ne. 'LM' .and. &
& which .ne. 'SM' .and. &
& which .ne. 'LR' .and. &
& which .ne. 'SR' .and. &
& which .ne. 'LI' .and. &
& which .ne. 'SI') then
ierr = -5
else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
ierr = -6
else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
ierr = -7
else if (mode .lt. 1 .or. mode .gt. 4) then
ierr = -10
else if (mode .eq. 1 .and. bmat .eq. 'G') then
ierr = -11
else if (ishift .lt. 0 .or. ishift .gt. 1) then
ierr = -12
end if
!
!        %------------%
!        | Error Exit |
!        %------------%
!
if (ierr .ne. 0) then
info = ierr
ido  = 99
go to 9000
end if
!
!        %------------------------%
!        | Set default parameters |
!        %------------------------%
!
if (nb .le. 0)           nb = 1
if (tol .le. zero)       tol = dlamch('EpsMach')
!
!        %----------------------------------------------%
!        | NP is the number of additional steps to      |
!        | extend the length NEV Lanczos factorization. |
!        | NEV0 is the local variable designating the   |
!        | size of the invariant subspace desired.      |
!        %----------------------------------------------%
!
np     = ncv - nev
nev0   = nev
!
!        %-----------------------------%
!        | Zero out internal workspace |
!        %-----------------------------%
!
do 10 j = 1, 3*ncv**2 + 6*ncv
workl(j) = zero
10 continue
!
!        %-------------------------------------------------------------%
!        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        |
!        | etc... and the remaining workspace.                         |
!        | Also update pointer to be used on output.                   |
!        | Memory is laid out as follows:                              |
!        | workl(1:ncv*ncv) := generated Hessenberg matrix             |
!        | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary        |
!        |                                   parts of ritz values      |
!        | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds        |
!        | workl(ncv*ncv+3*ncv+1:2*ncv*ncv+3*ncv) := rotation matrix Q |
!        | workl(2*ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) := workspace       |
!        | The final workspace is needed by subroutine dneigh called   |
!        | by dnaup2. Subroutine dneigh calls LAPACK routines for      |
!        | calculating eigenvalues and the last row of the eigenvector |
!        | matrix.                                                     |
!        %-------------------------------------------------------------%
!
ldh    = ncv
ldq    = ncv
ih     = 1
ritzr  = ih     + ldh*ncv
ritzi  = ritzr  + ncv
bounds = ritzi  + ncv
iq     = bounds + ncv
iw     = iq     + ldq*ncv
next   = iw     + ncv**2 + 3*ncv
!
ipntr(4) = next
ipntr(5) = ih
ipntr(6) = ritzr
ipntr(7) = ritzi
ipntr(8) = bounds
ipntr(14) = iw
!
end if
!
!     %-------------------------------------------------------%
!     | Carry out the Implicitly restarted Arnoldi Iteration. |
!     %-------------------------------------------------------%
!
call dnaup2 &
& ( ido, bmat, n, which, nev0, np, tol, resid, mode, &
& ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritzr), &
& workl(ritzi), workl(bounds), workl(iq), ldq, workl(iw), &
& ipntr, workd, info )
!
!     %--------------------------------------------------%
!     | ido .ne. 99 implies use of reverse communication |
!     | to compute operations involving OP or shifts.    |
!     %--------------------------------------------------%
!
if (ido .eq. 3) iparam(8) = np
if (ido .ne. 99) go to 9000
!
iparam(3) = mxiter
iparam(5) = np
!
!     %------------------------------------%
!     | Exit if there was an informational |
!     | error within dnaup2.               |
!     %------------------------------------%
!
if (info .lt. 0) go to 9000
if (info .eq. 2) info = 3
!
!        %--------------------------------------------------------%
!        | Version Number & Version Date are defined in version.h |
!        %--------------------------------------------------------%
!
9000 continue
!
return
!
!     %---------------%
!     | End of dnaupd |
!     %---------------%
!
end
!\BeginDoc
!
!\Name: dnaup2
!
!\Description:
!  Intermediate level interface called by dnaupd.
!
!\Usage:
!  call dnaup2
!     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
!       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS,
!       Q, LDQ, WORKL, IPNTR, WORKD, INFO )
!
!\Arguments
!
!  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dnaupd.
!  MODE, ISHIFT, MXITER: see the definition of IPARAM in dnaupd.
!
!  NP      Integer.  (INPUT/OUTPUT)
!          Contains the number of implicit shifts to apply during
!          each Arnoldi iteration.
!          If ISHIFT=1, NP is adjusted dynamically at each iteration
!          to accelerate convergence and prevent stagnation.
!          This is also roughly equal to the number of matrix-vector
!          products (involving the operator OP) per Arnoldi iteration.
!          The logic for adjusting is contained within the current
!          subroutine.
!          If ISHIFT=0, NP is the number of shifts the user needs
!          to provide via reverse comunication. 0 < NP < NCV-NEV.
!          NP may be less than NCV-NEV for two reasons. The first, is
!          to keep complex conjugate pairs of "wanted" Ritz values
!          together. The second, is that a leading block of the current
!          upper Hessenberg matrix has split off and contains "unwanted"
!          Ritz values.
!          Upon termination of the IRA iteration, NP contains the number
!          of "converged" wanted Ritz values.
!
!  IUPD    Integer.  (INPUT)
!          IUPD .EQ. 0: use explicit restart instead implicit update.
!          IUPD .NE. 0: use implicit update.
!
!  V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT)
!          The Arnoldi basis vectors are returned in the first NEV
!          columns of V.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Double precision (NEV+NP) by (NEV+NP) array.  (OUTPUT)
!          H is used to store the generated upper Hessenberg matrix
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RITZR,  Double precision arrays of length NEV+NP.  (OUTPUT)
!  RITZI   RITZR(1:NEV) (resp. RITZI(1:NEV)) contains the real (resp.
!          imaginary) part of the computed Ritz values of OP.
!
!  BOUNDS  Double precision array of length NEV+NP.  (OUTPUT)
!          BOUNDS(1:NEV) contain the error bounds corresponding to
!          the computed Ritz values.
!
!  Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
!          Private (replicated) work array used to accumulate the
!          rotation in the shift application step.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Double precision work array of length at least
!          (NEV+NP)**2 + 3*(NEV+NP).  (INPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  It is used in shifts calculation, shifts
!          application and convergence checking.
!
!          On exit, the last 3*(NEV+NP) locations of WORKL contain
!          the Ritz values (real,imaginary) and associated Ritz
!          estimates of the current Hessenberg matrix.  They are
!          listed in the same order as returned from dneigh.
!
!          If ISHIFT .EQ. O and IDO .EQ. 3, the first 2*NP locations
!          of WORKL are used in reverse communication to hold the user
!          supplied shifts.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD for
!          vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in the
!                    shift-and-invert mode.  X is the current operand.
!          -------------------------------------------------------------
!
!  WORKD   Double precision work array of length 3*N.  (WORKSPACE)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The user should not use WORKD
!          as temporary workspace during the iteration !!!!!!!!!!
!          See Data Distribution Note in DNAUPD.
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =     0: Normal return.
!          =     1: Maximum number of iterations taken.
!                   All possible eigenvalues of OP has been found.
!                   NP returns the number of converged Ritz values.
!          =     2: No shifts could be applied.
!          =    -8: Error return from LAPACK eigenvalue calculation;
!                   This should never happen.
!          =    -9: Starting vector is zero.
!          = -9999: Could not build an Arnoldi factorization.
!                   Size that was built in returned in NP.
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
!     dgetv0  ARPACK initial vector generation routine.
!     dnaitr  ARPACK Arnoldi factorization routine.
!     dnapps  ARPACK application of implicit shifts routine.
!     dnconv  ARPACK convergence of Ritz values routine.
!     dneigh  ARPACK compute Ritz values and error bounds routine.
!     dngets  ARPACK reorder Ritz values and error bounds routine.
!     dsortc  ARPACK sorting routine.
!     ivout   ARPACK utility routine that prints integers.
!     second  ARPACK utility routine for timing.
!     dmout   ARPACK utility routine that prints matrices
!     dvout   ARPACK utility routine that prints vectors.
!     dlamch  LAPACK routine that determines machine constants.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dcopy   Level 1 BLAS that copies one vector to another .
!     ddot    Level 1 BLAS that computes the scalar product of two vectors.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     dswap   Level 1 BLAS that swaps two vectors.
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
! FILE: naup2.F   SID: 2.8   DATE OF SID: 10/17/00   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dnaup2 &
& ( ido, bmat, n, which, nev, np, tol, resid, mode, &
& ishift, mxiter, v, ldv, h, ldh, ritzr, ritzi, bounds, &
& q, ldq, workl, ipntr, workd, info )
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
implicit none
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
character  bmat*1, which*2
integer    ido, info, ishift, mode, ldh, ldq, ldv, mxiter, &
& n, nev, np
Double precision &
& tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
integer    ipntr(13)
Double precision &
& bounds(nev+np), h(ldh,nev+np), q(ldq,nev+np), resid(n), &
& ritzi(nev+np), ritzr(nev+np), v(ldv,nev+np), &
& workd(3*n), workl( (nev+np)*(nev+np+3) )
!
!     %------------%
!     | Parameters |
!     %------------%
!
Double precision &
& zero
parameter (zero = 0.0D+0)
integer    dnazero
integer(4) ione
parameter (dnazero = 0, ione = 1)

!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
character  wprime*2
logical    cnorm , getv0, initv, update, ushift, dnatrue
integer    ierr  , iter , j    , kplusp, nconv, &
& nevbef, nev0 , np0  , nptemp, numcnv
Double precision &
& rnorm , temp , eps23
save       cnorm , getv0, initv, update, ushift, &
& rnorm , iter , eps23, kplusp,  nconv , &
& nevbef, nev0 , np0  , numcnv
!
!     %-----------------------%
!     | Local array arguments |
!     %-----------------------%
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
external   dcopy , dgetv0, dnaitr, dnconv, dneigh, &
& dngets, dnapps
!     &           dngets, dnapps, second
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
!
Double precision &
& ddot, dnrm2, dlapy2, dlamch
external   ddot, dnrm2, dlapy2, dlamch
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
intrinsic    min, max, abs, sqrt
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
dnatrue = .true.
!
if (ido .eq. 0) then
!
!        %-------------------------------------%
!        | Get the machine dependent constant. |
!        %-------------------------------------%
!
eps23 = dlamch('Epsilon-Machine')
eps23 = eps23**(2.0D+0 / 3.0D+0)
!
nev0   = nev
np0    = np
!
!        %-------------------------------------%
!        | kplusp is the bound on the largest  |
!        |        Lanczos factorization built. |
!        | nconv is the current number of      |
!        |        "converged" eigenvlues.      |
!        | iter is the counter on the current  |
!        |      iteration step.                |
!        %-------------------------------------%
!
kplusp = nev + np
nconv  = 0
iter   = 0
!
!        %---------------------------------------%
!        | Set flags for computing the first NEV |
!        | steps of the Arnoldi factorization.   |
!        %---------------------------------------%
!
getv0    = .true.
update   = .false.
ushift   = .false.
cnorm    = .false.
!
if (info .ne. 0) then
!
!           %--------------------------------------------%
!           | User provides the initial residual vector. |
!           %--------------------------------------------%
!
initv = .true.
info  = 0
else
initv = .false.
end if
end if
!
!     %---------------------------------------------%
!     | Get a possibly random starting vector and   |
!     | force it into the range of the operator OP. |
!     %---------------------------------------------%
!
!   10 continue
!
if (getv0) then
call dgetv0 (ido, bmat, initv, n, ione, v, ldv, resid, rnorm, &
& ipntr, workd, info)
!
if (ido .ne. 99) go to 9000
!
if (abs(rnorm) .le. zero) then
!         if (rnorm .eq. zero) then
!
!           %-----------------------------------------%
!           | The initial vector is zero. Error exit. |
!           %-----------------------------------------%
!
info = -9
go to 1100
end if
getv0 = .false.
ido  = 0
end if
!
!     %-----------------------------------%
!     | Back from reverse communication : |
!     | continue with update step         |
!     %-----------------------------------%
!
if (update) go to 20
!
!     %-------------------------------------------%
!     | Back from computing user specified shifts |
!     %-------------------------------------------%
!
if (ushift) go to 50
!
!     %-------------------------------------%
!     | Back from computing residual norm   |
!     | at the end of the current iteration |
!     %-------------------------------------%
!
if (cnorm)  go to 100
!
!     %----------------------------------------------------------%
!     | Compute the first NEV steps of the Arnoldi factorization |
!     %----------------------------------------------------------%
!
call dnaitr (ido, bmat, n, dnazero, nev, mode, resid, rnorm, &
& v, ldv, h, ldh, ipntr, workd, info)
!
!     %---------------------------------------------------%
!     | ido .ne. 99 implies use of reverse communication  |
!     | to compute operations involving OP and possibly B |
!     %---------------------------------------------------%
!
if (ido .ne. 99) go to 9000
!
if (info .gt. 0) then
np   = info
mxiter = iter
info = -9999
go to 1200
end if
!
!     %--------------------------------------------------------------%
!     |                                                              |
!     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
!     |           Each iteration implicitly restarts the Arnoldi     |
!     |           factorization in place.                            |
!     |                                                              |
!     %--------------------------------------------------------------%
!
1000 continue
!
iter = iter + 1
!
!        %-----------------------------------------------------------%
!        | Compute NP additional steps of the Arnoldi factorization. |
!        | Adjust NP since NEV might have been updated by last call  |
!        | to the shift application routine dnapps.                  |
!        %-----------------------------------------------------------%
!
np  = kplusp - nev
!
!        %-----------------------------------------------------------%
!        | Compute NP additional steps of the Arnoldi factorization. |
!        %-----------------------------------------------------------%
!
ido = 0
20 continue
update = .true.
!
call dnaitr (ido  , bmat, n  , nev, np , mode , resid, &
& rnorm, v   , ldv, h  , ldh, ipntr, workd, &
& info)
!
!        %---------------------------------------------------%
!        | ido .ne. 99 implies use of reverse communication  |
!        | to compute operations involving OP and possibly B |
!        %---------------------------------------------------%
!
if (ido .ne. 99) go to 9000
!
if (info .gt. 0) then
np = info
mxiter = iter
info = -9999
go to 1200
end if
update = .false.
!
!        %--------------------------------------------------------%
!        | Compute the eigenvalues and corresponding error bounds |
!        | of the current upper Hessenberg matrix.                |
!        %--------------------------------------------------------%
!
call dneigh (rnorm, kplusp, h, ldh, ritzr, ritzi, bounds, &
& q, ldq, workl, ierr)
!
if (ierr .ne. 0) then
info = -8
go to 1200
end if
!
!        %----------------------------------------------------%
!        | Make a copy of eigenvalues and corresponding error |
!        | bounds obtained from dneigh.                       |
!        %----------------------------------------------------%
!
call dcopy(kplusp, ritzr, 1, workl(kplusp**2+1), 1)
call dcopy(kplusp, ritzi, 1, workl(kplusp**2+kplusp+1), 1)
call dcopy(kplusp, bounds, 1, workl(kplusp**2+2*kplusp+1), 1)
!
!        %---------------------------------------------------%
!        | Select the wanted Ritz values and their bounds    |
!        | to be used in the convergence test.               |
!        | The wanted part of the spectrum and corresponding |
!        | error bounds are in the last NEV loc. of RITZR,   |
!        | RITZI and BOUNDS respectively. The variables NEV  |
!        | and NP may be updated if the NEV-th wanted Ritz   |
!        | value has a non zero imaginary part. In this case |
!        | NEV is increased by one and NP decreased by one.  |
!        | NOTE: The last two arguments of dngets are no     |
!        | longer used as of version 2.1.                    |
!        %---------------------------------------------------%
!
nev = nev0
np = np0
numcnv = nev
call dngets (ishift, which, nev, np, ritzr, ritzi, &
& bounds)
if (nev .eq. nev0+1) numcnv = nev0+1
!
!        %-------------------%
!        | Convergence test. |
!        %-------------------%
!
call dcopy (nev, bounds(np+1), 1, workl(2*np+1), 1)
call dnconv (nev, ritzr(np+1), ritzi(np+1), workl(2*np+1), &
& tol, nconv)
!
!        %---------------------------------------------------------%
!        | Count the number of unwanted Ritz values that have zero |
!        | Ritz estimates. If any Ritz estimates are equal to zero |
!        | then a leading block of H of order equal to at least    |
!        | the number of Ritz values with zero Ritz estimates has  |
!        | split off. None of these Ritz values may be removed by  |
!        | shifting. Decrease NP the number of shifts to apply. If |
!        | no shifts may be applied, then prepare to exit          |
!        %---------------------------------------------------------%
!
nptemp = np
do 30 j=1, nptemp
!            if (bounds(j) .eq. zero) then
if (abs( bounds(j) ) .le. zero) then
np = np - 1
nev = nev + 1
end if
30 continue
!
if ( (nconv .ge. numcnv) .or. &
& (iter .gt. mxiter) .or. &
& (np .eq. 0) ) then
!
!
!           %------------------------------------------------%
!           | Prepare to exit. Put the converged Ritz values |
!           | and corresponding bounds in RITZ(1:NCONV) and  |
!           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
!           | careful when NCONV > NP                        |
!           %------------------------------------------------%
!
!           %------------------------------------------%
!           |  Use h( 3,1 ) as storage to communicate  |
!           |  rnorm to _neupd if needed               |
!           %------------------------------------------%

h(3,1) = rnorm
!
!           %----------------------------------------------%
!           | To be consistent with dngets, we first do a  |
!           | pre-processing sort in order to keep complex |
!           | conjugate pairs together.  This is similar   |
!           | to the pre-processing sort used in dngets    |
!           | except that the sort is done in the opposite |
!           | order.                                       |
!           %----------------------------------------------%
!
if (which .eq. 'LM') wprime = 'SR'
if (which .eq. 'SM') wprime = 'LR'
if (which .eq. 'LR') wprime = 'SM'
if (which .eq. 'SR') wprime = 'LM'
if (which .eq. 'LI') wprime = 'SM'
if (which .eq. 'SI') wprime = 'LM'
!
call dsortc (wprime, dnatrue, kplusp, ritzr, ritzi, bounds)
!
!           %----------------------------------------------%
!           | Now sort Ritz values so that converged Ritz  |
!           | values appear within the first NEV locations |
!           | of ritzr, ritzi and bounds, and the most     |
!           | desired one appears at the front.            |
!           %----------------------------------------------%
!
if (which .eq. 'LM') wprime = 'SM'
if (which .eq. 'SM') wprime = 'LM'
if (which .eq. 'LR') wprime = 'SR'
if (which .eq. 'SR') wprime = 'LR'
if (which .eq. 'LI') wprime = 'SI'
if (which .eq. 'SI') wprime = 'LI'
!
call dsortc(wprime, dnatrue, kplusp, ritzr, ritzi, bounds)
!
!           %--------------------------------------------------%
!           | Scale the Ritz estimate of each Ritz value       |
!           | by 1 / max(eps23,magnitude of the Ritz value).   |
!           %--------------------------------------------------%
!
do 35 j = 1, numcnv
temp = max(eps23,dlapy2(ritzr(j), &
& ritzi(j)))
bounds(j) = bounds(j)/temp
35 continue
!
!           %----------------------------------------------------%
!           | Sort the Ritz values according to the scaled Ritz  |
!           | esitmates.  This will push all the converged ones  |
!           | towards the front of ritzr, ritzi, bounds          |
!           | (in the case when NCONV < NEV.)                    |
!           %----------------------------------------------------%
!
wprime = 'LR'
call dsortc(wprime, dnatrue, numcnv, bounds, ritzr, ritzi)
!
!           %----------------------------------------------%
!           | Scale the Ritz estimate back to its original |
!           | value.                                       |
!           %----------------------------------------------%
!
do 40 j = 1, numcnv
temp = max(eps23, dlapy2(ritzr(j), &
& ritzi(j)))
bounds(j) = bounds(j)*temp
40 continue
!
!           %------------------------------------------------%
!           | Sort the converged Ritz values again so that   |
!           | the "threshold" value appears at the front of  |
!           | ritzr, ritzi and bound.                        |
!           %------------------------------------------------%
!
call dsortc(which, dnatrue, nconv, ritzr, ritzi, bounds)
!
!
!           %------------------------------------%
!           | Max iterations have been exceeded. |
!           %------------------------------------%
!
if (iter .gt. mxiter .and. nconv .lt. numcnv) info = 1
!
!           %---------------------%
!           | No shifts to apply. |
!           %---------------------%
!
if (np .eq. 0 .and. nconv .lt. numcnv) info = 2
!
np = nconv
go to 1100
!
else if ( (nconv .lt. numcnv) .and. (ishift .eq. 1) ) then
!
!           %-------------------------------------------------%
!           | Do not have all the requested eigenvalues yet.  |
!           | To prevent possible stagnation, adjust the size |
!           | of NEV.                                         |
!           %-------------------------------------------------%
!
nevbef = nev
nev = nev + min(nconv, np/2)
if (nev .eq. 1 .and. kplusp .ge. 6) then
nev = kplusp / 2
else if (nev .eq. 1 .and. kplusp .gt. 3) then
nev = 2
end if
np = kplusp - nev
!
!           %---------------------------------------%
!           | If the size of NEV was just increased |
!           | resort the eigenvalues.               |
!           %---------------------------------------%
!
if (nevbef .lt. nev) &
& call dngets (ishift, which, nev, np, ritzr, ritzi, &
& bounds)
!
end if
!
if (ishift .eq. 0) then
!
!           %-------------------------------------------------------%
!           | User specified shifts: reverse comminucation to       |
!           | compute the shifts. They are returned in the first    |
!           | 2*NP locations of WORKL.                              |
!           %-------------------------------------------------------%
!
ushift = .true.
ido = 3
go to 9000
end if
!
50 continue
!
!        %------------------------------------%
!        | Back from reverse communication;   |
!        | User specified shifts are returned |
!        | in WORKL(1:2*NP)                   |
!        %------------------------------------%
!
ushift = .false.
!
if ( ishift .eq. 0 ) then
!
!            %----------------------------------%
!            | Move the NP shifts from WORKL to |
!            | RITZR, RITZI to free up WORKL    |
!            | for non-exact shift case.        |
!            %----------------------------------%
!
call dcopy (np, workl,       1, ritzr(1), 1)
call dcopy (np, workl(np+1), 1, ritzi(1), 1)
end if
!
!
!        %---------------------------------------------------------%
!        | Apply the NP implicit shifts by QR bulge chasing.       |
!        | Each shift is applied to the whole upper Hessenberg     |
!        | matrix H.                                               |
!        | The first 2*N locations of WORKD are used as workspace. |
!        %---------------------------------------------------------%
!
call dnapps (n, nev, np, ritzr, ritzi, v, ldv, &
& h, ldh, resid, q, ldq, workl, workd)
!
!        %---------------------------------------------%
!        | Compute the B-norm of the updated residual. |
!        | Keep B*RESID in WORKD(1:N) to be used in    |
!        | the first step of the next call to dnaitr.  |
!        %---------------------------------------------%
!
cnorm = .true.
if (bmat .eq. 'G') then
call dcopy (n, resid, 1, workd(n+1), 1)
ipntr(1) = n + 1
ipntr(2) = 1
ido = 2
!
!           %----------------------------------%
!           | Exit in order to compute B*RESID |
!           %----------------------------------%
!
go to 9000
else if (bmat .eq. 'I') then
call dcopy (n, resid, 1, workd(1), 1)
end if
!
100 continue
!
!        %----------------------------------%
!        | Back from reverse communication; |
!        | WORKD(1:N) := B*RESID            |
!        %----------------------------------%
!
if (bmat .eq. 'G') then
rnorm = ddot (n, resid, 1, workd, 1)
rnorm = sqrt(abs(rnorm))
else if (bmat .eq. 'I') then
rnorm = dnrm2(n, resid, 1)
end if
cnorm = .false.
!
go to 1000
!
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
!
1100 continue
!
mxiter = iter
nev = numcnv
!
1200 continue
ido = 99
!
!     %------------%
!     | Error Exit |
!     %------------%
!
9000 continue
!
!     %---------------%
!     | End of dnaup2 |
!     %---------------%
!
return
end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dnapps
!
!\Description:
!  Given the Arnoldi factorization
!
!     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
!
!  apply NP implicit shifts resulting in
!
!     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
!
!  where Q is an orthogonal matrix which is the product of rotations
!  and reflections resulting from the NP bulge chage sweeps.
!  The updated Arnoldi factorization becomes:
!
!     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
!
!\Usage:
!  call dnapps
!     ( N, KEV, NP, SHIFTR, SHIFTI, V, LDV, H, LDH, RESID, Q, LDQ,
!       WORKL, WORKD )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Problem size, i.e. size of matrix A.
!
!  KEV     Integer.  (INPUT/OUTPUT)
!          KEV+NP is the size of the input matrix H.
!          KEV is the size of the updated matrix HNEW.  KEV is only
!          updated on ouput when fewer than NP shifts are applied in
!          order to keep the conjugate pair together.
!
!  NP      Integer.  (INPUT)
!          Number of implicit shifts to be applied.
!
!  SHIFTR, Double precision array of length NP.  (INPUT)
!  SHIFTI  Real and imaginary part of the shifts to be applied.
!          Upon, entry to dnapps, the shifts must be sorted so that the
!          conjugate pairs are in consecutive locations.
!
!  V       Double precision N by (KEV+NP) array.  (INPUT/OUTPUT)
!          On INPUT, V contains the current KEV+NP Arnoldi vectors.
!          On OUTPUT, V contains the updated KEV Arnoldi vectors
!          in the first KEV columns of V.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Double precision (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
!          On INPUT, H contains the current KEV+NP by KEV+NP upper
!          Hessenber matrix of the Arnoldi factorization.
!          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
!          matrix in the KEV leading submatrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RESID   Double precision array of length N.  (INPUT/OUTPUT)
!          On INPUT, RESID contains the the residual vector r_{k+p}.
!          On OUTPUT, RESID is the update residual vector rnew_{k}
!          in the first KEV locations.
!
!  Q       Double precision KEV+NP by KEV+NP work array.  (WORKSPACE)
!          Work array used to accumulate the rotations and reflections
!          during the bulge chase sweep.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Double precision work array of length (KEV+NP).  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
!
!  WORKD   Double precision work array of length 2*N.  (WORKSPACE)
!          Distributed array used in the application of the accumulated
!          orthogonal matrix Q.
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
!
!\Routines called:
!     ivout   ARPACK utility routine that prints integers.
!     second  ARPACK utility routine for timing.
!     dmout   ARPACK utility routine that prints matrices.
!     dvout   ARPACK utility routine that prints vectors.
!     dlabad  LAPACK routine that computes machine constants.
!     dlacpy  LAPACK matrix copy routine.
!     dlamch  LAPACK routine that determines machine constants.
!     dlanhs  LAPACK routine that computes various norms of a matrix.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dlarf   LAPACK routine that applies Householder reflection to
!             a matrix.
!     dlarfg  LAPACK Householder reflection construction routine.
!     dlartg  LAPACK Givens rotation construction routine.
!     dlaset  LAPACK matrix initialization routine.
!     dgemv   Level 2 BLAS routine for matrix vector multiplication.
!     daxpy   Level 1 BLAS that computes a vector triad.
!     dcopy   Level 1 BLAS that copies one vector to another .
!     dscal   Level 1 BLAS that scales a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.4'
!
!\SCCS Information: @(#)
! FILE: napps.F   SID: 2.4   DATE OF SID: 3/28/97   RELEASE: 2
!
!\Remarks
!  1. In this version, each shift is applied to all the sublocks of
!     the Hessenberg matrix H and not just to the submatrix that it
!     comes from. Deflation as in LAPACK routine dlahqr (QR algorithm
!     for upper Hessenberg matrices ) is used.
!     The subdiagonals of H are enforced to be non-negative.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dnapps &
& ( n, kev, np, shiftr, shifti, v, ldv, h, ldh, resid, q, ldq, &
& workl, workd )

implicit none
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
integer    kev, ldh, ldq, ldv, n, np
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
Double precision &
& h(ldh,kev+np), resid(n), shifti(np), shiftr(np), &
& v(ldv,kev+np), q(ldq,kev+np), workd(2*n), workl(kev+np)
!
!     %------------%
!     | Parameters |
!     %------------%
!
Double precision &
& one, zero
parameter (one = 1.0D+0, zero = 0.0D+0)
integer(4)  ione
parameter (ione = 1)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
integer    i, iend, ir, istart, j, jj, kplusp, nr
logical    cconj, first
Double precision &
& c, f, g, h11, h12, h21, h22, h32, ovfl, r, s, sigmai, &
& sigmar, smlnum, ulp, unfl, u(3), t, tau, tst1
save       first, ovfl, smlnum, ulp, unfl
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
external   daxpy, dcopy, dscal, dlacpy, dlarfg, dlarf, &
& dlaset, dlabad, dlartg
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
!
Double precision &
& dlamch, dlanhs, dlapy2
external   dlamch, dlanhs, dlapy2
!
!     %----------------------%
!     | Intrinsics Functions |
!     %----------------------%
!
intrinsic  abs, max, min
!
!     %----------------%
!     | Data statments |
!     %----------------%
!
data       first / .true. /
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
if (first) then
!
!        %-----------------------------------------------%
!        | Set machine-dependent constants for the       |
!        | stopping criterion. If norm(H) <= sqrt(OVFL), |
!        | overflow should not occur.                    |
!        | REFERENCE: LAPACK subroutine dlahqr           |
!        %-----------------------------------------------%
!
unfl = dlamch( 'safe minimum' )
ovfl = one / unfl
call dlabad( unfl, ovfl )
ulp = dlamch( 'precision' )
smlnum = unfl*( n / ulp )
first = .false.
end if
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
!      call second (t0)
kplusp = kev + np
!
!     %--------------------------------------------%
!     | Initialize Q to the identity to accumulate |
!     | the rotations and reflections              |
!     %--------------------------------------------%
!
call dlaset ('All', kplusp, kplusp, zero, one, q, ldq)
!
!     %----------------------------------------------%
!     | Quick return if there are no shifts to apply |
!     %----------------------------------------------%
!
if (np .eq. 0) go to 9000
!
!     %----------------------------------------------%
!     | Chase the bulge with the application of each |
!     | implicit shift. Each shift is applied to the |
!     | whole matrix including each block.           |
!     %----------------------------------------------%
!
cconj = .false.
do 110 jj = 1, np
sigmar = shiftr(jj)
sigmai = shifti(jj)
!
!        %-------------------------------------------------%
!        | The following set of conditionals is necessary  |
!        | in order that complex conjugate pairs of shifts |
!        | are applied together or not at all.             |
!        %-------------------------------------------------%
!
if ( cconj ) then
!
!           %-----------------------------------------%
!           | cconj = .true. means the previous shift |
!           | had non-zero imaginary part.            |
!           %-----------------------------------------%
!
cconj = .false.
go to 110
else if ( jj .lt. np .and. abs( sigmai ) .gt. zero ) then
!
!           %------------------------------------%
!           | Start of a complex conjugate pair. |
!           %------------------------------------%
!
cconj = .true.
else if ( jj .eq. np .and. abs( sigmai ) .gt. zero ) then
!
!           %----------------------------------------------%
!           | The last shift has a nonzero imaginary part. |
!           | Don't apply it; thus the order of the        |
!           | compressed H is order KEV+1 since only np-1  |
!           | were applied.                                |
!           %----------------------------------------------%
!
kev = kev + 1
go to 110
end if
istart = 1
20 continue
!
!        %--------------------------------------------------%
!        | if sigmai = 0 then                               |
!        |    Apply the jj-th shift ...                     |
!        | else                                             |
!        |    Apply the jj-th and (jj+1)-th together ...    |
!        |    (Note that jj < np at this point in the code) |
!        | end                                              |
!        | to the current block of H. The next do loop      |
!        | determines the current block ;                   |
!        %--------------------------------------------------%
!
do 30 i = istart, kplusp-1
!
!           %----------------------------------------%
!           | Check for splitting and deflation. Use |
!           | a standard test as in the QR algorithm |
!           | REFERENCE: LAPACK subroutine dlahqr    |
!           %----------------------------------------%
!
tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
!            if( tst1.eq.zero )
if( abs(tst1) .le. zero ) &
& tst1 = dlanhs( '1', kplusp-jj+1, h, ldh, workl )
if( abs( h( i+1,i ) ).le.max( ulp*tst1, smlnum ) ) then
iend = i
h(i+1,i) = zero
go to 40
end if
30 continue
iend = kplusp
40 continue
!
!        %------------------------------------------------%
!        | No reason to apply a shift to block of order 1 |
!        %------------------------------------------------%
!
if ( istart .eq. iend ) go to 100
!
!        %------------------------------------------------------%
!        | If istart + 1 = iend then no reason to apply a       |
!        | complex conjugate pair of shifts on a 2 by 2 matrix. |
!        %------------------------------------------------------%
!
if ( istart + 1 .eq. iend .and. abs( sigmai ) .gt. zero ) &
& go to 100
!
h11 = h(istart,istart)
h21 = h(istart+1,istart)
if ( abs( sigmai ) .le. zero ) then
!
!           %---------------------------------------------%
!           | Real-valued shift ==> apply single shift QR |
!           %---------------------------------------------%
!
f = h11 - sigmar
g = h21
!
do 80 i = istart, iend-1
!
!              %-----------------------------------------------------%
!              | Contruct the plane rotation G to zero out the bulge |
!              %-----------------------------------------------------%
!
call dlartg (f, g, c, s, r)
if (i .gt. istart) then
!
!                 %-------------------------------------------%
!                 | The following ensures that h(1:iend-1,1), |
!                 | the first iend-2 off diagonal of elements |
!                 | H, remain non negative.                   |
!                 %-------------------------------------------%
!
if (r .lt. zero) then
r = -r
c = -c
s = -s
end if
h(i,i-1) = r
h(i+1,i-1) = zero
end if
!
!              %---------------------------------------------%
!              | Apply rotation to the left of H;  H <- G'*H |
!              %---------------------------------------------%
!
do 50 j = i, kplusp
t        =  c*h(i,j) + s*h(i+1,j)
h(i+1,j) = -s*h(i,j) + c*h(i+1,j)
h(i,j)   = t
50 continue
!
!              %---------------------------------------------%
!              | Apply rotation to the right of H;  H <- H*G |
!              %---------------------------------------------%
!
do 60 j = 1, min(i+2,iend)
t        =  c*h(j,i) + s*h(j,i+1)
h(j,i+1) = -s*h(j,i) + c*h(j,i+1)
h(j,i)   = t
60 continue
!
!              %----------------------------------------------------%
!              | Accumulate the rotation in the matrix Q;  Q <- Q*G |
!              %----------------------------------------------------%
!
do 70 j = 1, min( i+jj, kplusp )
t        =   c*q(j,i) + s*q(j,i+1)
q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
q(j,i)   = t
70 continue
!
!              %---------------------------%
!              | Prepare for next rotation |
!              %---------------------------%
!
if (i .lt. iend-1) then
f = h(i+1,i)
g = h(i+2,i)
end if
80 continue
!
!           %-----------------------------------%
!           | Finished applying the real shift. |
!           %-----------------------------------%
!
else
!
!           %----------------------------------------------------%
!           | Complex conjugate shifts ==> apply double shift QR |
!           %----------------------------------------------------%
!
h12 = h(istart,istart+1)
h22 = h(istart+1,istart+1)
h32 = h(istart+2,istart+1)
!
!           %---------------------------------------------------------%
!           | Compute 1st column of (H - shift*I)*(H - conj(shift)*I) |
!           %---------------------------------------------------------%
!
s    = 2.0*sigmar
t = dlapy2 ( sigmar, sigmai )
u(1) = ( h11 * (h11 - s) + t * t ) / h21 + h12
u(2) = h11 + h22 - s
u(3) = h32
!
do 90 i = istart, iend-1
!
nr = min ( 3, iend-i+1 )
!
!              %-----------------------------------------------------%
!              | Construct Householder reflector G to zero out u(1). |
!              | G is of the form I - tau*( 1 u )' * ( 1 u' ).       |
!              %-----------------------------------------------------%
!
call dlarfg ( nr, u(1), u(2), 1, tau )
!
if (i .gt. istart) then
h(i,i-1)   = u(1)
h(i+1,i-1) = zero
if (i .lt. iend-1) h(i+2,i-1) = zero
end if
u(1) = one
!
!              %--------------------------------------%
!              | Apply the reflector to the left of H |
!              %--------------------------------------%
!
call dlarf ('Left', nr, kplusp-i+1, u, 1, tau, &
& h(i,i), ldh, workl)
!
!              %---------------------------------------%
!              | Apply the reflector to the right of H |
!              %---------------------------------------%
!
ir = min ( i+3, iend )
call dlarf ('Right', ir, nr, u, 1, tau, &
& h(1,i), ldh, workl)
!
!              %-----------------------------------------------------%
!              | Accumulate the reflector in the matrix Q;  Q <- Q*G |
!              %-----------------------------------------------------%
!
call dlarf ('Right', kplusp, nr, u, 1, tau, &
& q(1,i), ldq, workl)
!
!              %----------------------------%
!              | Prepare for next reflector |
!              %----------------------------%
!
if (i .lt. iend-1) then
u(1) = h(i+1,i)
u(2) = h(i+2,i)
if (i .lt. iend-2) u(3) = h(i+3,i)
end if
!
90 continue
!
!           %--------------------------------------------%
!           | Finished applying a complex pair of shifts |
!           | to the current block                       |
!           %--------------------------------------------%
!
end if
!
100 continue
!
!        %---------------------------------------------------------%
!        | Apply the same shift to the next block if there is any. |
!        %---------------------------------------------------------%
!
istart = iend + 1
if (iend .lt. kplusp) go to 20
!
!        %---------------------------------------------%
!        | Loop back to the top to get the next shift. |
!        %---------------------------------------------%
!
110 continue
!
!     %--------------------------------------------------%
!     | Perform a similarity transformation that makes   |
!     | sure that H will have non negative sub diagonals |
!     %--------------------------------------------------%
!
do 120 j=1,kev
if ( h(j+1,j) .lt. zero ) then
call dscal( kplusp-j+1, -one, h(j+1,j), ldh)
call dscal( min(j+2, kplusp), -one, h(1,j+1), &
& ione)
call dscal( min(j+np+1,kplusp), -one, q(1,j+1), &
& ione)
!v              call dscal( min(j+2, kplusp), -one, h(1,j+1), 1 )
!v              call dscal( min(j+np+1,kplusp), -one, q(1,j+1), 1 )
end if
120 continue
!
do 130 i = 1, kev
!
!        %--------------------------------------------%
!        | Final check for splitting and deflation.   |
!        | Use a standard test as in the QR algorithm |
!        | REFERENCE: LAPACK subroutine dlahqr        |
!        %--------------------------------------------%
!
tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
!         if( tst1.eq.zero )
if( abs(tst1) .le. zero) &
& tst1 = dlanhs( '1', kev, h, ldh, workl )
if( h( i+1,i ) .le. max( ulp*tst1, smlnum ) ) &
& h(i+1,i) = zero
130 continue
!
!     %-------------------------------------------------%
!     | Compute the (kev+1)-st column of (V*Q) and      |
!     | temporarily store the result in WORKD(N+1:2*N). |
!     | This is needed in the residual update since we  |
!     | cannot GUARANTEE that the corresponding entry   |
!     | of H would be zero as in exact arithmetic.      |
!     %-------------------------------------------------%
!
if (h(kev+1,kev) .gt. zero) &
& call dgemv ('N', n, kplusp, one, v, ldv, q(1,kev+1), 1, zero, &
& workd(n+1), 1)
!
!     %----------------------------------------------------------%
!     | Compute column 1 to kev of (V*Q) in backward order       |
!     | taking advantage of the upper Hessenberg structure of Q. |
!     %----------------------------------------------------------%
!
do 140 i = 1, kev
call dgemv ('N', n, kplusp-i+1, one, v, ldv, &
& q(1,kev-i+1), 1, zero, workd(1), 1)
call dcopy (n, workd, 1, v(1,kplusp-i+1), 1)
140 continue
!
!     %-------------------------------------------------%
!     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
!     %-------------------------------------------------%
!
call dlacpy ('A', n, kev, v(1,kplusp-kev+1), ldv, v, ldv)
!
!     %--------------------------------------------------------------%
!     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
!     %--------------------------------------------------------------%
!
if (h(kev+1,kev) .gt. zero) &
& call dcopy (n, workd(n+1), 1, v(1,kev+1), 1)
!
!     %-------------------------------------%
!     | Update the residual vector:         |
!     |    r <- sigmak*r + betak*v(:,kev+1) |
!     | where                               |
!     |    sigmak = (e_{kplusp}'*Q)*e_{kev} |
!     |    betak = e_{kev+1}'*H*e_{kev}     |
!     %-------------------------------------%
!
call dscal (n, q(kplusp,kev), resid(1), ione)
if (h(kev+1,kev) .gt. zero) &
& call daxpy (n, h(kev+1,kev), v(1,kev+1), 1, resid, 1)
!
9000 continue
!
return
!
!     %---------------%
!     | End of dnapps |
!     %---------------%
!
end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dnconv
!
!\Description:
!  Convergence testing for the nonsymmetric Arnoldi eigenvalue routine.
!
!\Usage:
!  call dnconv
!     ( N, RITZR, RITZI, BOUNDS, TOL, NCONV )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Number of Ritz values to check for convergence.
!
!  RITZR,  Double precision arrays of length N.  (INPUT)
!  RITZI   Real and imaginary parts of the Ritz values to be checked
!          for convergence.

!  BOUNDS  Double precision array of length N.  (INPUT)
!          Ritz estimates for the Ritz values in RITZR and RITZI.
!
!  TOL     Double precision scalar.  (INPUT)
!          Desired backward error for a Ritz value to be considered
!          "converged".
!
!  NCONV   Integer scalar.  (OUTPUT)
!          Number of "converged" Ritz values.
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
!\Routines called:
!     second  ARPACK utility routine for timing.
!     dlamch  LAPACK routine that determines machine constants.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.1'
!
!\SCCS Information: @(#)
! FILE: nconv.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
!
!\Remarks
!     1. xxxx
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dnconv (n, ritzr, ritzi, bounds, tol, nconv)

implicit none
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
integer    n, nconv
Double precision &
& tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%

Double precision &
& ritzr(n), ritzi(n), bounds(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
integer    i
Double precision &
& temp, eps23
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
Double precision &
& dlapy2, dlamch
external   dlapy2, dlamch

!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------------------%
!     | Convergence test: unlike in the symmetric code, I am not    |
!     | using things like refined error bounds and gap condition    |
!     | because I don't know the exact equivalent concept.          |
!     |                                                             |
!     | Instead the i-th Ritz value is considered "converged" when: |
!     |                                                             |
!     |     bounds(i) .le. ( TOL * | ritz | )                       |
!     |                                                             |
!     | for some appropriate choice of norm.                        |
!     %-------------------------------------------------------------%
!
!     %---------------------------------%
!     | Get machine dependent constant. |
!     %---------------------------------%
!
eps23 = dlamch('Epsilon-Machine')
eps23 = eps23**(2.0D+0 / 3.0D+0)
!
nconv  = 0
do 20 i = 1, n
temp = max( eps23, dlapy2( ritzr(i), ritzi(i) ) )
if (bounds(i) .le. tol*temp)   nconv = nconv + 1
20 continue
!
return
!
!     %---------------%
!     | End of dnconv |
!     %---------------%
!
end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsortc
!
!\Description:
!  Sorts the complex array in XREAL and XIMAG into the order
!  specified by WHICH and optionally applies the permutation to the
!  real array Y. It is assumed that if an element of XIMAG is
!  nonzero, then its negative is also an element. In other words,
!  both members of a complex conjugate pair are to be sorted and the
!  pairs are kept adjacent to each other.
!
!\Usage:
!  call dsortc
!     ( WHICH, APPLY, N, XREAL, XIMAG, Y )
!
!\Arguments
!  WHICH   Character*2.  (Input)
!          'LM' -> sort XREAL,XIMAG into increasing order of magnitude.
!          'SM' -> sort XREAL,XIMAG into decreasing order of magnitude.
!          'LR' -> sort XREAL into increasing order of algebraic.
!          'SR' -> sort XREAL into decreasing order of algebraic.
!          'LI' -> sort XIMAG into increasing order of magnitude.
!          'SI' -> sort XIMAG into decreasing order of magnitude.
!          NOTE: If an element of XIMAG is non-zero, then its negative
!                is also an element.
!
!  APPLY   Logical.  (Input)
!          APPLY = .TRUE.  -> apply the sorted order to array Y.
!          APPLY = .FALSE. -> do not apply the sorted order to array Y.
!
!  N       Integer.  (INPUT)
!          Size of the arrays.
!
!  XREAL,  Double precision array of length N.  (INPUT/OUTPUT)
!  XIMAG   Real and imaginary part of the array to be sorted.
!
!  Y       Double precision array of length N.  (INPUT/OUTPUT)
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.1'
!               Adapted from the sort routine in LANSO.
!
!\SCCS Information: @(#)
! FILE: sortc.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dsortc (which, apply, n, xreal, ximag, y)
!
implicit none
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
character  which*2
logical    apply
integer    n
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
Double precision &
& xreal(0:n-1), ximag(0:n-1), y(0:n-1)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
integer    i, igap, j
Double precision &
& temp, temp1, temp2
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
Double precision &
& dlapy2
external   dlapy2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
igap = n / 2
!
if (which .eq. 'LM') then
!
!        %------------------------------------------------------%
!        | Sort XREAL,XIMAG into increasing order of magnitude. |
!        %------------------------------------------------------%
!
10 continue
if (igap .eq. 0) go to 9000
!
do 30 i = igap, n-1
j = i-igap
20 continue
!
if (j.lt.0) go to 30
!
temp1 = dlapy2(xreal(j),ximag(j))
temp2 = dlapy2(xreal(j+igap),ximag(j+igap))
!
if (temp1.gt.temp2) then
temp = xreal(j)
xreal(j) = xreal(j+igap)
xreal(j+igap) = temp
!
temp = ximag(j)
ximag(j) = ximag(j+igap)
ximag(j+igap) = temp
!
if (apply) then
temp = y(j)
y(j) = y(j+igap)
y(j+igap) = temp
end if
else
go to 30
end if
j = j-igap
go to 20
30 continue
igap = igap / 2
go to 10
!
else if (which .eq. 'SM') then
!
!        %------------------------------------------------------%
!        | Sort XREAL,XIMAG into decreasing order of magnitude. |
!        %------------------------------------------------------%
!
40 continue
if (igap .eq. 0) go to 9000
!
do 60 i = igap, n-1
j = i-igap
50 continue
!
if (j .lt. 0) go to 60
!
temp1 = dlapy2(xreal(j),ximag(j))
temp2 = dlapy2(xreal(j+igap),ximag(j+igap))
!
if (temp1.lt.temp2) then
temp = xreal(j)
xreal(j) = xreal(j+igap)
xreal(j+igap) = temp
!
temp = ximag(j)
ximag(j) = ximag(j+igap)
ximag(j+igap) = temp
!
if (apply) then
temp = y(j)
y(j) = y(j+igap)
y(j+igap) = temp
end if
else
go to 60
endif
j = j-igap
go to 50
60 continue
igap = igap / 2
go to 40
!
else if (which .eq. 'LR') then
!
!        %------------------------------------------------%
!        | Sort XREAL into increasing order of algebraic. |
!        %------------------------------------------------%
!
70 continue
if (igap .eq. 0) go to 9000
!
do 90 i = igap, n-1
j = i-igap
80 continue
!
if (j.lt.0) go to 90
!
if (xreal(j).gt.xreal(j+igap)) then
temp = xreal(j)
xreal(j) = xreal(j+igap)
xreal(j+igap) = temp
!
temp = ximag(j)
ximag(j) = ximag(j+igap)
ximag(j+igap) = temp
!
if (apply) then
temp = y(j)
y(j) = y(j+igap)
y(j+igap) = temp
end if
else
go to 90
endif
j = j-igap
go to 80
90 continue
igap = igap / 2
go to 70
!
else if (which .eq. 'SR') then
!
!        %------------------------------------------------%
!        | Sort XREAL into decreasing order of algebraic. |
!        %------------------------------------------------%
!
100 continue
if (igap .eq. 0) go to 9000
do 120 i = igap, n-1
j = i-igap
110 continue
!
if (j.lt.0) go to 120
!
if (xreal(j).lt.xreal(j+igap)) then
temp = xreal(j)
xreal(j) = xreal(j+igap)
xreal(j+igap) = temp
!
temp = ximag(j)
ximag(j) = ximag(j+igap)
ximag(j+igap) = temp
!
if (apply) then
temp = y(j)
y(j) = y(j+igap)
y(j+igap) = temp
end if
else
go to 120
endif
j = j-igap
go to 110
120 continue
igap = igap / 2
go to 100
!
else if (which .eq. 'LI') then
!
!        %------------------------------------------------%
!        | Sort XIMAG into increasing order of magnitude. |
!        %------------------------------------------------%
!
130 continue
if (igap .eq. 0) go to 9000
do 150 i = igap, n-1
j = i-igap
140 continue
!
if (j.lt.0) go to 150
!
if (abs(ximag(j)).gt.abs(ximag(j+igap))) then
temp = xreal(j)
xreal(j) = xreal(j+igap)
xreal(j+igap) = temp
!
temp = ximag(j)
ximag(j) = ximag(j+igap)
ximag(j+igap) = temp
!
if (apply) then
temp = y(j)
y(j) = y(j+igap)
y(j+igap) = temp
end if
else
go to 150
endif
j = j-igap
go to 140
150 continue
igap = igap / 2
go to 130
!
else if (which .eq. 'SI') then
!
!        %------------------------------------------------%
!        | Sort XIMAG into decreasing order of magnitude. |
!        %------------------------------------------------%
!
160 continue
if (igap .eq. 0) go to 9000
do 180 i = igap, n-1
j = i-igap
170 continue
!
if (j.lt.0) go to 180
!
if (abs(ximag(j)).lt.abs(ximag(j+igap))) then
temp = xreal(j)
xreal(j) = xreal(j+igap)
xreal(j+igap) = temp
!
temp = ximag(j)
ximag(j) = ximag(j+igap)
ximag(j+igap) = temp
!
if (apply) then
temp = y(j)
y(j) = y(j+igap)
y(j+igap) = temp
end if
else
go to 180
endif
j = j-igap
go to 170
180 continue
igap = igap / 2
go to 160
end if
!
9000 continue
return
!
!     %---------------%
!     | End of dsortc |
!     %---------------%
!
end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dneigh
!
!\Description:
!  Compute the eigenvalues of the current upper Hessenberg matrix
!  and the corresponding Ritz estimates given the current residual norm.
!
!\Usage:
!  call dneigh
!     ( RNORM, N, H, LDH, RITZR, RITZI, BOUNDS, Q, LDQ, WORKL, IERR )
!
!\Arguments
!  RNORM   Double precision scalar.  (INPUT)
!          Residual norm corresponding to the current upper Hessenberg
!          matrix H.
!
!  N       Integer.  (INPUT)
!          Size of the matrix H.
!
!  H       Double precision N by N array.  (INPUT)
!          H contains the current upper Hessenberg matrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RITZR,  Double precision arrays of length N.  (OUTPUT)
!  RITZI   On output, RITZR(1:N) (resp. RITZI(1:N)) contains the real
!          (respectively imaginary) parts of the eigenvalues of H.
!
!  BOUNDS  Double precision array of length N.  (OUTPUT)
!          On output, BOUNDS contains the Ritz estimates associated with
!          the eigenvalues RITZR and RITZI.  This is equal to RNORM
!          times the last components of the eigenvectors corresponding
!          to the eigenvalues in RITZR and RITZI.
!
!  Q       Double precision N by N array.  (WORKSPACE)
!          Workspace needed to store the eigenvectors of H.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Double precision work array of length N**2 + 3*N.  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  This is needed to keep the full Schur form
!          of H and also in the calculation of the eigenvectors of H.
!
!  IERR    Integer.  (OUTPUT)
!          Error exit flag from dlaqrb or dtrevc.
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
!\Routines called:
!     dlaqrb  ARPACK routine to compute the real Schur form of an
!             upper Hessenberg matrix and last row of the Schur vectors.
!     second  ARPACK utility routine for timing.
!     dmout   ARPACK utility routine that prints matrices
!     dvout   ARPACK utility routine that prints vectors.
!     dlacpy  LAPACK matrix copy routine.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dtrevc  LAPACK routine to compute the eigenvectors of a matrix
!             in upper quasi-triangular form
!     dgemv   Level 2 BLAS routine for matrix vector multiplication.
!     dcopy   Level 1 BLAS that copies one vector to another .
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     dscal   Level 1 BLAS that scales a vector.
!
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.1'
!
!\SCCS Information: @(#)
! FILE: neigh.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
!
!\Remarks
!     None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dneigh (rnorm, n, h, ldh, ritzr, ritzi, bounds, &
& q, ldq, workl, ierr)

implicit none
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
integer    ierr, n, ldh, ldq
Double precision &
& rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
Double precision &
& bounds(n), h(ldh,n), q(ldq,n), ritzi(n), ritzr(n), &
& workl(n*(n+3))
!
!     %------------%
!     | Parameters |
!     %------------%
!
Double precision &
& one, zero
parameter (one = 1.0D+0, zero = 0.0D+0)
!
integer(4) &
& ione
parameter (ione = 1)

!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
!     % dneightrue = .true. to adress logical(4)/(8)
!     % initialization for spam/spam64
!     % analogue for dneighone
logical    select(1), dneightrue
integer    i, iconj,  dneighone
Double precision &
& temp, vl(1)
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
external   dcopy, dlacpy, dlaqrb, dtrevc
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
Double precision &
& dlapy2, dnrm2
external   dlapy2, dnrm2
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
intrinsic  abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
dneightrue = .true.
dneighone = 1
!
!
!     %-----------------------------------------------------------%
!     | 1. Compute the eigenvalues, the last components of the    |
!     |    corresponding Schur vectors and the full Schur form T  |
!     |    of the current upper Hessenberg matrix H.              |
!     | dlaqrb returns the full Schur form of H in WORKL(1:N**2)  |
!     | and the last components of the Schur vectors in BOUNDS.   |
!     %-----------------------------------------------------------%
!
!x      call dlacpy ('All', n, n, h, ldh, workl, n)
call dlacpy ('All', n, n, h(1,1), ldh, workl, n)
call dlaqrb (dneightrue, n, dneighone, n, workl, n, ritzr, ritzi, &
& bounds, ierr)
if (ierr .ne. 0) go to 9000
!
!
!     %-----------------------------------------------------------%
!     | 2. Compute the eigenvectors of the full Schur form T and  |
!     |    apply the last components of the Schur vectors to get  |
!     |    the last components of the corresponding eigenvectors. |
!     | Remember that if the i-th and (i+1)-st eigenvalues are    |
!     | complex conjugate pairs, then the real & imaginary part   |
!     | of the eigenvector components are split across adjacent   |
!     | columns of Q.                                             |
!     %-----------------------------------------------------------%
!
call dtrevc ('R', 'A', select, n, workl, n, vl, n, q, ldq, &
& n, n, workl(n*n+1), ierr)
!
if (ierr .ne. 0) go to 9000
!
!     %------------------------------------------------%
!     | Scale the returning eigenvectors so that their |
!     | euclidean norms are all one. LAPACK subroutine |
!     | dtrevc returns each eigenvector normalized so  |
!     | that the element of largest magnitude has      |
!     | magnitude 1; here the magnitude of a complex   |
!     | number (x,y) is taken to be |x| + |y|.         |
!     %------------------------------------------------%
!
iconj = 0
do 10 i=1, n
if ( abs( ritzi(i) ) .le. zero ) then
!
!           %----------------------%
!           | Real eigenvalue case |
!           %----------------------%
!
temp = dnrm2( n, q(1,i), 1 )
call dscal ( n, one / temp, q(1,i), ione )
else
!
!           %-------------------------------------------%
!           | Complex conjugate pair case. Note that    |
!           | since the real and imaginary part of      |
!           | the eigenvector are stored in consecutive |
!           | columns, we further normalize by the      |
!           | square root of two.                       |
!           %-------------------------------------------%
!
if (iconj .eq. 0) then
temp = dlapy2( dnrm2( n, q(1,i), 1 ), &
& dnrm2( n, q(1,i+1), 1 ) )
call dscal ( n, one / temp, q(1,i), &
& ione )
call dscal ( n, one / temp, q(1,i+1), &
& ione )
iconj = 1
else
iconj = 0
end if
end if
10 continue
!
call dgemv('T', n, n, one, q, ldq, bounds(1), 1, zero, &
& workl(1), 1)
!
!
!     %----------------------------%
!     | Compute the Ritz estimates |
!     %----------------------------%
!
iconj = 0
do 20 i = 1, n
if ( abs( ritzi(i) ) .le. zero ) then
!
!           %----------------------%
!           | Real eigenvalue case |
!           %----------------------%
!
bounds(i) = rnorm * abs( workl(i) )
else
!
!           %-------------------------------------------%
!           | Complex conjugate pair case. Note that    |
!           | since the real and imaginary part of      |
!           | the eigenvector are stored in consecutive |
!           | columns, we need to take the magnitude    |
!           | of the last components of the two vectors |
!           %-------------------------------------------%
!
if (iconj .eq. 0) then
bounds(i) = rnorm * dlapy2( workl(i), workl(i+1) )
bounds(i+1) = bounds(i)
iconj = 1
else
iconj = 0
end if
end if
20 continue
!
9000 continue
return
!
!     %---------------%
!     | End of dneigh |
!     %---------------%
!
end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dlaqrb
!
!\Description:
!  Compute the eigenvalues and the Schur decomposition of an upper
!  Hessenberg submatrix in rows and columns ILO to IHI.  Only the
!  last component of the Schur vectors are computed.
!
!  This is mostly a modification of the LAPACK routine dlahqr.
!
!\Usage:
!  call dlaqrb
!     ( WANTT, N, ILO, IHI, H, LDH, WR, WI,  Z, INFO )
!
!\Arguments
!  WANTT   Logical variable.  (INPUT)
!          = .TRUE. : the full Schur form T is required;
!          = .FALSE.: only eigenvalues are required.
!
!  N       Integer.  (INPUT)
!          The order of the matrix H.  N >= 0.
!
!  ILO     Integer.  (INPUT)
!  IHI     Integer.  (INPUT)
!          It is assumed that H is already upper quasi-triangular in
!          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
!          ILO = 1). SLAQRB works primarily with the Hessenberg
!          submatrix in rows and columns ILO to IHI, but applies
!          transformations to all of H if WANTT is .TRUE..
!          1 <= ILO <= max(1,IHI); IHI <= N.
!
!  H       Double precision array, dimension (LDH,N).  (INPUT/OUTPUT)
!          On entry, the upper Hessenberg matrix H.
!          On exit, if WANTT is .TRUE., H is upper quasi-triangular in
!          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in
!          standard form. If WANTT is .FALSE., the contents of H are
!          unspecified on exit.
!
!  LDH     Integer.  (INPUT)
!          The leading dimension of the array H. LDH >= max(1,N).
!
!  WR      Double precision array, dimension (N).  (OUTPUT)
!  WI      Double precision array, dimension (N).  (OUTPUT)
!          The real and imaginary parts, respectively, of the computed
!          eigenvalues ILO to IHI are stored in the corresponding
!          elements of WR and WI. If two eigenvalues are computed as a
!          complex conjugate pair, they are stored in consecutive
!          elements of WR and WI, say the i-th and (i+1)th, with
!          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
!          eigenvalues are stored in the same order as on the diagonal
!          of the Schur form returned in H, with WR(i) = H(i,i), and, if
!          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
!          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
!
!  Z       Double precision array, dimension (N).  (OUTPUT)
!          On exit Z contains the last components of the Schur vectors.
!
!  INFO    Integer.  (OUPUT)
!          = 0: successful exit
!          > 0: SLAQRB failed to compute all the eigenvalues ILO to IHI
!               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
!               elements i+1:ihi of WR and WI contain those eigenvalues
!               which have been successfully computed.
!
!\Remarks
!  1. None.
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     dlabad  LAPACK routine that computes machine constants.
!     dlamch  LAPACK routine that determines machine constants.
!     dlanhs  LAPACK routine that computes various norms of a matrix.
!     dlanv2  LAPACK routine that computes the Schur factorization of
!             2 by 2 nonsymmetric matrix in standard form.
!     dlarfg  LAPACK Householder reflection construction routine.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     drot    Level 1 BLAS that applies a rotation to a 2 by 2 matrix.

!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.4'
!               Modified from the LAPACK routine dlahqr so that only the
!               last component of the Schur vectors are computed.
!
!\SCCS Information: @(#)
! FILE: laqrb.F   SID: 2.2   DATE OF SID: 8/27/96   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dlaqrb ( wantt, n, ilo, ihi, h, ldh, wr, wi, &
& z, info )
!
implicit none
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
logical    wantt
integer    ihi, ilo, info, ldh, n
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
Double precision &
& h( ldh, * ), wi( * ), wr( * ), z( * )
!
!     %------------%
!     | Parameters |
!     %------------%
!
Double precision &
& zero, one, dat1, dat2
parameter (zero = 0.0D+0, one = 1.0D+0, dat1 = 7.5D-1, &
& dat2 = -4.375D-1)
!
integer(4) ione
parameter (ione = 1)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
integer    i, i1, i2, itn, its, j, k, l, m, nh, nr
integer    dlqthree
parameter (dlqthree = 3)
Double precision &
& cs, h00, h10, h11, h12, h21, h22, h33, h33s, &
& h43h34, h44, h44s, ovfl, s, smlnum, sn, sum, &
& t1, t2, t3, tst1, ulp, unfl, v1, v2, v3
Double precision &
& v( 3 ), work( 1 )
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
!
Double precision &
& dlamch, dlanhs
external   dlamch, dlanhs
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
external   dcopy, dlabad, dlanv2, dlarfg, drot
!
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
info = 0
!
!     % if i2 not initialized; Warning: ‘i2’ may be used uninitialized
!     % in this function [-Wmaybe-uninitialized]
i2 = 0
!
!     %--------------------------%
!     | Quick return if possible |
!     %--------------------------%
!
if( n.eq.0 ) &
& return
if( ilo.eq.ihi ) then
wr( ilo ) = h( ilo, ilo )
wi( ilo ) = zero
return
end if
!
!     %---------------------------------------------%
!     | Initialize the vector of last components of |
!     | the Schur vectors for accumulation.         |
!     %---------------------------------------------%
!
do 5 j = 1, n-1
z(j) = zero
5 continue
z(n) = one
!
nh = ihi - ilo + 1
!
!     %-------------------------------------------------------------%
!     | Set machine-dependent constants for the stopping criterion. |
!     | If norm(H) <= sqrt(OVFL), overflow should not occur.        |
!     %-------------------------------------------------------------%
!
unfl = dlamch( 'safe minimum' )
ovfl = one / unfl
call dlabad( unfl, ovfl )
ulp = dlamch( 'precision' )
smlnum = unfl*( nh / ulp )
!
!     %---------------------------------------------------------------%
!     | I1 and I2 are the indices of the first row and last column    |
!     | of H to which transformations must be applied. If eigenvalues |
!     | only are computed, I1 and I2 are set inside the main loop.    |
!     | Zero out H(J+2,J) = ZERO for J=1:N if WANTT = .TRUE.          |
!     | else H(J+2,J) for J=ILO:IHI-ILO-1 if WANTT = .FALSE.          |
!     %---------------------------------------------------------------%
!
if( wantt ) then
i1 = 1
i2 = n
do 8 i=1,i2-2
h(i1+i+1,i) = zero
8 continue
else
do 9 i=1, ihi-ilo-1
h(ilo+i+1,ilo+i-1) = zero
9 continue
end if
!
!     %---------------------------------------------------%
!     | ITN is the total number of QR iterations allowed. |
!     %---------------------------------------------------%
!
itn = 30*nh
!
!     ------------------------------------------------------------------
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
!     H(L,L-1) is negligible so that the matrix splits.
!     ------------------------------------------------------------------
!
i = ihi
10 continue
l = ilo
if( i.lt.ilo ) &
& go to 150

!     %--------------------------------------------------------------%
!     | Perform QR iterations on rows and columns ILO to I until a   |
!     | submatrix of order 1 or 2 splits off at the bottom because a |
!     | subdiagonal element has become negligible.                   |
!     %--------------------------------------------------------------%

do 130 its = 0, itn
!
!        %----------------------------------------------%
!        | Look for a single small subdiagonal element. |
!        %----------------------------------------------%
!
do 20 k = i, l + 1, -1
tst1 = abs( h( k-1, k-1 ) ) + abs( h( k, k ) )
if( abs(tst1) .le. zero) &
& tst1 = dlanhs( '1', i-l+1, h( l, l ), ldh, work )
if( abs( h( k, k-1 ) ).le.max( ulp*tst1, smlnum ) ) &
& go to 30
20 continue
30 continue
l = k
if( l.gt.ilo ) then
!
!           %------------------------%
!           | H(L,L-1) is negligible |
!           %------------------------%
!
h( l, l-1 ) = zero
end if
!
!        %-------------------------------------------------------------%
!        | Exit from loop if a submatrix of order 1 or 2 has split off |
!        %-------------------------------------------------------------%
!
if( l.ge.i-1 ) &
& go to 140
!
!        %---------------------------------------------------------%
!        | Now the active submatrix is in rows and columns L to I. |
!        | If eigenvalues only are being computed, only the active |
!        | submatrix need be transformed.                          |
!        %---------------------------------------------------------%
!
if( .not.wantt ) then
i1 = l
i2 = i
end if
!
if( its.eq.10 .or. its.eq.20 ) then
!
!           %-------------------%
!           | Exceptional shift |
!           %-------------------%
!
s = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
h44 = dat1*s
h33 = h44
h43h34 = dat2*s*s
!
else
!
!           %-----------------------------------------%
!           | Prepare to use Wilkinson's double shift |
!           %-----------------------------------------%
!
h44 = h( i, i )
h33 = h( i-1, i-1 )
h43h34 = h( i, i-1 )*h( i-1, i )
end if
!
!        %-----------------------------------------------------%
!        | Look for two consecutive small subdiagonal elements |
!        %-----------------------------------------------------%
!
do 40 m = i - 2, l, -1
!
!           %---------------------------------------------------------%
!           | Determine the effect of starting the double-shift QR    |
!           | iteration at row M, and see if this would make H(M,M-1) |
!           | negligible.                                             |
!           %---------------------------------------------------------%
!
h11 = h( m, m )
h22 = h( m+1, m+1 )
h21 = h( m+1, m )
h12 = h( m, m+1 )
h44s = h44 - h11
h33s = h33 - h11
v1 = ( h33s*h44s-h43h34 ) / h21 + h12
v2 = h22 - h11 - h33s - h44s
v3 = h( m+2, m+1 )
s = abs( v1 ) + abs( v2 ) + abs( v3 )
v1 = v1 / s
v2 = v2 / s
v3 = v3 / s
v( 1 ) = v1
v( 2 ) = v2
v( 3 ) = v3
if( m.eq.l ) &
& go to 50
h00 = h( m-1, m-1 )
h10 = h( m, m-1 )
tst1 = abs( v1 )*( abs( h00 )+abs( h11 )+abs( h22 ) )
if( abs( h10 )*( abs( v2 )+abs( v3 ) ).le.ulp*tst1 ) &
& go to 50
40 continue
50 continue
!
!        %----------------------%
!        | Double-shift QR step |
!        %----------------------%
!
do 120 k = m, i - 1
!
!           ------------------------------------------------------------
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix. NR is the order of G.
!           ------------------------------------------------------------
!
nr = min( dlqthree, i-k+1 )
if( k.gt.m ) &
& call dcopy( nr, h( k, k-1 ), 1, v(1), 1 )
call dlarfg( nr, v( 1 ), v( 2 ), 1, t1 )
if( k.gt.m ) then
h( k, k-1 ) = v( 1 )
h( k+1, k-1 ) = zero
if( k.lt.i-1 ) &
& h( k+2, k-1 ) = zero
else if( m.gt.l ) then
h( k, k-1 ) = -h( k, k-1 )
end if
v2 = v( 2 )
t2 = t1*v2
if( nr.eq.3 ) then
v3 = v( 3 )
t3 = t1*v3
!
!              %------------------------------------------------%
!              | Apply G from the left to transform the rows of |
!              | the matrix in columns K to I2.                 |
!              %------------------------------------------------%
!
do 60 j = k, i2
sum = h( k, j ) + v2*h( k+1, j ) + v3*h( k+2, j )
h( k, j ) = h( k, j ) - sum*t1
h( k+1, j ) = h( k+1, j ) - sum*t2
h( k+2, j ) = h( k+2, j ) - sum*t3
60 continue
!
!              %----------------------------------------------------%
!              | Apply G from the right to transform the columns of |
!              | the matrix in rows I1 to min(K+3,I).               |
!              %----------------------------------------------------%
!
do 70 j = i1, min( k+3, i )
sum = h( j, k ) + v2*h( j, k+1 ) + v3*h( j, k+2 )
h( j, k ) = h( j, k ) - sum*t1
h( j, k+1 ) = h( j, k+1 ) - sum*t2
h( j, k+2 ) = h( j, k+2 ) - sum*t3
70 continue
!
!              %----------------------------------%
!              | Accumulate transformations for Z |
!              %----------------------------------%
!
sum      = z( k ) + v2*z( k+1 ) + v3*z( k+2 )
z( k )   = z( k ) - sum*t1
z( k+1 ) = z( k+1 ) - sum*t2
z( k+2 ) = z( k+2 ) - sum*t3

else if( nr.eq.2 ) then
!
!              %------------------------------------------------%
!              | Apply G from the left to transform the rows of |
!              | the matrix in columns K to I2.                 |
!              %------------------------------------------------%
!
do 90 j = k, i2
sum = h( k, j ) + v2*h( k+1, j )
h( k, j ) = h( k, j ) - sum*t1
h( k+1, j ) = h( k+1, j ) - sum*t2
90 continue
!
!              %----------------------------------------------------%
!              | Apply G from the right to transform the columns of |
!              | the matrix in rows I1 to min(K+3,I).               |
!              %----------------------------------------------------%
!
do 100 j = i1, i
sum = h( j, k ) + v2*h( j, k+1 )
h( j, k ) = h( j, k ) - sum*t1
h( j, k+1 ) = h( j, k+1 ) - sum*t2
100 continue
!
!              %----------------------------------%
!              | Accumulate transformations for Z |
!              %----------------------------------%
!
sum      = z( k ) + v2*z( k+1 )
z( k )   = z( k ) - sum*t1
z( k+1 ) = z( k+1 ) - sum*t2
end if
120 continue

130 continue
!
!     %-------------------------------------------------------%
!     | Failure to converge in remaining number of iterations |
!     %-------------------------------------------------------%
!
info = i
return

140 continue

if( l.eq.i ) then
!
!        %------------------------------------------------------%
!        | H(I,I-1) is negligible: one eigenvalue has converged |
!        %------------------------------------------------------%
!
wr( i ) = h( i, i )
wi( i ) = zero

else if( l.eq.i-1 ) then
!
!        %--------------------------------------------------------%
!        | H(I-1,I-2) is negligible;                              |
!        | a pair of eigenvalues have converged.                  |
!        |                                                        |
!        | Transform the 2-by-2 submatrix to standard Schur form, |
!        | and compute and store the eigenvalues.                 |
!        %--------------------------------------------------------%
!
call dlanv2( h( i-1, i-1 ), h( i-1, i ), h( i, i-1 ), &
& h( i, i ), wr( i-1 ), wi( i-1 ), wr( i ), wi( i ), &
& cs, sn )

if( wantt ) then
!
!           %-----------------------------------------------------%
!           | Apply the transformation to the rest of H and to Z, |
!           | as required.                                        |
!           %-----------------------------------------------------%
!
if( i2.gt.i ) &
& call drot( i2-i, h( i-1, i+1 ), ldh, h( i, i+1 ), ldh, &
& cs, sn )
call drot( i-i1-1, h( i1, i-1 ), ione, h( i1, i ), ione, &
& cs, sn )
sum      = cs*z( i-1 ) + sn*z( i )
z( i )   = cs*z( i )   - sn*z( i-1 )
z( i-1 ) = sum
end if
end if
!
!     %---------------------------------------------------------%
!     | Decrement number of remaining iterations, and return to |
!     | start of the main loop with new value of I.             |
!     %---------------------------------------------------------%
!
itn = itn - its
i = l - 1
go to 10

150 continue
return
!
!     %---------------%
!     | End of dlaqrb |
!     %---------------%
!
end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dngets
!
!\Description:
!  Given the eigenvalues of the upper Hessenberg matrix H,
!  computes the NP shifts AMU that are zeros of the polynomial of
!  degree NP which filters out components of the unwanted eigenvectors
!  corresponding to the AMU's based on some given criteria.
!
!  NOTE: call this even in the case of user specified shifts in order
!  to sort the eigenvalues, and error bounds of H for later use.
!
!\Usage:
!  call dngets
!     ( ISHIFT, WHICH, KEV, NP, RITZR, RITZI, BOUNDS, SHIFTR, SHIFTI )
!
!\Arguments
!  ISHIFT  Integer.  (INPUT)
!          Method for selecting the implicit shifts at each iteration.
!          ISHIFT = 0: user specified shifts
!          ISHIFT = 1: exact shift with respect to the matrix H.
!
!  WHICH   Character*2.  (INPUT)
!          Shift selection criteria.
!          'LM' -> want the KEV eigenvalues of largest magnitude.
!          'SM' -> want the KEV eigenvalues of smallest magnitude.
!          'LR' -> want the KEV eigenvalues of largest real part.
!          'SR' -> want the KEV eigenvalues of smallest real part.
!          'LI' -> want the KEV eigenvalues of largest imaginary part.
!          'SI' -> want the KEV eigenvalues of smallest imaginary part.
!
!  KEV      Integer.  (INPUT/OUTPUT)
!           INPUT: KEV+NP is the size of the matrix H.
!           OUTPUT: Possibly increases KEV by one to keep complex conjugate
!           pairs together.
!
!  NP       Integer.  (INPUT/OUTPUT)
!           Number of implicit shifts to be computed.
!           OUTPUT: Possibly decreases NP by one to keep complex conjugate
!           pairs together.
!
!  RITZR,  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
!  RITZI   On INPUT, RITZR and RITZI contain the real and imaginary
!          parts of the eigenvalues of H.
!          On OUTPUT, RITZR and RITZI are sorted so that the unwanted
!          eigenvalues are in the first NP locations and the wanted
!          portion is in the last KEV locations.  When exact shifts are
!          selected, the unwanted part corresponds to the shifts to
!          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
!          are further sorted so that the ones with largest Ritz values
!          are first.
!
!  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
!          Error bounds corresponding to the ordering in RITZ.
!
!  SHIFTR, SHIFTI  *** USE deprecated as of version 2.1. ***
!
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
!\Routines called:
!     dsortc  ARPACK sorting routine.
!     dcopy   Level 1 BLAS that copies one vector to another .
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.1'
!
!\SCCS Information: @(#)
! FILE: ngets.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
!
!\Remarks
!     1. xxxx
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dngets ( ishift, which, kev, np, ritzr, ritzi, bounds)
!     &                    shiftr, shifti )
!
implicit none
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
character  which*2
integer    ishift, kev, np
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
Double precision &
& bounds(kev+np), ritzr(kev+np), ritzi(kev+np)
!
!     %------------%
!     | Parameters |
!     %------------%
!
Double precision &
& zero
parameter (zero = 0.0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
logical    dngetstrue
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
external   dcopy, dsortc
!
!     %----------------------%
!     | Intrinsics Functions |
!     %----------------------%
!
intrinsic  abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
dngetstrue = .true.
!
!
!     %----------------------------------------------------%
!     | LM, SM, LR, SR, LI, SI case.                       |
!     | Sort the eigenvalues of H into the desired order   |
!     | and apply the resulting order to BOUNDS.           |
!     | The eigenvalues are sorted so that the wanted part |
!     | are always in the last KEV locations.              |
!     | We first do a pre-processing sort in order to keep |
!     | complex conjugate pairs together                   |
!     %----------------------------------------------------%
!
if (which .eq. 'LM') then
call dsortc ('LR', dngetstrue, kev+np, ritzr, ritzi, bounds)
else if (which .eq. 'SM') then
call dsortc ('SR', dngetstrue, kev+np, ritzr, ritzi, bounds)
else if (which .eq. 'LR') then
call dsortc ('LM', dngetstrue, kev+np, ritzr, ritzi, bounds)
else if (which .eq. 'SR') then
call dsortc ('SM', dngetstrue, kev+np, ritzr, ritzi, bounds)
else if (which .eq. 'LI') then
call dsortc ('LM', dngetstrue, kev+np, ritzr, ritzi, bounds)
else if (which .eq. 'SI') then
call dsortc ('SM', dngetstrue, kev+np, ritzr, ritzi, bounds)
end if
!
call dsortc (which, dngetstrue, kev+np, ritzr, ritzi, bounds)
!
!     %-------------------------------------------------------%
!     | Increase KEV by one if the ( ritzr(np),ritzi(np) )    |
!     | = ( ritzr(np+1),-ritzi(np+1) ) and ritz(np) .ne. zero |
!     | Accordingly decrease NP by one. In other words keep   |
!     | complex conjugate pairs together.                     |
!     %-------------------------------------------------------%
!
if (       (  abs(ritzr(np+1) - ritzr(np)) .le. zero) &
& .and. (  abs(ritzi(np+1) + ritzi(np)) .le. zero) ) then

np = np - 1
kev = kev + 1
end if
!
if ( ishift .eq. 1 ) then
!
!        %-------------------------------------------------------%
!        | Sort the unwanted Ritz values used as shifts so that  |
!        | the ones with largest Ritz estimates are first        |
!        | This will tend to minimize the effects of the         |
!        | forward instability of the iteration when they shifts |
!        | are applied in subroutine dnapps.                     |
!        | Be careful and use 'SR' since we want to sort BOUNDS! |
!        %-------------------------------------------------------%
!
call dsortc ( 'SR', dngetstrue, np, bounds, ritzr, ritzi )
end if
!
return
!
!     %---------------%
!     | End of dngets |
!     %---------------%
!
end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dnaitr
!
!\Description:
!  Reverse communication interface for applying NP additional steps to
!  a K step nonsymmetric Arnoldi factorization.
!
!  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
!
!          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
!
!  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
!
!          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
!
!  where OP and B are as in dnaupd.  The B-norm of r_{k+p} is also
!  computed and returned.
!
!\Usage:
!  call dnaitr
!     ( IDO, BMAT, N, K, NP, NB, RESID, RNORM, V, LDV, H, LDH,
!       IPNTR, WORKD, INFO )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y.
!                    This is for the restart phase to force the new
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y,
!                    IPNTR(3) is the pointer into WORK for B * X.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y.
!          IDO = 99: done
!          -------------------------------------------------------------
!          When the routine is used in the "shift-and-invert" mode, the
!          vector B * Q is already available and do not need to be
!          recompute in forming OP * Q.
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.  See dnaupd.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*M**x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  K       Integer.  (INPUT)
!          Current size of V and H.
!
!  NP      Integer.  (INPUT)
!          Number of additional Arnoldi steps to take.
!
!  NB      Integer.  (INPUT)
!          Blocksize to be used in the recurrence.
!          Only work for NB = 1 right now.  The goal is to have a
!          program that implement both the block and non-block method.
!
!  RESID   Double precision array of length N.  (INPUT/OUTPUT)
!          On INPUT:  RESID contains the residual vector r_{k}.
!          On OUTPUT: RESID contains the residual vector r_{k+p}.
!
!  RNORM   Double precision scalar.  (INPUT/OUTPUT)
!          B-norm of the starting residual on input.
!          B-norm of the updated residual r_{k+p} on output.
!
!  V       Double precision N by K+NP array.  (INPUT/OUTPUT)
!          On INPUT:  V contains the Arnoldi vectors in the first K
!          columns.
!          On OUTPUT: V contains the new NP Arnoldi vectors in the next
!          NP columns.  The first K columns are unchanged.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Double precision (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
!          H is used to store the generated upper Hessenberg matrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORK for
!          vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in the
!                    shift-and-invert mode.  X is the current operand.
!          -------------------------------------------------------------
!
!  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The calling program should not
!          use WORKD as temporary workspace during the iteration !!!!!!
!          On input, WORKD(1:N) = B*RESID and is used to save some
!          computation at the first step.
!
!  INFO    Integer.  (OUTPUT)
!          = 0: Normal exit.
!          > 0: Size of the spanning invariant subspace of OP found.
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
!     dgetv0  ARPACK routine to generate the initial vector.
!     ivout   ARPACK utility routine that prints integers.
!     second  ARPACK utility routine for timing.
!     dmout   ARPACK utility routine that prints matrices
!     dvout   ARPACK utility routine that prints vectors.
!     dlabad  LAPACK routine that computes machine constants.
!     dlamch  LAPACK routine that determines machine constants.
!     dlascl  LAPACK routine for careful scaling of a matrix.
!     dlanhs  LAPACK routine that computes various norms of a matrix.
!     dgemv   Level 2 BLAS routine for matrix vector multiplication.
!     daxpy   Level 1 BLAS that computes a vector triad.
!     dscal   Level 1 BLAS that scales a vector.
!     dcopy   Level 1 BLAS that copies one vector to another .
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
!\Revision history:
!     xx/xx/92: Version ' 2.4'
!
!\SCCS Information: @(#)
! FILE: naitr.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
!
!\Remarks
!  The algorithm implemented is:
!
!  restart = .false.
!  Given V_{k} = [v_{1}, ..., v_{k}], r_{k};
!  r_{k} contains the initial residual vector even for k = 0;
!  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already
!  computed by the calling program.
!
!  betaj = rnorm ; p_{k+1} = B*r_{k} ;
!  For  j = k+1, ..., k+np  Do
!     1) if ( betaj < tol ) stop or restart depending on j.
!        ( At present tol is zero )
!        if ( restart ) generate a new starting vector.
!     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];
!        p_{j} = p_{j}/betaj
!     3) r_{j} = OP*v_{j} where OP is defined as in dnaupd
!        For shift-invert mode p_{j} = B*v_{j} is already available.
!        wnorm = || OP*v_{j} ||
!     4) Compute the j-th step residual vector.
!        w_{j} =  V_{j}^T * B * OP * v_{j}
!        r_{j} =  OP*v_{j} - V_{j} * w_{j}
!        H(:,j) = w_{j};
!        H(j,j-1) = rnorm
!        rnorm = || r_(j) ||
!        If (rnorm > 0.717*wnorm) accept step and go back to 1)
!     5) Re-orthogonalization step:
!        s = V_{j}'*B*r_{j}
!        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
!        alphaj = alphaj + s_{j};
!     6) Iterative refinement step:
!        If (rnorm1 > 0.717*rnorm) then
!           rnorm = rnorm1
!           accept step and go back to 1)
!        Else
!           rnorm = rnorm1
!           If this is the first time in step 6), go to 5)
!           Else r_{j} lies in the span of V_{j} numerically.
!              Set r_{j} = 0 and rnorm = 0; go to 1)
!        EndIf
!  End Do
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dnaitr &
& (ido, bmat, n, k, np, nb, resid, rnorm, v, ldv, h, ldh, &
& ipntr, workd, info)

implicit none
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
character  bmat*1
integer    ido, info, k, ldh, ldv, n, nb, np
Double precision &
& rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
integer    ipntr(3)
Double precision &
& h(ldh,k+np), resid(n), v(ldv,k+np), workd(3*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
Double precision &
& one, zero
parameter (one = 1.0D+0, zero = 0.0D+0)
logical    fls
parameter (fls = .false. )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
logical    first, orth1, orth2, rstart, step3, step4
integer    ierr, i, infol, ipj, irj, ivj, iter, itry, j, &
& jj
integer    dnaitrone
parameter (dnaitrone = 1)
!
integer(4) &
& ione
parameter (ione = 1)
!
Double precision &
& betaj, ovfl, temp1, rnorm1, smlnum, tst1, ulp, unfl, &
& wnorm
save       first, orth1, orth2, rstart, step3, step4, &
& ierr, ipj, irj, ivj, iter, itry, j,  ovfl, &
& betaj, rnorm1, smlnum, ulp, unfl, wnorm
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
external   daxpy, dcopy, dscal, dgemv, dgetv0, dlabad
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
Double precision &
& ddot, dnrm2, dlanhs, dlamch
external   ddot, dnrm2, dlanhs, dlamch
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
intrinsic    abs, sqrt
!
!     %-----------------%
!     | Data statements |
!     %-----------------%
!
data      first / .true. /
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
if (nb .gt. 1000) then
goto 9000
end if
!
if (first) then
!
!        %-----------------------------------------%
!        | Set machine-dependent constants for the |
!        | the splitting and deflation criterion.  |
!        | If norm(H) <= sqrt(OVFL),               |
!        | overflow should not occur.              |
!        | REFERENCE: LAPACK subroutine dlahqr     |
!        %-----------------------------------------%
!
unfl = dlamch( 'safe minimum' )
ovfl = one / unfl
call dlabad( unfl, ovfl )
ulp = dlamch( 'precision' )
smlnum = unfl*( n / ulp )
first = .false.
end if
!
if (ido .eq. 0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
!         call second (t0)
!m         msglvl = mnaitr
!
!        %------------------------------%
!        | Initial call to this routine |
!        %------------------------------%
!
info   = 0
step3  = .false.
step4  = .false.
rstart = .false.
orth1  = .false.
orth2  = .false.
j      = k + 1
ipj    = 1
irj    = ipj   + n
ivj    = irj   + n
end if
!
!     %-------------------------------------------------%
!     | When in reverse communication mode one of:      |
!     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
!     | will be .true. when ....                        |
!     | STEP3: return from computing OP*v_{j}.          |
!     | STEP4: return from computing B-norm of OP*v_{j} |
!     | ORTH1: return from computing B-norm of r_{j+1}  |
!     | ORTH2: return from computing B-norm of          |
!     |        correction to the residual vector.       |
!     | RSTART: return from OP computations needed by   |
!     |         dgetv0.                                 |
!     %-------------------------------------------------%
!
if (step3)  go to 50
if (step4)  go to 60
if (orth1)  go to 70
if (orth2)  go to 90
if (rstart) go to 30
!
!     %-----------------------------%
!     | Else this is the first step |
!     %-----------------------------%
!
!     %--------------------------------------------------------------%
!     |                                                              |
!     |        A R N O L D I     I T E R A T I O N     L O O P       |
!     |                                                              |
!     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
!     %--------------------------------------------------------------%

1000 continue
!
!        %---------------------------------------------------%
!        | STEP 1: Check if the B norm of j-th residual      |
!        | vector is zero. Equivalent to determing whether   |
!        | an exact j-step Arnoldi factorization is present. |
!        %---------------------------------------------------%
!
betaj = rnorm
if (rnorm .gt. zero) go to 40
!
!           %---------------------------------------------------%
!           | Invariant subspace found, generate a new starting |
!           | vector which is orthogonal to the current Arnoldi |
!           | basis and continue the iteration.                 |
!           %---------------------------------------------------%
!           %---------------------------------------------%
!           | ITRY is the loop variable that controls the |
!           | maximum amount of times that a restart is   |
!           | attempted. NRSTRT is used by stat.h         |
!           %---------------------------------------------%
!
betaj  = zero
!p            nrstrt = nrstrt + 1
itry   = 1
20 continue
rstart = .true.
ido    = 0
30 continue
!
!           %--------------------------------------%
!           | If in reverse communication mode and |
!           | RSTART = .true. flow returns here.   |
!           %--------------------------------------%
!
call dgetv0 (ido, bmat, fls, n, j, v, ldv, &
& resid, rnorm, ipntr, workd, ierr)
if (ido .ne. 99) go to 9000
if (ierr .lt. 0) then
itry = itry + 1
if (itry .le. 3) go to 20
!
!              %------------------------------------------------%
!              | Give up after several restart attempts.        |
!              | Set INFO to the size of the invariant subspace |
!              | which spans OP and exit.                       |
!              %------------------------------------------------%
!
info = j - 1
ido = 99
go to 9000
end if
!
40 continue
!
!        %---------------------------------------------------------%
!        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
!        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
!        | when reciprocating a small RNORM, test against lower    |
!        | machine bound.                                          |
!        %---------------------------------------------------------%
!
call dcopy (n, resid, 1, v(1,j), 1)
if (rnorm .ge. unfl) then
temp1 = one / rnorm
call dscal (n, temp1, v(1,j), ione)
call dscal (n, temp1, workd(ipj), ione)
else
!
!            %-----------------------------------------%
!            | To scale both v_{j} and p_{j} carefully |
!            | use LAPACK routine SLASCL               |
!            %-----------------------------------------%
!
call dlascl ('General', i, i, rnorm, one, n, 1, &
& v(1,j), n, infol)
call dlascl ('General', i, i, rnorm, one, n, 1, &
& workd(ipj), n, infol)
end if
!
!        %------------------------------------------------------%
!        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
!        | Note that this is not quite yet r_{j}. See STEP 4    |
!        %------------------------------------------------------%
!
step3 = .true.
call dcopy (n, v(1,j), 1, workd(ivj), 1)
ipntr(1) = ivj
ipntr(2) = irj
ipntr(3) = ipj
ido = 1
!
!        %-----------------------------------%
!        | Exit in order to compute OP*v_{j} |
!        %-----------------------------------%
!
go to 9000
50 continue
!
!        %----------------------------------%
!        | Back from reverse communication; |
!        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   |
!        | if step3 = .true.                |
!        %----------------------------------%
!
step3 = .false.
!
!        %------------------------------------------%
!        | Put another copy of OP*v_{j} into RESID. |
!        %------------------------------------------%
!
call dcopy (n, workd(irj), 1, resid(1), 1)
!
!        %---------------------------------------%
!        | STEP 4:  Finish extending the Arnoldi |
!        |          factorization to length j.   |
!        %---------------------------------------%
!
if (bmat .eq. 'G') then
step4 = .true.
ipntr(1) = irj
ipntr(2) = ipj
ido = 2
!
!           %-------------------------------------%
!           | Exit in order to compute B*OP*v_{j} |
!           %-------------------------------------%
!
go to 9000
else if (bmat .eq. 'I') then
call dcopy (n, resid, 1, workd(ipj), 1)
end if
60 continue
!
!        %----------------------------------%
!        | Back from reverse communication; |
!        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} |
!        | if step4 = .true.                |
!        %----------------------------------%
!
step4 = .false.
!
!        %-------------------------------------%
!        | The following is needed for STEP 5. |
!        | Compute the B-norm of OP*v_{j}.     |
!        %-------------------------------------%
!
if (bmat .eq. 'G') then
wnorm = ddot (n, resid, 1, workd(ipj), 1)
wnorm = sqrt(abs(wnorm))
else if (bmat .eq. 'I') then
wnorm = dnrm2(n, resid, 1)
end if
!
!        %-----------------------------------------%
!        | Compute the j-th residual corresponding |
!        | to the j step factorization.            |
!        | Use Classical Gram Schmidt and compute: |
!        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
!        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
!        %-----------------------------------------%
!
!
!        %------------------------------------------%
!        | Compute the j Fourier coefficients w_{j} |
!        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
!        %------------------------------------------%
!
call dgemv ('T', n, j, one, v, ldv, workd(ipj), 1, &
& zero, h(1,j), 1)
!
!        %--------------------------------------%
!        | Orthogonalize r_{j} against V_{j}.   |
!        | RESID contains OP*v_{j}. See STEP 3. |
!        %--------------------------------------%
!
call dgemv ('N', n, j, -one, v, ldv, h(1,j), 1, &
& one, resid(1), 1)
!
if (j .gt. 1) h(j,j-1) = betaj
!
orth1 = .true.
!
if (bmat .eq. 'G') then
call dcopy (n, resid, 1, workd(irj), 1)
ipntr(1) = irj
ipntr(2) = ipj
ido = 2
!
!           %----------------------------------%
!           | Exit in order to compute B*r_{j} |
!           %----------------------------------%
!
go to 9000
else if (bmat .eq. 'I') then
call dcopy (n, resid, 1, workd(ipj), 1)
end if
70 continue
!
!        %---------------------------------------------------%
!        | Back from reverse communication if ORTH1 = .true. |
!        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
!        %---------------------------------------------------%
!
orth1 = .false.
!
!        %------------------------------%
!        | Compute the B-norm of r_{j}. |
!        %------------------------------%
!
if (bmat .eq. 'G') then
rnorm = ddot (n, resid, 1, workd(ipj), 1)
rnorm = sqrt(abs(rnorm))
else if (bmat .eq. 'I') then
rnorm = dnrm2(n, resid, 1)
end if
!
!        %-----------------------------------------------------------%
!        | STEP 5: Re-orthogonalization / Iterative refinement phase |
!        | Maximum NITER_ITREF tries.                                |
!        |                                                           |
!        |          s      = V_{j}^T * B * r_{j}                     |
!        |          r_{j}  = r_{j} - V_{j}*s                         |
!        |          alphaj = alphaj + s_{j}                          |
!        |                                                           |
!        | The stopping criteria used for iterative refinement is    |
!        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
!        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
!        | Determine if we need to correct the residual. The goal is |
!        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
!        | The following test determines whether the sine of the     |
!        | angle between  OP*x and the computed residual is less     |
!        | than or equal to 0.717.                                   |
!        %-----------------------------------------------------------%
!
if (rnorm .gt. 0.717*wnorm) go to 100
iter  = 0
!
!        %---------------------------------------------------%
!        | Enter the Iterative refinement phase. If further  |
!        | refinement is necessary, loop back here. The loop |
!        | variable is ITER. Perform a step of Classical     |
!        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
!        %---------------------------------------------------%
!
80 continue
!
!        %----------------------------------------------------%
!        | Compute V_{j}^T * B * r_{j}.                       |
!        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
!        %----------------------------------------------------%
!
call dgemv ('T', n, j, one, v, ldv, workd(ipj), 1, &
& zero, workd(irj), 1)
!
!        %---------------------------------------------%
!        | Compute the correction to the residual:     |
!        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
!        | The correction to H is v(:,1:J)*H(1:J,1:J)  |
!        | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
!        %---------------------------------------------%
!
call dgemv ('N', n, j, -one, v, ldv, workd(irj), 1, &
& one, resid(1), 1)
call daxpy (j, one, workd(irj), 1, h(1,j), 1)
!
orth2 = .true.
if (bmat .eq. 'G') then
call dcopy (n, resid, 1, workd(irj), 1)
ipntr(1) = irj
ipntr(2) = ipj
ido = 2
!
!           %-----------------------------------%
!           | Exit in order to compute B*r_{j}. |
!           | r_{j} is the corrected residual.  |
!           %-----------------------------------%
!
go to 9000
else if (bmat .eq. 'I') then
call dcopy (n, resid, 1, workd(ipj), 1)
end if
90 continue
!
!        %---------------------------------------------------%
!        | Back from reverse communication if ORTH2 = .true. |
!        %---------------------------------------------------%
!
!        %-----------------------------------------------------%
!        | Compute the B-norm of the corrected residual r_{j}. |
!        %-----------------------------------------------------%
!
if (bmat .eq. 'G') then
rnorm1 = ddot (n, resid, 1, workd(ipj), 1)
rnorm1 = sqrt(abs(rnorm1))
else if (bmat .eq. 'I') then
rnorm1 = dnrm2(n, resid, 1)
end if
!
!        %-----------------------------------------%
!        | Determine if we need to perform another |
!        | step of re-orthogonalization.           |
!        %-----------------------------------------%
!
if (rnorm1 .gt. 0.717*rnorm) then
!
!           %---------------------------------------%
!           | No need for further refinement.       |
!           | The cosine of the angle between the   |
!           | corrected residual vector and the old |
!           | residual vector is greater than 0.717 |
!           | In other words the corrected residual |
!           | and the old residual vector share an  |
!           | angle of less than arcCOS(0.717)      |
!           %---------------------------------------%
!
rnorm = rnorm1
!
else
!
!           %-------------------------------------------%
!           | Another step of iterative refinement step |
!           | is required. NITREF is used by stat.h     |
!           %-------------------------------------------%
!
rnorm  = rnorm1
iter   = iter + 1
if (iter .le. 1) go to 80
!
!           %-------------------------------------------------%
!           | Otherwise RESID is numerically in the span of V |
!           %-------------------------------------------------%
!
do 95 jj = 1, n
resid(jj) = zero
95 continue
rnorm = zero
end if
!
!        %----------------------------------------------%
!        | Branch here directly if iterative refinement |
!        | wasn't necessary or after at most NITER_REF  |
!        | steps of iterative refinement.               |
!        %----------------------------------------------%
!
100 continue
!
rstart = .false.
orth2  = .false.
!
!        %------------------------------------%
!        | STEP 6: Update  j = j+1;  Continue |
!        %------------------------------------%
!
j = j + 1
if (j .gt. k+np) then
ido = 99
do 110 i = max(dnaitrone,k), k+np-1
!
!              %--------------------------------------------%
!              | Check for splitting and deflation.         |
!              | Use a standard test as in the QR algorithm |
!              | REFERENCE: LAPACK subroutine dlahqr        |
!              %--------------------------------------------%
!
tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
if( abs(tst1) .le. zero ) &
& tst1 = dlanhs( '1', k+np, h, ldh, workd(n+1) )
if( abs( h( i+1,i ) ).le.max( ulp*tst1, smlnum ) ) &
& h(i+1,i) = zero
110 continue
!
go to 9000
end if
!
!        %--------------------------------------------------------%
!        | Loop back to extend the factorization by another step. |
!        %--------------------------------------------------------%
!
go to 1000
!
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
!
9000 continue
return
!
!     %---------------%
!     | End of dnaitr |
!     %---------------%
!
end
!
