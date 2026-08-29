c     Modifications:
c     2023-10-17: eliminated all kind=4 occurences.
c
c      
c dneupd.f
c dnaupd.f
c dnaup2.f
c dnapps.f
c dnconv.f
c dsortc.f
c dneigh.f
c dlaqrb.f
c dngets.f
c dgetv0.f
c dnaitr.f
c


c\BeginDoc
c
c\Name: dneupd
c
c\Description:
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  basis is always computed.  There is an additional storage cost of n*nev
c  if both are requested (in this case a separate array Z must be supplied).
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are derived from approximate eigenvalues and eigenvectors of
c  of the linear operator OP prescribed by the MODE selection in the
c  call to DNAUPD .  DNAUPD  must be called before this routine is called.
c  These approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such
c  in the comments that follow.  The computed orthonormal basis for the
c  invariant subspace corresponding to these Ritz values is referred to as a
c  Schur basis.
c
c  See documentation in the header of the subroutine DNAUPD  for
c  definition of OP as well as other terms and the relation of computed
c  Ritz values and Ritz vectors of OP with respect to the given problem
c  A*z = lambda*B*z.  For a brief description, see definitions of
c  IPARAM(7), MODE and WHICH in the documentation of DNAUPD .
c
c\Usage:
c  call dneupd
c     ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, WORKEV, BMAT,
c       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL,
c       LWORKL, INFO )
c
c\Arguments:
c  RVEC    LOGICAL  (INPUT)
c          Specifies whether a basis for the invariant subspace corresponding
c          to the converged Ritz value approximations for the eigenproblem
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
c                                See Remarks below.
c
c  HOWMNY  Character*1  (INPUT)
c          Specifies the form of the basis for the invariant subspace
c          corresponding to the converged Ritz values that is to be computed.
c
c          = 'A': Compute NEV Ritz vectors;
c          = 'P': Compute NEV Schur vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE..
c          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
c
c  DR      Double precision  array of dimension NEV+1.  (OUTPUT)
c          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains
c          the real part of the Ritz  approximations to the eigenvalues of
c          A*z = lambda*B*z.
c          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit:
c          DR contains the real part of the Ritz values of OP computed by
c          DNAUPD . A further computation must be performed by the user
c          to transform the Ritz values computed for OP by DNAUPD  to those
c          of the original system A*z = lambda*B*z. See remark 3 below.
c
c  DI      Double precision  array of dimension NEV+1.  (OUTPUT)
c          On exit, DI contains the imaginary part of the Ritz value
c          approximations to the eigenvalues of A*z = lambda*B*z associated
c          with DR.
c
c          NOTE: When Ritz values are complex, they will come in complex
c                conjugate pairs.  If eigenvectors are requested, the
c                corresponding Ritz vectors will also come in conjugate
c                pairs and the real and imaginary parts of these are
c                represented in two consecutive columns of the array Z
c                (see below).
c
c  Z       Double precision  N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
c          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
c          Z represent approximate eigenvectors (Ritz vectors) corresponding
c          to the NCONV=IPARAM(5) Ritz values for eigensystem
c          A*z = lambda*B*z.
c
c          The complex Ritz vector associated with the Ritz value
c          with positive imaginary part is stored in two consecutive
c          columns.  The first column holds the real part of the Ritz
c          vector and the second column holds the imaginary part.  The
c          Ritz vector associated with the Ritz value with negative
c          imaginary part is simply the complex conjugate of the Ritz vector
c          associated with the positive imaginary part.
c
c          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi
c          basis array V computed by DNAUPD .  In this case the Arnoldi basis
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
c
c  SIGMAR  Double precision   (INPUT)
c          If IPARAM(7) = 3 or 4, represents the real part of the shift.
c          Not referenced if IPARAM(7) = 1 or 2.
c
c  SIGMAI  Double precision   (INPUT)
c          If IPARAM(7) = 3 or 4, represents the imaginary part of the shift.
c          Not referenced if IPARAM(7) = 1 or 2. See remark 3 below.
c
c  WORKEV  Double precision  work array of dimension 3*NCV.  (WORKSPACE)
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to DNAUPD  that was just completed.               ****
c
c  NOTE: The remaining arguments
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
c           WORKD, WORKL, LWORKL, INFO
c
c         must be passed directly to DNEUPD  following the last call
c         to DNAUPD .  These arguments MUST NOT BE MODIFIED between
c         the the last call to DNAUPD  and the call to DNEUPD .
c
c  Three of these parameters (V, WORKL, INFO) are also output parameters:
c
c  V       Double precision  N by NCV array.  (INPUT/OUTPUT)
c
c          Upon INPUT: the NCV columns of V contain the Arnoldi basis
c                      vectors for OP as constructed by DNAUPD  .
c
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       contain approximate Schur vectors that span the
c                       desired invariant subspace.  See Remark 2 below.
c
c          NOTE: If the array Z has been set equal to first NEV+1 columns
c          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
c          Arnoldi basis held by V has been overwritten by the desired
c          Ritz vectors.  If a separate array Z has been passed then
c          the first NCONV=IPARAM(5) columns of V will contain approximate
c          Schur vectors that span the desired invariant subspace.
c
c  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:ncv*ncv+3*ncv) contains information obtained in
c          dnaupd .  They are not changed by dneupd .
c          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the
c          real and imaginary part of the untransformed Ritz values,
c          the upper quasi-triangular matrix for H, and the
c          associated matrix representation of the invariant subspace for H.
c
c          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
c          of the above information computed by dneupd .
c          -------------------------------------------------------------
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     dneupd  if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c
c          =  0: Normal exit.
c
c          =  1: The Schur form computed by LAPACK routine dlahqr
c                could not be reordered by LAPACK routine dtrsen .
c                Re-enter subroutine dneupd  with IPARAM(5)=NCV and
c                increase the size of the arrays DR and DI to have
c                dimension at least dimension NCV and allocate at least NCV
c                columns for Z. NOTE: Not necessary if Z and V share
c                the same space. Please notify the authors if this error
c                occurs.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from calculation of a real Schur form.
c                Informational error from LAPACK routine dlahqr .
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine dtrevc .
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
c          = -14: DNAUPD  did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: DNEUPD got a different count of the number of converged
c                 Ritz values than DNAUPD got.  This indicates the user
c                 probably made an error in passing data from DNAUPD to
c                 DNEUPD or that the data was modified before entering
c                 DNEUPD
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
c     Real Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     dmout    ARPACK utility routine that prints matrices
c     dvout    ARPACK utility routine that prints vectors.
c     dgeqr2   LAPACK routine that computes the QR factorization of
c             a matrix.
c     dlacpy   LAPACK matrix copy routine.
c     dlahqr   LAPACK routine to compute the real Schur form of an
c             upper Hessenberg matrix.
c     dlamch   LAPACK routine that determines machine constants.
c     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlaset   LAPACK matrix initialization routine.
c     dorm2r   LAPACK routine that applies an orthogonal matrix in
c             factored form.
c     dtrevc   LAPACK routine to compute the eigenvectors of a matrix
c             in upper quasi-triangular form.
c     dtrsen   LAPACK routine that re-orders the Schur form.
c     dtrmm    Level 3 BLAS matrix times an upper triangular matrix.
c     dger     Level 2 BLAS rank one update to a matrix.
c     dcopy    Level 1 BLAS that copies one vector to another .
c     ddot     Level 1 BLAS that computes the scalar product of two vectors.
c     dnrm2    Level 1 BLAS that computes the norm of a vector.
c     dscal    Level 1 BLAS that scales a vector.
c
c\Remarks
c
c  1. Currently only HOWMNY = 'A' and 'P' are implemented.
c
c     Let trans(X) denote the transpose of X.
c
c  2. Schur vectors are an orthogonal representation for the basis of
c     Ritz vectors. Thus, their numerical properties are often superior.
c     If RVEC = .TRUE. then the relationship
c             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
c     trans(V(:,1:IPARAM(5))) * V(:,1:IPARAM(5)) = I are approximately
c     satisfied. Here T is the leading submatrix of order IPARAM(5) of the
c     real upper quasi-triangular matrix stored workl(ipntr(12)). That is,
c     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
c     each 2-by-2 diagonal block has its diagonal elements equal and its
c     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
c     diagonal block is a complex conjugate pair of Ritz values. The real
c     Ritz values are stored on the diagonal of T.
c
c  3. If IPARAM(7) = 3 or 4 and SIGMAI is not equal zero, then the user must
c     form the IPARAM(5) Rayleigh quotients in order to transform the Ritz
c     values computed by DNAUPD  for OP to those of A*z = lambda*B*z.
c     Set RVEC = .true. and HOWMNY = 'A', and
c     compute
c           trans(Z(:,I)) * A * Z(:,I) if DI(I) = 0.
c     If DI(I) is not equal to zero and DI(I+1) = - D(I),
c     then the desired real and imaginary parts of the Ritz value are
c           trans(Z(:,I)) * A * Z(:,I) +  trans(Z(:,I+1)) * A * Z(:,I+1),
c           trans(Z(:,I)) * A * Z(:,I+1) -  trans(Z(:,I+1)) * A * Z(:,I),
c     respectively.
c     Another possibility is to set RVEC = .true. and HOWMNY = 'P' and
c     compute trans(V(:,1:IPARAM(5))) * A * V(:,1:IPARAM(5)) and then an upper
c     quasi-triangular matrix of order IPARAM(5) is computed. See remark
c     2 above.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Chao Yang                    Houston, Texas
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: neupd.F   SID: 2.7   DATE OF SID: 09/20/00   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
      subroutine dneupd (rvec , howmny, select, dr    , di,
     &                   z    , ldz   , sigmar, sigmai, workev,
     &                   bmat , n     , which , nev   , tol,
     &                   resid, ncv   , v     , ldv   , iparam,
     &                   ipntr, workd , workl , lworkl, info)
c
      implicit none
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision
     &           sigmar, sigmai, tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(8), ipntr(14)
      logical    select(ncv)
      Double precision
     &           dr(nev+1)    , di(nev+1), resid(n)  ,
     &           v(ldv,ncv)   , z(ldz,*) , workd(3*n),
     &           workl(lworkl), workev(3*ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0 , zero = 0.0D-5 )
c
      integer(4)
     &           ione
      parameter (ione = 1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  type*6
      integer    bounds, ierr  , ih    , ihbds   ,
     &           iheigr, iheigi, iconj , nconv   ,
     &           invsub, iuptri, iwev  , iwork(1),
     &           j     , k     , ldh   , ldq     ,
     &           mode  , outncv, ritzr   ,
     &           ritzi , wri   , wrr   , irr     ,
     &           iri   , ibd   , ishift, numcnv  ,
     &           np    , jj
      logical    reord
      Double precision
     &           conds  , rnorm, sep  , temp,
     &           vl(1,1), temp1, eps23
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy  , dger   , dgeqr2 , dlacpy ,
     &           dlahqr , dlaset , dorm2r ,
     &           dtrevc , dtrmm  , dtrsen , dscal
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlapy2 , dnrm2 , dlamch , ddot
      external   dlapy2 , dnrm2 , dlamch , ddot
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, min, sqrt
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
c
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch ('Epsilon-Machine')
      eps23 = eps23**(2.0D+0  / 3.0D+0 )
c
c     %--------------%
c     | Quick return |
c     %--------------%
c
      ierr = 0
c
      if (nconv .le. 0) then
         ierr = -14
      else if (n .le. 0) then
         ierr = -1
      else if (nev .le. 0) then
         ierr = -2
      else if (ncv .le. nev+1 .or.  ncv .gt. n) then
         ierr = -3
      else if (which .ne. 'LM' .and.
     &        which .ne. 'SM' .and.
     &        which .ne. 'LR' .and.
     &        which .ne. 'SR' .and.
     &        which .ne. 'LI' .and.
     &        which .ne. 'SI') then
         ierr = -5
      else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
         ierr = -6
      else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
         ierr = -7
      else if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec ) then
         ierr = -13
      else if (howmny .eq. 'S' ) then
         ierr = -12
      end if
c
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
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
c
c     %--------------------------------------------------------%
c     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q   |
c     | etc... and the remaining workspace.                    |
c     | Also update pointer to be used on output.              |
c     | Memory is laid out as follows:                         |
c     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
c     | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary   |
c     |                                   parts of ritz values |
c     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds   |
c     %--------------------------------------------------------%
c
c     %-----------------------------------------------------------%
c     | The following is used and set by DNEUPD .                  |
c     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
c     |                             real part of the Ritz values. |
c     | workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed |
c     |                        imaginary part of the Ritz values. |
c     | workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed |
c     |                           error bounds of the Ritz values |
c     | workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper |
c     |                             quasi-triangular matrix for H |
c     | workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the    |
c     |       associated matrix representation of the invariant   |
c     |       subspace for H.                                     |
c     | GRAND total of NCV * ( 3 * NCV + 6 ) locations.           |
c     %-----------------------------------------------------------%
c
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
c
c     %-----------------------------------------%
c     | irr points to the REAL part of the Ritz |
c     |     values computed by _neigh before    |
c     |     exiting _naup2.                     |
c     | iri points to the IMAGINARY part of the |
c     |     Ritz values computed by _neigh      |
c     |     before exiting _naup2.              |
c     | ibd points to the Ritz estimates        |
c     |     computed by _neigh before exiting   |
c     |     _naup2.                             |
c     %-----------------------------------------%
c
      irr = ipntr(14)+ncv*ncv
      iri = irr+ncv
      ibd = iri+ncv
c
c     %------------------------------------%
c     | RNORM is B-norm of the RESID(1:N). |
c     %------------------------------------%
c
      rnorm = workl(ih+2)
      workl(ih+2) = zero
c
c
      if (rvec) then
c
         reord = .false.
c
c        %---------------------------------------------------%
c        | Use the temporary bounds array to store indices   |
c        | These will be used to mark the select array later |
c        %---------------------------------------------------%
c
         do 10 j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
   10    continue
c
c        %-------------------------------------%
c        | Select the wanted Ritz values.      |
c        | Sort the Ritz values so that the    |
c        | wanted ones appear at the tailing   |
c        | NEV positions of workl(irr) and     |
c        | workl(iri).  Move the corresponding |
c        | error estimates in workl(bound)     |
c        | accordingly.                        |
c        %-------------------------------------%
c
         np     = ncv - nev
         ishift = 0
         call dngets (ishift       , which     , nev       ,
     &                np           , workl(irr), workl(iri),
     &                workl(bounds))
c
c
c        %-----------------------------------------------------%
c        | Record indices of the converged wanted Ritz values  |
c        | Mark the select array for possible reordering       |
c        %-----------------------------------------------------%
c
         numcnv = 0
         do 11 j = 1,ncv
            temp1 = max(eps23,
     &                 dlapy2 ( workl(irr+ncv-j), workl(iri+ncv-j) ))
            jj = INT(workl(bounds + ncv - j))
            if (numcnv .lt. nconv .and.
     &          workl(ibd+jj-1) .le. tol*temp1) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj .gt. nev) reord = .true.
            endif
   11    continue
c
c        %-----------------------------------------------------------%
c        | Check the count (numcnv) of converged Ritz values with    |
c        | the number (nconv) reported by dnaupd.  If these two      |
c        | are different then there has probably been an error       |
c        | caused by incorrect passing of the dnaupd data.           |
c        %-----------------------------------------------------------%
c
c
         if (numcnv .ne. nconv) then
            info = -15
            go to 9000
         end if
c
c        %-----------------------------------------------------------%
c        | Call LAPACK routine dlahqr  to compute the real Schur form |
c        | of the upper Hessenberg matrix returned by DNAUPD .        |
c        | Make a copy of the upper Hessenberg matrix.               |
c        | Initialize the Schur vector matrix Q to the identity.     |
c        %-----------------------------------------------------------%
c
         call dcopy (ldh*ncv, workl(ih), ione,
     &               workl(iuptri), ione)
         call dlaset ('All', ncv, ncv,
     &                zero , one, workl(invsub),
     &                ldq)
         call dlahqr (.true., .true.       , ncv,
     &                ione, ncv, workl(iuptri),
     &                ldh   , workl(iheigr), workl(iheigi),
     &                ione, ncv, workl(invsub),
     &                ldq   , ierr)
         call dcopy (ncv, workl(invsub+ncv-1), ldq,
     &               workl(ihbds), 1)
c
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
c
         if (reord) then
c
c           %-----------------------------------------------------%
c           | Reorder the computed upper quasi-triangular matrix. |
c           %-----------------------------------------------------%
c
            call dtrsen ('None'       , 'V'          ,
     &                   select       , ncv          ,
     &                   workl(iuptri), ldh          ,
     &                   workl(invsub), ldq          ,
     &                   workl(iheigr), workl(iheigi),
     &                   nconv        , conds        ,
     &                   sep          , workl(ihbds) ,
     &                   ncv          , iwork        ,
     &                   1            , ierr)
c
            if (ierr .eq. 1) then
               info = 1
               go to 9000
            end if
c
         end if
c
c        %---------------------------------------%
c        | Copy the last row of the Schur vector |
c        | into workl(ihbds).  This will be used |
c        | to compute the Ritz estimates of      |
c        | converged Ritz values.                |
c        %---------------------------------------%
c
         call dcopy (ncv, workl(invsub+ncv-1), ldq, 
     &               workl(ihbds), ione)
c
c        %----------------------------------------------------%
c        | Place the computed eigenvalues of H into DR and DI |
c        | if a spectral transformation was not used.         |
c        %----------------------------------------------------%
c
         if (type .eq. 'REGULR') then
            call dcopy (nconv, workl(iheigr), 1, dr(ione), 1)
            call dcopy (nconv, workl(iheigi), 1, di(ione), 1)
         end if
c
c        %----------------------------------------------------------%
c        | Compute the QR factorization of the matrix representing  |
c        | the wanted invariant subspace located in the first NCONV |
c        | columns of workl(invsub,ldq).                            |
c        %----------------------------------------------------------%
c
         call dgeqr2 (ncv, nconv , workl(invsub),
     &               ldq, workev, workev(ncv+1),
     &               ierr)
c
c        %---------------------------------------------------------%
c        | * Postmultiply V by Q using dorm2r .                     |
c        | * Copy the first NCONV columns of VQ into Z.            |
c        | * Postmultiply Z by R.                                  |
c        | The N by NCONV matrix Z is now a matrix representation  |
c        | of the approximate invariant subspace associated with   |
c        | the Ritz values in workl(iheigr) and workl(iheigi)      |
c        | The first NCONV columns of V are now approximate Schur  |
c        | vectors associated with the real upper quasi-triangular |
c        | matrix of order NCONV in workl(iuptri)                  |
c        %---------------------------------------------------------%
c
         call dorm2r ('Right', 'Notranspose', n            ,
     &                ncv   , nconv        , workl(invsub),
     &                ldq   , workev       , v            ,
     &                ldv   , workd(n+1)   , ierr)
         call dlacpy ('All', n, nconv, v(1,1), ldv, z, ldz)
c
         do 20 j=1, nconv
c
c           %---------------------------------------------------%
c           | Perform both a column and row scaling if the      |
c           | diagonal element of workl(invsub,ldq) is negative |
c           | I'm lazy and don't take advantage of the upper    |
c           | quasi-triangular form of workl(iuptri,ldq)        |
c           | Note that since Q is orthogonal, R is a diagonal  |
c           | matrix consisting of plus or minus ones           |
c           %---------------------------------------------------%
c
            if (workl(invsub+(j-1)*ldq+j-1) .lt. zero) then
               call dscal (nconv, -one, workl(iuptri+j-1),
     &                     ldq)
               call dscal (nconv, -one, workl(iuptri+(j-1)*ldq),
     &                     ione)
            end if
c
 20      continue
c
         if (howmny .eq. 'A') then
c
c           %--------------------------------------------%
c           | Compute the NCONV wanted eigenvectors of T |
c           | located in workl(iuptri,ldq).              |
c           %--------------------------------------------%
c
            do 30 j=1, ncv
               if (j .le. nconv) then
                  select(j) = .true.
               else
                  select(j) = .false.
               end if
 30         continue
c
            call dtrevc ('Right', 'Select'     , select       ,
     &                   ncv    , workl(iuptri), ldq          ,
     &                   vl, ione, workl(invsub),
     &                   ldq    , ncv          , outncv       ,
     &                   workev(1) , ierr)
c
            if (ierr .ne. 0) then
                info = -9
                go to 9000
            end if
c
c           %------------------------------------------------%
c           | Scale the returning eigenvectors so that their |
c           | Euclidean norms are all one. LAPACK subroutine |
c           | dtrevc  returns each eigenvector normalized so  |
c           | that the element of largest magnitude has      |
c           | magnitude 1;                                   |
c           %------------------------------------------------%
c
            iconj = 0
            do 40 j=1, nconv
c
               if ( abs( workl(iheigi+j-1) ) .le. zero ) then
c
c                 %----------------------%
c                 | real eigenvalue case |
c                 %----------------------%
c
                  temp = dnrm2 ( ncv, workl(invsub+(j-1)*ldq),
     &                          ione)
                  call dscal ( ncv, one / temp,
     &                 workl(invsub+(j-1)*ldq), ione )
c
               else
c
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 | columns, we further normalize by the      |
c                 | square root of two.                       |
c                 %-------------------------------------------%
c
                  if (iconj .eq. 0) then
                     temp = dlapy2 (dnrm2 (ncv,
     &                                   workl(invsub+(j-1)*ldq),
     &                                   1),
     &                             dnrm2 (ncv,
     &                                   workl(invsub+j*ldq),
     &                                   1))
                     call dscal (ncv, one/temp,
     &                           workl(invsub+(j-1)*ldq),
     &                           ione)
                     call dscal (ncv, one/temp,
     &                           workl(invsub+j*ldq),
     &                           ione )
                     iconj = 1
                  else
                     iconj = 0
                  end if
c
               end if
c
 40         continue
c
            call dgemv ('T', ncv, nconv, one, workl(invsub),
     &               ldq, workl(ihbds), ione,
     &               zero, workev(ione), 1)
c
            iconj = 0
            do 45 j=1, nconv
               if (workl(iheigi+j-1) .le. zero) then
c
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 %-------------------------------------------%
c
                  if (iconj .eq. 0) then
                     workev(j) = dlapy2 (workev(j), workev(j+1))
                     workev(j+1) = workev(j)
                     iconj = 1
                  else
                     iconj = 0
                  end if
               end if
 45         continue
c
c           %---------------------------------------%
c           | Copy Ritz estimates into workl(ihbds) |
c           %---------------------------------------%
c
            call dcopy (nconv, workev, 1, workl(ihbds), 1)
c
c           %---------------------------------------------------------%
c           | Compute the QR factorization of the eigenvector matrix  |
c           | associated with leading portion of T in the first NCONV |
c           | columns of workl(invsub,ldq).                           |
c           %---------------------------------------------------------%
c
            call dgeqr2 (ncv, nconv , workl(invsub),
     &                   ldq, workev, workev(ncv+1),
     &                   ierr)
c
c           %----------------------------------------------%
c           | * Postmultiply Z by Q.                       |
c           | * Postmultiply Z by R.                       |
c           | The N by NCONV matrix Z is now contains the  |
c           | Ritz vectors associated with the Ritz values |
c           | in workl(iheigr) and workl(iheigi).          |
c           %----------------------------------------------%
c
            call dorm2r ('Right', 'Notranspose', n            ,
     &                   ncv  , nconv        , workl(invsub),
     &                   ldq  , workev       , z            ,
     &                   ldz  , workd(n+1)   , ierr)
c
            call dtrmm ('Right'   , 'Upper'       , 'No transpose',
     &                  'Non-unit', n            , nconv         ,
     &                  one       , workl(invsub), ldq           ,
     &                  z         , ldz)
c
         end if
c
      else
c
c        %------------------------------------------------------%
c        | An approximate invariant subspace is not needed.     |
c        | Place the Ritz values computed DNAUPD  into DR and DI |
c        %------------------------------------------------------%
c
         call dcopy (nconv, workl(ritzr), 1, dr(ione), 1)
         call dcopy (nconv, workl(ritzi), 1, di(ione), 1)
         call dcopy (nconv, workl(ritzr), 1, workl(iheigr), 1)
         call dcopy (nconv, workl(ritzi), 1, workl(iheigi), 1)
         call dcopy (nconv, workl(bounds), 1, workl(ihbds), 1)
      end if
c
c     %------------------------------------------------%
c     | Transform the Ritz values and possibly vectors |
c     | and corresponding error bounds of OP to those  |
c     | of A*x = lambda*B*x.                           |
c     %------------------------------------------------%
c
      if (type .eq. 'REGULR') then
c
         if (rvec)
     &      call dscal (ncv, rnorm, workl(ihbds), ione)
c
      else
c
c        %---------------------------------------%
c        |   A spectral transformation was used. |
c        | * Determine the Ritz estimates of the |
c        |   Ritz values in the original system. |
c        %---------------------------------------%
c
         if (type .eq. 'SHIFTI') then
c
            if (rvec)
     &         call dscal (ncv, rnorm, workl(ihbds), ione)
c
            do 50 k=1, ncv
               temp = dlapy2 ( workl(iheigr+k-1),
     &                        workl(iheigi+k-1) )
               workl(ihbds+k-1) = abs( workl(ihbds+k-1) )
     &                          / temp / temp
 50         continue
c
         else if (type .eq. 'REALPT') then
c
            do 60 k=1, ncv
 60         continue
c
         else if (type .eq. 'IMAGPT') then
c
            do 70 k=1, ncv
 70         continue
c
         end if
c
c        %-----------------------------------------------------------%
c        | *  Transform the Ritz values back to the original system. |
c        |    For TYPE = 'SHIFTI' the transformation is              |
c        |             lambda = 1/theta + sigma                      |
c        |    For TYPE = 'REALPT' or 'IMAGPT' the user must from     |
c        |    Rayleigh quotients or a projection. See remark 3 above.|
c        | NOTES:                                                    |
c        | *The Ritz vectors are not affected by the transformation. |
c        %-----------------------------------------------------------%
c
         if (type .eq. 'SHIFTI') then
c
            do 80 k=1, ncv
               temp = dlapy2 ( workl(iheigr+k-1),
     &                        workl(iheigi+k-1) )
               workl(iheigr+k-1) = workl(iheigr+k-1)/temp/temp
     &                           + sigmar
               workl(iheigi+k-1) = -workl(iheigi+k-1)/temp/temp
     &                           + sigmai
 80         continue
c
            call dcopy (nconv, workl(iheigr), 1, dr(ione), 1)
            call dcopy (nconv, workl(iheigi), 1, di(ione), 1)
c
         else if (type .eq. 'REALPT' .or. type .eq. 'IMAGPT') then
c
            call dcopy (nconv, workl(iheigr), 1, dr(ione), 1)
            call dcopy (nconv, workl(iheigi), 1, di(ione), 1)
c
         end if
c
      end if
c
c     %-------------------------------------------------%
c     | Eigenvector Purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 2.                                   |
c     %-------------------------------------------------%
c
      if (rvec .and. howmny .eq. 'A' .and. type .eq. 'SHIFTI') then
c
c        %------------------------------------------------%
c        | Purify the computed Ritz vectors by adding a   |
c        | little bit of the residual vector:             |
c        |                      T                         |
c        |          resid(:)*( e    s ) / theta           |
c        |                      NCV                       |
c        | where H s = s theta. Remember that when theta  |
c        | has nonzero imaginary part, the corresponding  |
c        | Ritz vector is stored across two columns of Z. |
c        %------------------------------------------------%
c
         iconj = 0
         do 110 j=1, nconv
c            if (workl(iheigi+j-1) .eq. zero) then
            if (abs(workl(iheigi+j-1)) .le. zero) then
               workev(j) =  workl(invsub+(j-1)*ldq+ncv-1) /
     &                      workl(iheigr+j-1)
            else if (iconj .eq. 0) then
               temp = dlapy2 ( workl(iheigr+j-1), workl(iheigi+j-1) )
               workev(j) = ( workl(invsub+(j-1)*ldq+ncv-1) *
     &                       workl(iheigr+j-1) +
     &                       workl(invsub+j*ldq+ncv-1) *
     &                       workl(iheigi+j-1) ) / temp / temp
               workev(j+1) = ( workl(invsub+j*ldq+ncv-1) *
     &                         workl(iheigr+j-1) -
     &                         workl(invsub+(j-1)*ldq+ncv-1) *
     &                         workl(iheigi+j-1) ) / temp / temp
               iconj = 1
            else
               iconj = 0
            end if
 110     continue
c
c        %---------------------------------------%
c        | Perform a rank one update to Z and    |
c        | purify all the Ritz vectors together. |
c        %---------------------------------------%
c
         call dger (n, nconv, one, resid, 1, workev, 1, z, ldz)
c
      end if
c
 9000 continue
c
      return
c
c     %---------------%
c     | End of DNEUPD  |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: dnaupd
c
c\Description:
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  iteration. This subroutine computes approximations to a few eigenpairs
c  of a linear operator "OP" with respect to a semi-inner product defined by
c  a symmetric positive semi-definite real matrix B. B may be the identity
c  matrix. NOTE: If the linear operator "OP" is real and symmetric
c  with respect to the real positive semi-definite symmetric matrix B,
c  i.e. B*OP = (OP`)*B, then subroutine dsaupd should be used instead.
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  dnaupd is usually called iteratively to solve one of the
c  following problems:
c
c  Mode 1:  A*x = lambda*x.
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
c           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M.
c           ===> shift-and-invert mode (in real arithmetic)
c           If OP*x = amu*x, then
c           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
c           Note: If sigma is real, i.e. imaginary part of sigma is zero;
c                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M
c                 amu == 1/(lambda-sigma).
c
c  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
c           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M.
c           ===> shift-and-invert mode (in real arithmetic)
c           If OP*x = amu*x, then
c           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
c
c  Both mode 3 and 4 give the same enhancement to eigenvalues close to
c  the (complex) shift sigma.  However, as lambda goes to infinity,
c  the operator OP in mode 4 dampens the eigenvalues more strongly than
c  does OP defined in mode 3.
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call dnaupd
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first
c          call to dnaupd.  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          dnaupd with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3 and 4, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute the IPARAM(8) real and imaginary parts
c                    of the shifts where INPTR(14) is the pointer
c                    into WORKL for placing the shifts. See Remark
c                    5 below.
c          IDO = 99: done
c          -------------------------------------------------------------
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          'LM' -> want the NEV eigenvalues of largest magnitude.
c          'SM' -> want the NEV eigenvalues of smallest magnitude.
c          'LR' -> want the NEV eigenvalues of largest real part.
c          'SR' -> want the NEV eigenvalues of smallest real part.
c          'LI' -> want the NEV eigenvalues of largest imaginary part.
c          'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c  NEV     Integer.  (INPUT/OUTPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
c
c  TOL     Double precision scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
c          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
c          DEFAULT = DLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH).
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT:
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V. NCV must satisfy the two
c          inequalities 2 <= NCV-NEV and NCV <= N.
c          This will indicate how many Arnoldi vectors are generated
c          at each iteration.  After the startup phase in which NEV
c          Arnoldi vectors are generated, the algorithm generates
c          approximately NCV-NEV Arnoldi vectors at each subsequent update
c          iteration. Most of the cost in generating each Arnoldi vector is
c          in the matrix-vector operation OP*x.
c          NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of Ritz
c          values are kept together. (See remark 4 below)
c
c  V       Double precision array N by NCV.  (OUTPUT)
c          Contains the final set of Arnoldi basis vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are provided by the user via
c                      reverse communication.  The real and imaginary
c                      parts of the NCV eigenvalues of the Hessenberg
c                      matrix H are returned in the part of the WORKL
c                      array corresponding to RITZR and RITZI. See remark
c                      5 below.
c          ISHIFT = 1: exact shifts with respect to the current
c                      Hessenberg matrix H.  This is equivalent to
c                      restarting the iteration with a starting vector
c                      that is a linear combination of approximate Schur
c                      vectors associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced.
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed.
c          On OUTPUT: actual number of Arnoldi update iterations taken.
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used.
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4; See under \Description of dnaupd for the
c          four modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), dnaupd returns NP, the number
c          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
c          5 below.
c
crc          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
crc          OUTPUT: NUMOP  = total number of OP*x operations,
crc                  NUMOPB = total number of B*x operations if BMAT='G',
crc                  NUMREO = total number of steps of re-orthogonalization.
c
c  IPNTR   Integer array of length 14.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by NCV upper Hessenberg matrix
c                    H in WORKL.
c          IPNTR(6): pointer to the real part of the ritz value array
c                    RITZR in WORKL.
c          IPNTR(7): pointer to the imaginary part of the ritz value array
c                    RITZI in WORKL.
c          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
c                    with the Ritz values located in RITZR and RITZI in WORKL.
c
c          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
c
c          Note: IPNTR(9:13) is only referenced by dneupd. See Remark 2 below.
c
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     dneupd if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration. Upon termination
c          WORKD(1:N) contains B*RESID(1:N). If an invariant subspace
c          associated with the converged Ritz values is desired, see remark
c          2 below, subroutine dneupd uses this output.
c          See Data Distribution Note below.
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least 3*NCV**2 + 6*NCV.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the
c                Implicitly restarted Arnoldi iteration. One possibility
c                is to increase the size of NCV relative to NEV.
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iteration
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation;
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization.
c
c\Remarks
c  1. The computed Ritz values are approximate eigenvalues of OP. The
c     selection of WHICH should be made with this in mind when
c     Mode = 3 and 4.  After convergence, approximate eigenvalues of the
c     original problem may be obtained with the ARPACK subroutine dneupd.
c
c  2. If a basis for the invariant subspace corresponding to the converged Ritz
c     values is needed, the user must call dneupd immediately following
c     completion of dnaupd. This is new starting with release 2 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL`
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular
c     linear systems should be solved with L and L` rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L`z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requrement is that NCV > NEV + 2.
c     However, it is recommended that NCV .ge. 2*NEV+1.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.  The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically.
c     See Chapter 8 of Reference 2 for further information.
c
c  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
c     NP = IPARAM(8) real and imaginary parts of the shifts in locations
c         real part                  imaginary part
c         -----------------------    --------------
c     1   WORKL(IPNTR(14))           WORKL(IPNTR(14)+NP)
c     2   WORKL(IPNTR(14)+1)         WORKL(IPNTR(14)+NP+1)
c                        .                          .
c                        .                          .
c                        .                          .
c     NP  WORKL(IPNTR(14)+NP-1)      WORKL(IPNTR(14)+2*NP-1).
c
c     Only complex conjugate pairs of shifts may be applied and the pairs
c     must be placed in consecutive locations. The real part of the
c     eigenvalues of the current upper Hessenberg matrix are located in
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the imaginary part
c     in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). They are ordered
c     according to the order defined by WHICH. The complex conjugate
c     pairs are kept together and the associated Ritz estimates are located in
c     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note:
c
c  Fortran-D syntax:
c  ================
c  Double precision resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c  decompose  d1(n), d2(n,ncv)
c  align      resid(i) with d1(i)
c  align      v(i,j)   with d2(i,j)
c  align      workd(i) with d1(i)     range (1:n)
c  align      workd(i) with d1(i-n)   range (n+1:2*n)
c  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
c  distribute d1(block), d2(block,:)
c  replicated workl(lworkl)
c
c  Cray MPP syntax:
c  ===============
c  Double precision  resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
c  shared     resid(block), v(block,:), workd(block,:)
c  replicated workl(lworkl)
c
c  CM2/CM5 syntax:
c  ==============
c
c-----------------------------------------------------------------------
c
c     include   'ex-nonsym.doc'
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
c     Real Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     dnaup2  ARPACK routine that implements the Implicitly Restarted
c             Arnoldi Iteration.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     12/16/93: Version '1.1'
c
c\SCCS Information: @(#)
c FILE: naupd.F   SID: 2.10   DATE OF SID: 08/23/02   RELEASE: 2
c
c\Remarks
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
c
      implicit none
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(8), ipntr(14)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           zero
      parameter (zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    bounds, ierr, ih, iq, ishift, iupd, iw,
     &           ldh, ldq, levec, mode, mxiter, nb,
     &           nev0, next, np, ritzi, ritzr, j
      save       bounds, ih, iq, ishift, iupd, iw, ldh, ldq,
     &           levec, mode,  mxiter, nb, nev0, next,
     &           np, ritzi, ritzr
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dnaup2
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlamch
      external   dlamch
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (levec .gt. 10000000) then
        goto 9000
      end if
c
      if (ido .eq. 0) then
c
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
c        %----------------%
c        | Error checking |
c        %----------------%
c
         ierr   = 0
         ishift = iparam(1)
c         levec  = iparam(2)
         mxiter = iparam(3)
c         nb     = iparam(4)
         nb     = 1
c
c        %--------------------------------------------%
c        | Revision 2 performs only implicit restart. |
c        %--------------------------------------------%
c
         iupd   = 1
         mode   = iparam(7)
c
         if (n .le. 0) then
            ierr = -1
         else if (nev .le. 0) then
            ierr = -2
         else if (ncv .le. nev+1 .or.  ncv .gt. n) then
            ierr = -3
         else if (mxiter .le.          0) then
            ierr = 4
         else if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LR' .and.
     &       which .ne. 'SR' .and.
     &       which .ne. 'LI' .and.
     &       which .ne. 'SI') then
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
c
c        %------------%
c        | Error Exit |
c        %------------%
c
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
c
c        %------------------------%
c        | Set default parameters |
c        %------------------------%
c
         if (nb .le. 0)           nb = 1
         if (tol .le. zero)       tol = dlamch('EpsMach')
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        | NEV0 is the local variable designating the   |
c        | size of the invariant subspace desired.      |
c        %----------------------------------------------%
c
         np     = ncv - nev
         nev0   = nev
c
c        %-----------------------------%
c        | Zero out internal workspace |
c        %-----------------------------%
c
         do 10 j = 1, 3*ncv**2 + 6*ncv
            workl(j) = zero
  10     continue
c
c        %-------------------------------------------------------------%
c        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        |
c        | etc... and the remaining workspace.                         |
c        | Also update pointer to be used on output.                   |
c        | Memory is laid out as follows:                              |
c        | workl(1:ncv*ncv) := generated Hessenberg matrix             |
c        | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary        |
c        |                                   parts of ritz values      |
c        | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds        |
c        | workl(ncv*ncv+3*ncv+1:2*ncv*ncv+3*ncv) := rotation matrix Q |
c        | workl(2*ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) := workspace       |
c        | The final workspace is needed by subroutine dneigh called   |
c        | by dnaup2. Subroutine dneigh calls LAPACK routines for      |
c        | calculating eigenvalues and the last row of the eigenvector |
c        | matrix.                                                     |
c        %-------------------------------------------------------------%
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritzr  = ih     + ldh*ncv
         ritzi  = ritzr  + ncv
         bounds = ritzi  + ncv
         iq     = bounds + ncv
         iw     = iq     + ldq*ncv
         next   = iw     + ncv**2 + 3*ncv
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritzr
         ipntr(7) = ritzi
         ipntr(8) = bounds
         ipntr(14) = iw
c
      end if
c
c     %-------------------------------------------------------%
c     | Carry out the Implicitly restarted Arnoldi Iteration. |
c     %-------------------------------------------------------%
c
      call dnaup2
     &   ( ido, bmat, n, which, nev0, np, tol, resid, mode,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritzr),
     &     workl(ritzi), workl(bounds), workl(iq), ldq, workl(iw),
     &     ipntr, workd, info )
c
c     %--------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication |
c     | to compute operations involving OP or shifts.    |
c     %--------------------------------------------------%
c
      if (ido .eq. 3) iparam(8) = np
      if (ido .ne. 99) go to 9000
c
      iparam(3) = mxiter
      iparam(5) = np
c
c     %------------------------------------%
c     | Exit if there was an informational |
c     | error within dnaup2.               |
c     %------------------------------------%
c
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
c
c        %--------------------------------------------------------%
c        | Version Number & Version Date are defined in version.h |
c        %--------------------------------------------------------%
c
 9000 continue
c
      return
c
c     %---------------%
c     | End of dnaupd |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: dnaup2
c
c\Description:
c  Intermediate level interface called by dnaupd.
c
c\Usage:
c  call dnaup2
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS,
c       Q, LDQ, WORKL, IPNTR, WORKD, INFO )
c
c\Arguments
c
c  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dnaupd.
c  MODE, ISHIFT, MXITER: see the definition of IPARAM in dnaupd.
c
c  NP      Integer.  (INPUT/OUTPUT)
c          Contains the number of implicit shifts to apply during
c          each Arnoldi iteration.
c          If ISHIFT=1, NP is adjusted dynamically at each iteration
c          to accelerate convergence and prevent stagnation.
c          This is also roughly equal to the number of matrix-vector
c          products (involving the operator OP) per Arnoldi iteration.
c          The logic for adjusting is contained within the current
c          subroutine.
c          If ISHIFT=0, NP is the number of shifts the user needs
c          to provide via reverse comunication. 0 < NP < NCV-NEV.
c          NP may be less than NCV-NEV for two reasons. The first, is
c          to keep complex conjugate pairs of "wanted" Ritz values
c          together. The second, is that a leading block of the current
c          upper Hessenberg matrix has split off and contains "unwanted"
c          Ritz values.
c          Upon termination of the IRA iteration, NP contains the number
c          of "converged" wanted Ritz values.
c
c  IUPD    Integer.  (INPUT)
c          IUPD .EQ. 0: use explicit restart instead implicit update.
c          IUPD .NE. 0: use implicit update.
c
c  V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT)
c          The Arnoldi basis vectors are returned in the first NEV
c          columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Double precision (NEV+NP) by (NEV+NP) array.  (OUTPUT)
c          H is used to store the generated upper Hessenberg matrix
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RITZR,  Double precision arrays of length NEV+NP.  (OUTPUT)
c  RITZI   RITZR(1:NEV) (resp. RITZI(1:NEV)) contains the real (resp.
c          imaginary) part of the computed Ritz values of OP.
c
c  BOUNDS  Double precision array of length NEV+NP.  (OUTPUT)
c          BOUNDS(1:NEV) contain the error bounds corresponding to
c          the computed Ritz values.
c
c  Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
c          Private (replicated) work array used to accumulate the
c          rotation in the shift application step.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length at least
c          (NEV+NP)**2 + 3*(NEV+NP).  (INPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  It is used in shifts calculation, shifts
c          application and convergence checking.
c
c          On exit, the last 3*(NEV+NP) locations of WORKL contain
c          the Ritz values (real,imaginary) and associated Ritz
c          estimates of the current Hessenberg matrix.  They are
c          listed in the same order as returned from dneigh.
c
c          If ISHIFT .EQ. O and IDO .EQ. 3, the first 2*NP locations
c          of WORKL are used in reverse communication to hold the user
c          supplied shifts.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD for
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c
c  WORKD   Double precision work array of length 3*N.  (WORKSPACE)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note in DNAUPD.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =     0: Normal return.
c          =     1: Maximum number of iterations taken.
c                   All possible eigenvalues of OP has been found.
c                   NP returns the number of converged Ritz values.
c          =     2: No shifts could be applied.
c          =    -8: Error return from LAPACK eigenvalue calculation;
c                   This should never happen.
c          =    -9: Starting vector is zero.
c          = -9999: Could not build an Arnoldi factorization.
c                   Size that was built in returned in NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     dgetv0  ARPACK initial vector generation routine.
c     dnaitr  ARPACK Arnoldi factorization routine.
c     dnapps  ARPACK application of implicit shifts routine.
c     dnconv  ARPACK convergence of Ritz values routine.
c     dneigh  ARPACK compute Ritz values and error bounds routine.
c     dngets  ARPACK reorder Ritz values and error bounds routine.
c     dsortc  ARPACK sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dmout   ARPACK utility routine that prints matrices
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     ddot    Level 1 BLAS that computes the scalar product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dswap   Level 1 BLAS that swaps two vectors.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: naup2.F   SID: 2.8   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnaup2
     &   ( ido, bmat, n, which, nev, np, tol, resid, mode,
     &     ishift, mxiter, v, ldv, h, ldh, ritzr, ritzi, bounds,
     &     q, ldq, workl, ipntr, workd, info )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      implicit none
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ishift, mode, ldh, ldq, ldv, mxiter,
     &           n, nev, np
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(13)
      Double precision
     &           bounds(nev+np), h(ldh,nev+np), q(ldq,nev+np), resid(n),
     &           ritzi(nev+np), ritzr(nev+np), v(ldv,nev+np),
     &           workd(3*n), workl( (nev+np)*(nev+np+3) )
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           zero
      parameter (zero = 0.0D+0)
      integer    dnazero
      integer(4) ione
      parameter (dnazero = 0, ione = 1)

c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  wprime*2
      logical    cnorm , getv0, initv, update, ushift, dnatrue
      integer    ierr  , iter , j    , kplusp, nconv,
     &           nevbef, nev0 , np0  , nptemp, numcnv
      Double precision
     &           rnorm , temp , eps23
      save       cnorm , getv0, initv, update, ushift,
     &           rnorm , iter , eps23, kplusp,  nconv ,
     &           nevbef, nev0 , np0  , numcnv
c
c     %-----------------------%
c     | Local array arguments |
c     %-----------------------%
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy , dgetv0, dnaitr, dnconv, dneigh,
     &           dngets, dnapps
c     &           dngets, dnapps, second
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
c
      Double precision
     &           ddot, dnrm2, dlapy2, dlamch
      external   ddot, dnrm2, dlapy2, dlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    min, max, abs, sqrt
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      dnatrue = .true.
c
      if (ido .eq. 0) then
c
c        %-------------------------------------%
c        | Get the machine dependent constant. |
c        %-------------------------------------%
c
         eps23 = dlamch('Epsilon-Machine')
         eps23 = eps23**(2.0D+0 / 3.0D+0)
c
         nev0   = nev
         np0    = np
c
c        %-------------------------------------%
c        | kplusp is the bound on the largest  |
c        |        Lanczos factorization built. |
c        | nconv is the current number of      |
c        |        "converged" eigenvlues.      |
c        | iter is the counter on the current  |
c        |      iteration step.                |
c        %-------------------------------------%
c
         kplusp = nev + np
         nconv  = 0
         iter   = 0
c
c        %---------------------------------------%
c        | Set flags for computing the first NEV |
c        | steps of the Arnoldi factorization.   |
c        %---------------------------------------%
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
c
         if (info .ne. 0) then
c
c           %--------------------------------------------%
c           | User provides the initial residual vector. |
c           %--------------------------------------------%
c
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
c
c     %---------------------------------------------%
c     | Get a possibly random starting vector and   |
c     | force it into the range of the operator OP. |
c     %---------------------------------------------%
c
c   10 continue
c
      if (getv0) then
         call dgetv0 (ido, bmat, initv, n, ione, v, ldv, resid, rnorm,
     &                ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (abs(rnorm) .le. zero) then
c         if (rnorm .eq. zero) then
c
c           %-----------------------------------------%
c           | The initial vector is zero. Error exit. |
c           %-----------------------------------------%
c
            info = -9
            go to 1100
         end if
         getv0 = .false.
         ido  = 0
      end if
c
c     %-----------------------------------%
c     | Back from reverse communication : |
c     | continue with update step         |
c     %-----------------------------------%
c
      if (update) go to 20
c
c     %-------------------------------------------%
c     | Back from computing user specified shifts |
c     %-------------------------------------------%
c
      if (ushift) go to 50
c
c     %-------------------------------------%
c     | Back from computing residual norm   |
c     | at the end of the current iteration |
c     %-------------------------------------%
c
      if (cnorm)  go to 100
c
c     %----------------------------------------------------------%
c     | Compute the first NEV steps of the Arnoldi factorization |
c     %----------------------------------------------------------%
c
      call dnaitr (ido, bmat, n, dnazero, nev, mode, resid, rnorm,
     &             v, ldv, h, ldh, ipntr, workd, info)
c
c     %---------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication  |
c     | to compute operations involving OP and possibly B |
c     %---------------------------------------------------%
c
      if (ido .ne. 99) go to 9000
c
      if (info .gt. 0) then
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
c
c     %--------------------------------------------------------------%
c     |                                                              |
c     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
c     |           Each iteration implicitly restarts the Arnoldi     |
c     |           factorization in place.                            |
c     |                                                              |
c     %--------------------------------------------------------------%
c
 1000 continue
c
         iter = iter + 1
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        | Adjust NP since NEV might have been updated by last call  |
c        | to the shift application routine dnapps.                  |
c        %-----------------------------------------------------------%
c
         np  = kplusp - nev
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        %-----------------------------------------------------------%
c
         ido = 0
   20    continue
         update = .true.
c
         call dnaitr (ido  , bmat, n  , nev, np , mode , resid,
     &                rnorm, v   , ldv, h  , ldh, ipntr, workd,
     &                info)
c
c        %---------------------------------------------------%
c        | ido .ne. 99 implies use of reverse communication  |
c        | to compute operations involving OP and possibly B |
c        %---------------------------------------------------%
c
         if (ido .ne. 99) go to 9000
c
         if (info .gt. 0) then
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
c
c        %--------------------------------------------------------%
c        | Compute the eigenvalues and corresponding error bounds |
c        | of the current upper Hessenberg matrix.                |
c        %--------------------------------------------------------%
c
         call dneigh (rnorm, kplusp, h, ldh, ritzr, ritzi, bounds,
     &                q, ldq, workl, ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %----------------------------------------------------%
c        | Make a copy of eigenvalues and corresponding error |
c        | bounds obtained from dneigh.                       |
c        %----------------------------------------------------%
c
         call dcopy(kplusp, ritzr, 1, workl(kplusp**2+1), 1)
         call dcopy(kplusp, ritzi, 1, workl(kplusp**2+kplusp+1), 1)
         call dcopy(kplusp, bounds, 1, workl(kplusp**2+2*kplusp+1), 1)
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The wanted part of the spectrum and corresponding |
c        | error bounds are in the last NEV loc. of RITZR,   |
c        | RITZI and BOUNDS respectively. The variables NEV  |
c        | and NP may be updated if the NEV-th wanted Ritz   |
c        | value has a non zero imaginary part. In this case |
c        | NEV is increased by one and NP decreased by one.  |
c        | NOTE: The last two arguments of dngets are no     |
c        | longer used as of version 2.1.                    |
c        %---------------------------------------------------%
c
         nev = nev0
         np = np0
         numcnv = nev
         call dngets (ishift, which, nev, np, ritzr, ritzi,
     &                bounds)
         if (nev .eq. nev0+1) numcnv = nev0+1
c
c        %-------------------%
c        | Convergence test. |
c        %-------------------%
c
         call dcopy (nev, bounds(np+1), 1, workl(2*np+1), 1)
         call dnconv (nev, ritzr(np+1), ritzi(np+1), workl(2*np+1),
     &        tol, nconv)
c
c        %---------------------------------------------------------%
c        | Count the number of unwanted Ritz values that have zero |
c        | Ritz estimates. If any Ritz estimates are equal to zero |
c        | then a leading block of H of order equal to at least    |
c        | the number of Ritz values with zero Ritz estimates has  |
c        | split off. None of these Ritz values may be removed by  |
c        | shifting. Decrease NP the number of shifts to apply. If |
c        | no shifts may be applied, then prepare to exit          |
c        %---------------------------------------------------------%
c
         nptemp = np
         do 30 j=1, nptemp
c            if (bounds(j) .eq. zero) then
            if (abs( bounds(j) ) .le. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
c
         if ( (nconv .ge. numcnv) .or.
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
c
c
c           %------------------------------------------------%
c           | Prepare to exit. Put the converged Ritz values |
c           | and corresponding bounds in RITZ(1:NCONV) and  |
c           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
c           | careful when NCONV > NP                        |
c           %------------------------------------------------%
c
c           %------------------------------------------%
c           |  Use h( 3,1 ) as storage to communicate  |
c           |  rnorm to _neupd if needed               |
c           %------------------------------------------%

            h(3,1) = rnorm
c
c           %----------------------------------------------%
c           | To be consistent with dngets, we first do a  |
c           | pre-processing sort in order to keep complex |
c           | conjugate pairs together.  This is similar   |
c           | to the pre-processing sort used in dngets    |
c           | except that the sort is done in the opposite |
c           | order.                                       |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SR'
            if (which .eq. 'SM') wprime = 'LR'
            if (which .eq. 'LR') wprime = 'SM'
            if (which .eq. 'SR') wprime = 'LM'
            if (which .eq. 'LI') wprime = 'SM'
            if (which .eq. 'SI') wprime = 'LM'
c
            call dsortc (wprime, dnatrue, kplusp, ritzr, ritzi, bounds)
c
c           %----------------------------------------------%
c           | Now sort Ritz values so that converged Ritz  |
c           | values appear within the first NEV locations |
c           | of ritzr, ritzi and bounds, and the most     |
c           | desired one appears at the front.            |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SM'
            if (which .eq. 'SM') wprime = 'LM'
            if (which .eq. 'LR') wprime = 'SR'
            if (which .eq. 'SR') wprime = 'LR'
            if (which .eq. 'LI') wprime = 'SI'
            if (which .eq. 'SI') wprime = 'LI'
c
            call dsortc(wprime, dnatrue, kplusp, ritzr, ritzi, bounds)
c
c           %--------------------------------------------------%
c           | Scale the Ritz estimate of each Ritz value       |
c           | by 1 / max(eps23,magnitude of the Ritz value).   |
c           %--------------------------------------------------%
c
            do 35 j = 1, numcnv
                temp = max(eps23,dlapy2(ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)/temp
 35         continue
c
c           %----------------------------------------------------%
c           | Sort the Ritz values according to the scaled Ritz  |
c           | esitmates.  This will push all the converged ones  |
c           | towards the front of ritzr, ritzi, bounds          |
c           | (in the case when NCONV < NEV.)                    |
c           %----------------------------------------------------%
c
            wprime = 'LR'
            call dsortc(wprime, dnatrue, numcnv, bounds, ritzr, ritzi)
c
c           %----------------------------------------------%
c           | Scale the Ritz estimate back to its original |
c           | value.                                       |
c           %----------------------------------------------%
c
            do 40 j = 1, numcnv
                temp = max(eps23, dlapy2(ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)*temp
 40         continue
c
c           %------------------------------------------------%
c           | Sort the converged Ritz values again so that   |
c           | the "threshold" value appears at the front of  |
c           | ritzr, ritzi and bound.                        |
c           %------------------------------------------------%
c
            call dsortc(which, dnatrue, nconv, ritzr, ritzi, bounds)
c
c
c           %------------------------------------%
c           | Max iterations have been exceeded. |
c           %------------------------------------%
c
            if (iter .gt. mxiter .and. nconv .lt. numcnv) info = 1
c
c           %---------------------%
c           | No shifts to apply. |
c           %---------------------%
c
            if (np .eq. 0 .and. nconv .lt. numcnv) info = 2
c
            np = nconv
            go to 1100
c
         else if ( (nconv .lt. numcnv) .and. (ishift .eq. 1) ) then
c
c           %-------------------------------------------------%
c           | Do not have all the requested eigenvalues yet.  |
c           | To prevent possible stagnation, adjust the size |
c           | of NEV.                                         |
c           %-------------------------------------------------%
c
            nevbef = nev
            nev = nev + min(nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 3) then
               nev = 2
            end if
            np = kplusp - nev
c
c           %---------------------------------------%
c           | If the size of NEV was just increased |
c           | resort the eigenvalues.               |
c           %---------------------------------------%
c
            if (nevbef .lt. nev)
     &         call dngets (ishift, which, nev, np, ritzr, ritzi,
     &              bounds)
c
         end if
c
         if (ishift .eq. 0) then
c
c           %-------------------------------------------------------%
c           | User specified shifts: reverse comminucation to       |
c           | compute the shifts. They are returned in the first    |
c           | 2*NP locations of WORKL.                              |
c           %-------------------------------------------------------%
c
            ushift = .true.
            ido = 3
            go to 9000
         end if
c
   50    continue
c
c        %------------------------------------%
c        | Back from reverse communication;   |
c        | User specified shifts are returned |
c        | in WORKL(1:2*NP)                   |
c        %------------------------------------%
c
         ushift = .false.
c
         if ( ishift .eq. 0 ) then
c
c            %----------------------------------%
c            | Move the NP shifts from WORKL to |
c            | RITZR, RITZI to free up WORKL    |
c            | for non-exact shift case.        |
c            %----------------------------------%
c
             call dcopy (np, workl,       1, ritzr(1), 1)
             call dcopy (np, workl(np+1), 1, ritzi(1), 1)
         end if
c
c
c        %---------------------------------------------------------%
c        | Apply the NP implicit shifts by QR bulge chasing.       |
c        | Each shift is applied to the whole upper Hessenberg     |
c        | matrix H.                                               |
c        | The first 2*N locations of WORKD are used as workspace. |
c        %---------------------------------------------------------%
c
         call dnapps (n, nev, np, ritzr, ritzi, v, ldv,
     &                h, ldh, resid, q, ldq, workl, workd)
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to dnaitr.  |
c        %---------------------------------------------%
c
         cnorm = .true.
         if (bmat .eq. 'G') then
            call dcopy (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
c
c           %----------------------------------%
c           | Exit in order to compute B*RESID |
c           %----------------------------------%
c
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(1), 1)
         end if
c
  100    continue
c
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(1:N) := B*RESID            |
c        %----------------------------------%
c
         if (bmat .eq. 'G') then
            rnorm = ddot (n, resid, 1, workd, 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
         cnorm = .false.
c
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 1100 continue
c
      mxiter = iter
      nev = numcnv
c
 1200 continue
      ido = 99
c
c     %------------%
c     | Error Exit |
c     %------------%
c
 9000 continue
c
c     %---------------%
c     | End of dnaup2 |
c     %---------------%
c
      return
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dnapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP implicit shifts resulting in
c
c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
c
c  where Q is an orthogonal matrix which is the product of rotations
c  and reflections resulting from the NP bulge chage sweeps.
c  The updated Arnoldi factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call dnapps
c     ( N, KEV, NP, SHIFTR, SHIFTI, V, LDV, H, LDH, RESID, Q, LDQ,
c       WORKL, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. size of matrix A.
c
c  KEV     Integer.  (INPUT/OUTPUT)
c          KEV+NP is the size of the input matrix H.
c          KEV is the size of the updated matrix HNEW.  KEV is only
c          updated on ouput when fewer than NP shifts are applied in
c          order to keep the conjugate pair together.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.
c
c  SHIFTR, Double precision array of length NP.  (INPUT)
c  SHIFTI  Real and imaginary part of the shifts to be applied.
c          Upon, entry to dnapps, the shifts must be sorted so that the
c          conjugate pairs are in consecutive locations.
c
c  V       Double precision N by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, V contains the current KEV+NP Arnoldi vectors.
c          On OUTPUT, V contains the updated KEV Arnoldi vectors
c          in the first KEV columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Double precision (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, H contains the current KEV+NP by KEV+NP upper
c          Hessenber matrix of the Arnoldi factorization.
c          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
c          matrix in the KEV leading submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT, RESID contains the the residual vector r_{k+p}.
c          On OUTPUT, RESID is the update residual vector rnew_{k}
c          in the first KEV locations.
c
c  Q       Double precision KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations and reflections
c          during the bulge chase sweep.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length (KEV+NP).  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  WORKD   Double precision work array of length 2*N.  (WORKSPACE)
c          Distributed array used in the application of the accumulated
c          orthogonal matrix Q.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dmout   ARPACK utility routine that prints matrices.
c     dvout   ARPACK utility routine that prints vectors.
c     dlabad  LAPACK routine that computes machine constants.
c     dlacpy  LAPACK matrix copy routine.
c     dlamch  LAPACK routine that determines machine constants.
c     dlanhs  LAPACK routine that computes various norms of a matrix.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlarf   LAPACK routine that applies Householder reflection to
c             a matrix.
c     dlarfg  LAPACK Householder reflection construction routine.
c     dlartg  LAPACK Givens rotation construction routine.
c     dlaset  LAPACK matrix initialization routine.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     dscal   Level 1 BLAS that scales a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c
c\SCCS Information: @(#)
c FILE: napps.F   SID: 2.4   DATE OF SID: 3/28/97   RELEASE: 2
c
c\Remarks
c  1. In this version, each shift is applied to all the sublocks of
c     the Hessenberg matrix H and not just to the submatrix that it
c     comes from. Deflation as in LAPACK routine dlahqr (QR algorithm
c     for upper Hessenberg matrices ) is used.
c     The subdiagonals of H are enforced to be non-negative.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnapps
     &   ( n, kev, np, shiftr, shifti, v, ldv, h, ldh, resid, q, ldq,
     &     workl, workd )

      implicit none
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    kev, ldh, ldq, ldv, n, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           h(ldh,kev+np), resid(n), shifti(np), shiftr(np),
     &           v(ldv,kev+np), q(ldq,kev+np), workd(2*n), workl(kev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
      integer(4)  ione
      parameter (ione = 1)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      integer    i, iend, ir, istart, j, jj, kplusp, nr
      logical    cconj, first
      Double precision
     &           c, f, g, h11, h12, h21, h22, h32, ovfl, r, s, sigmai,
     &           sigmar, smlnum, ulp, unfl, u(3), t, tau, tst1
      save       first, ovfl, smlnum, ulp, unfl
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   daxpy, dcopy, dscal, dlacpy, dlarfg, dlarf,
     &           dlaset, dlabad, dlartg
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
c
      Double precision
     &           dlamch, dlanhs, dlapy2
      external   dlamch, dlanhs, dlapy2
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs, max, min
c
c     %----------------%
c     | Data statments |
c     %----------------%
c
      data       first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
c
c        %-----------------------------------------------%
c        | Set machine-dependent constants for the       |
c        | stopping criterion. If norm(H) <= sqrt(OVFL), |
c        | overflow should not occur.                    |
c        | REFERENCE: LAPACK subroutine dlahqr           |
c        %-----------------------------------------------%
c
         unfl = dlamch( 'safe minimum' )
         ovfl = one / unfl
         call dlabad( unfl, ovfl )
         ulp = dlamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
c      call second (t0)
      kplusp = kev + np
c
c     %--------------------------------------------%
c     | Initialize Q to the identity to accumulate |
c     | the rotations and reflections              |
c     %--------------------------------------------%
c
      call dlaset ('All', kplusp, kplusp, zero, one, q, ldq)
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c
c     %----------------------------------------------%
c     | Chase the bulge with the application of each |
c     | implicit shift. Each shift is applied to the |
c     | whole matrix including each block.           |
c     %----------------------------------------------%
c
      cconj = .false.
      do 110 jj = 1, np
         sigmar = shiftr(jj)
         sigmai = shifti(jj)
c
c        %-------------------------------------------------%
c        | The following set of conditionals is necessary  |
c        | in order that complex conjugate pairs of shifts |
c        | are applied together or not at all.             |
c        %-------------------------------------------------%
c
         if ( cconj ) then
c
c           %-----------------------------------------%
c           | cconj = .true. means the previous shift |
c           | had non-zero imaginary part.            |
c           %-----------------------------------------%
c
            cconj = .false.
            go to 110
         else if ( jj .lt. np .and. abs( sigmai ) .gt. zero ) then
c
c           %------------------------------------%
c           | Start of a complex conjugate pair. |
c           %------------------------------------%
c
            cconj = .true.
         else if ( jj .eq. np .and. abs( sigmai ) .gt. zero ) then
c
c           %----------------------------------------------%
c           | The last shift has a nonzero imaginary part. |
c           | Don't apply it; thus the order of the        |
c           | compressed H is order KEV+1 since only np-1  |
c           | were applied.                                |
c           %----------------------------------------------%
c
            kev = kev + 1
            go to 110
         end if
         istart = 1
   20    continue
c
c        %--------------------------------------------------%
c        | if sigmai = 0 then                               |
c        |    Apply the jj-th shift ...                     |
c        | else                                             |
c        |    Apply the jj-th and (jj+1)-th together ...    |
c        |    (Note that jj < np at this point in the code) |
c        | end                                              |
c        | to the current block of H. The next do loop      |
c        | determines the current block ;                   |
c        %--------------------------------------------------%
c
         do 30 i = istart, kplusp-1
c
c           %----------------------------------------%
c           | Check for splitting and deflation. Use |
c           | a standard test as in the QR algorithm |
c           | REFERENCE: LAPACK subroutine dlahqr    |
c           %----------------------------------------%
c
            tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
c            if( tst1.eq.zero )
            if( abs(tst1) .le. zero )
     &         tst1 = dlanhs( '1', kplusp-jj+1, h, ldh, workl )
            if( abs( h( i+1,i ) ).le.max( ulp*tst1, smlnum ) ) then
               iend = i
               h(i+1,i) = zero
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
c
c        %------------------------------------------------%
c        | No reason to apply a shift to block of order 1 |
c        %------------------------------------------------%
c
         if ( istart .eq. iend ) go to 100
c
c        %------------------------------------------------------%
c        | If istart + 1 = iend then no reason to apply a       |
c        | complex conjugate pair of shifts on a 2 by 2 matrix. |
c        %------------------------------------------------------%
c
         if ( istart + 1 .eq. iend .and. abs( sigmai ) .gt. zero )
     &      go to 100
c
         h11 = h(istart,istart)
         h21 = h(istart+1,istart)
         if ( abs( sigmai ) .le. zero ) then
c
c           %---------------------------------------------%
c           | Real-valued shift ==> apply single shift QR |
c           %---------------------------------------------%
c
            f = h11 - sigmar
            g = h21
c
            do 80 i = istart, iend-1
c
c              %-----------------------------------------------------%
c              | Contruct the plane rotation G to zero out the bulge |
c              %-----------------------------------------------------%
c
               call dlartg (f, g, c, s, r)
               if (i .gt. istart) then
c
c                 %-------------------------------------------%
c                 | The following ensures that h(1:iend-1,1), |
c                 | the first iend-2 off diagonal of elements |
c                 | H, remain non negative.                   |
c                 %-------------------------------------------%
c
                  if (r .lt. zero) then
                     r = -r
                     c = -c
                     s = -s
                  end if
                  h(i,i-1) = r
                  h(i+1,i-1) = zero
               end if
c
c              %---------------------------------------------%
c              | Apply rotation to the left of H;  H <- G'*H |
c              %---------------------------------------------%
c
               do 50 j = i, kplusp
                  t        =  c*h(i,j) + s*h(i+1,j)
                  h(i+1,j) = -s*h(i,j) + c*h(i+1,j)
                  h(i,j)   = t
   50          continue
c
c              %---------------------------------------------%
c              | Apply rotation to the right of H;  H <- H*G |
c              %---------------------------------------------%
c
               do 60 j = 1, min(i+2,iend)
                  t        =  c*h(j,i) + s*h(j,i+1)
                  h(j,i+1) = -s*h(j,i) + c*h(j,i+1)
                  h(j,i)   = t
   60          continue
c
c              %----------------------------------------------------%
c              | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c              %----------------------------------------------------%
c
               do 70 j = 1, min( i+jj, kplusp )
                  t        =   c*q(j,i) + s*q(j,i+1)
                  q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                  q(j,i)   = t
   70          continue
c
c              %---------------------------%
c              | Prepare for next rotation |
c              %---------------------------%
c
               if (i .lt. iend-1) then
                  f = h(i+1,i)
                  g = h(i+2,i)
               end if
   80       continue
c
c           %-----------------------------------%
c           | Finished applying the real shift. |
c           %-----------------------------------%
c
         else
c
c           %----------------------------------------------------%
c           | Complex conjugate shifts ==> apply double shift QR |
c           %----------------------------------------------------%
c
            h12 = h(istart,istart+1)
            h22 = h(istart+1,istart+1)
            h32 = h(istart+2,istart+1)
c
c           %---------------------------------------------------------%
c           | Compute 1st column of (H - shift*I)*(H - conj(shift)*I) |
c           %---------------------------------------------------------%
c
            s    = 2.0*sigmar
            t = dlapy2 ( sigmar, sigmai )
            u(1) = ( h11 * (h11 - s) + t * t ) / h21 + h12
            u(2) = h11 + h22 - s
            u(3) = h32
c
            do 90 i = istart, iend-1
c
               nr = min ( 3, iend-i+1 )
c
c              %-----------------------------------------------------%
c              | Construct Householder reflector G to zero out u(1). |
c              | G is of the form I - tau*( 1 u )' * ( 1 u' ).       |
c              %-----------------------------------------------------%
c
               call dlarfg ( nr, u(1), u(2), 1, tau )
c
               if (i .gt. istart) then
                  h(i,i-1)   = u(1)
                  h(i+1,i-1) = zero
                  if (i .lt. iend-1) h(i+2,i-1) = zero
               end if
               u(1) = one
c
c              %--------------------------------------%
c              | Apply the reflector to the left of H |
c              %--------------------------------------%
c
               call dlarf ('Left', nr, kplusp-i+1, u, 1, tau,
     &                     h(i,i), ldh, workl)
c
c              %---------------------------------------%
c              | Apply the reflector to the right of H |
c              %---------------------------------------%
c
               ir = min ( i+3, iend )
               call dlarf ('Right', ir, nr, u, 1, tau,
     &                     h(1,i), ldh, workl)
c
c              %-----------------------------------------------------%
c              | Accumulate the reflector in the matrix Q;  Q <- Q*G |
c              %-----------------------------------------------------%
c
               call dlarf ('Right', kplusp, nr, u, 1, tau,
     &                     q(1,i), ldq, workl)
c
c              %----------------------------%
c              | Prepare for next reflector |
c              %----------------------------%
c
               if (i .lt. iend-1) then
                  u(1) = h(i+1,i)
                  u(2) = h(i+2,i)
                  if (i .lt. iend-2) u(3) = h(i+3,i)
               end if
c
   90       continue
c
c           %--------------------------------------------%
c           | Finished applying a complex pair of shifts |
c           | to the current block                       |
c           %--------------------------------------------%
c
         end if
c
  100    continue
c
c        %---------------------------------------------------------%
c        | Apply the same shift to the next block if there is any. |
c        %---------------------------------------------------------%
c
         istart = iend + 1
         if (iend .lt. kplusp) go to 20
c
c        %---------------------------------------------%
c        | Loop back to the top to get the next shift. |
c        %---------------------------------------------%
c
  110 continue
c
c     %--------------------------------------------------%
c     | Perform a similarity transformation that makes   |
c     | sure that H will have non negative sub diagonals |
c     %--------------------------------------------------%
c
      do 120 j=1,kev
         if ( h(j+1,j) .lt. zero ) then
            call dscal( kplusp-j+1, -one, h(j+1,j), ldh)
            call dscal( min(j+2, kplusp), -one, h(1,j+1),
     &                  ione)
            call dscal( min(j+np+1,kplusp), -one, q(1,j+1),
     &                  ione)
cv              call dscal( min(j+2, kplusp), -one, h(1,j+1), 1 )
cv              call dscal( min(j+np+1,kplusp), -one, q(1,j+1), 1 )
         end if
 120  continue
c
      do 130 i = 1, kev
c
c        %--------------------------------------------%
c        | Final check for splitting and deflation.   |
c        | Use a standard test as in the QR algorithm |
c        | REFERENCE: LAPACK subroutine dlahqr        |
c        %--------------------------------------------%
c
         tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
c         if( tst1.eq.zero )
         if( abs(tst1) .le. zero)
     &       tst1 = dlanhs( '1', kev, h, ldh, workl )
         if( h( i+1,i ) .le. max( ulp*tst1, smlnum ) )
     &       h(i+1,i) = zero
 130  continue
c
c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is needed in the residual update since we  |
c     | cannot GUARANTEE that the corresponding entry   |
c     | of H would be zero as in exact arithmetic.      |
c     %-------------------------------------------------%
c
      if (h(kev+1,kev) .gt. zero)
     &    call dgemv ('N', n, kplusp, one, v, ldv, q(1,kev+1), 1, zero,
     &                workd(n+1), 1)
c
c     %----------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order       |
c     | taking advantage of the upper Hessenberg structure of Q. |
c     %----------------------------------------------------------%
c
      do 140 i = 1, kev
         call dgemv ('N', n, kplusp-i+1, one, v, ldv,
     &        q(1,kev-i+1), 1, zero, workd(1), 1)
         call dcopy (n, workd, 1, v(1,kplusp-i+1), 1)
  140 continue
c
c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%
c
      call dlacpy ('A', n, kev, v(1,kplusp-kev+1), ldv, v, ldv)
c
c     %--------------------------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
c     %--------------------------------------------------------------%
c
      if (h(kev+1,kev) .gt. zero)
     &   call dcopy (n, workd(n+1), 1, v(1,kev+1), 1)
c
c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kplusp}'*Q)*e_{kev} |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%
c
      call dscal (n, q(kplusp,kev), resid(1), ione)
      if (h(kev+1,kev) .gt. zero)
     &   call daxpy (n, h(kev+1,kev), v(1,kev+1), 1, resid, 1)
c
 9000 continue
c
      return
c
c     %---------------%
c     | End of dnapps |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dnconv
c
c\Description:
c  Convergence testing for the nonsymmetric Arnoldi eigenvalue routine.
c
c\Usage:
c  call dnconv
c     ( N, RITZR, RITZI, BOUNDS, TOL, NCONV )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Number of Ritz values to check for convergence.
c
c  RITZR,  Double precision arrays of length N.  (INPUT)
c  RITZI   Real and imaginary parts of the Ritz values to be checked
c          for convergence.

c  BOUNDS  Double precision array of length N.  (INPUT)
c          Ritz estimates for the Ritz values in RITZR and RITZI.
c
c  TOL     Double precision scalar.  (INPUT)
c          Desired backward error for a Ritz value to be considered
c          "converged".
c
c  NCONV   Integer scalar.  (OUTPUT)
c          Number of "converged" Ritz values.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     second  ARPACK utility routine for timing.
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c
c\SCCS Information: @(#)
c FILE: nconv.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     1. xxxx
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnconv (n, ritzr, ritzi, bounds, tol, nconv)

      implicit none
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    n, nconv
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%

      Double precision
     &           ritzr(n), ritzi(n), bounds(n)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i
      Double precision
     &           temp, eps23
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlapy2, dlamch
      external   dlapy2, dlamch

c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------------------------------------%
c     | Convergence test: unlike in the symmetric code, I am not    |
c     | using things like refined error bounds and gap condition    |
c     | because I don't know the exact equivalent concept.          |
c     |                                                             |
c     | Instead the i-th Ritz value is considered "converged" when: |
c     |                                                             |
c     |     bounds(i) .le. ( TOL * | ritz | )                       |
c     |                                                             |
c     | for some appropriate choice of norm.                        |
c     %-------------------------------------------------------------%
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch('Epsilon-Machine')
      eps23 = eps23**(2.0D+0 / 3.0D+0)
c
      nconv  = 0
      do 20 i = 1, n
         temp = max( eps23, dlapy2( ritzr(i), ritzi(i) ) )
         if (bounds(i) .le. tol*temp)   nconv = nconv + 1
   20 continue
c
      return
c
c     %---------------%
c     | End of dnconv |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsortc
c
c\Description:
c  Sorts the complex array in XREAL and XIMAG into the order
c  specified by WHICH and optionally applies the permutation to the
c  real array Y. It is assumed that if an element of XIMAG is
c  nonzero, then its negative is also an element. In other words,
c  both members of a complex conjugate pair are to be sorted and the
c  pairs are kept adjacent to each other.
c
c\Usage:
c  call dsortc
c     ( WHICH, APPLY, N, XREAL, XIMAG, Y )
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> sort XREAL,XIMAG into increasing order of magnitude.
c          'SM' -> sort XREAL,XIMAG into decreasing order of magnitude.
c          'LR' -> sort XREAL into increasing order of algebraic.
c          'SR' -> sort XREAL into decreasing order of algebraic.
c          'LI' -> sort XIMAG into increasing order of magnitude.
c          'SI' -> sort XIMAG into decreasing order of magnitude.
c          NOTE: If an element of XIMAG is non-zero, then its negative
c                is also an element.
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to array Y.
c          APPLY = .FALSE. -> do not apply the sorted order to array Y.
c
c  N       Integer.  (INPUT)
c          Size of the arrays.
c
c  XREAL,  Double precision array of length N.  (INPUT/OUTPUT)
c  XIMAG   Real and imaginary part of the array to be sorted.
c
c  Y       Double precision array of length N.  (INPUT/OUTPUT)
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c               Adapted from the sort routine in LANSO.
c
c\SCCS Information: @(#)
c FILE: sortc.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsortc (which, apply, n, xreal, ximag, y)
c
      implicit none
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  which*2
      logical    apply
      integer    n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           xreal(0:n-1), ximag(0:n-1), y(0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      Double precision
     &           temp, temp1, temp2
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlapy2
      external   dlapy2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      igap = n / 2
c
      if (which .eq. 'LM') then
c
c        %------------------------------------------------------%
c        | Sort XREAL,XIMAG into increasing order of magnitude. |
c        %------------------------------------------------------%
c
   10    continue
         if (igap .eq. 0) go to 9000
c
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            temp1 = dlapy2(xreal(j),ximag(j))
            temp2 = dlapy2(xreal(j+igap),ximag(j+igap))
c
            if (temp1.gt.temp2) then
                temp = xreal(j)
                xreal(j) = xreal(j+igap)
                xreal(j+igap) = temp
c
                temp = ximag(j)
                ximag(j) = ximag(j+igap)
                ximag(j+igap) = temp
c
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
   30    continue
         igap = igap / 2
         go to 10
c
      else if (which .eq. 'SM') then
c
c        %------------------------------------------------------%
c        | Sort XREAL,XIMAG into decreasing order of magnitude. |
c        %------------------------------------------------------%
c
   40    continue
         if (igap .eq. 0) go to 9000
c
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j .lt. 0) go to 60
c
            temp1 = dlapy2(xreal(j),ximag(j))
            temp2 = dlapy2(xreal(j+igap),ximag(j+igap))
c
            if (temp1.lt.temp2) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c
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
   60    continue
         igap = igap / 2
         go to 40
c
      else if (which .eq. 'LR') then
c
c        %------------------------------------------------%
c        | Sort XREAL into increasing order of algebraic. |
c        %------------------------------------------------%
c
   70    continue
         if (igap .eq. 0) go to 9000
c
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c
            if (xreal(j).gt.xreal(j+igap)) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c
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
   90    continue
         igap = igap / 2
         go to 70
c
      else if (which .eq. 'SR') then
c
c        %------------------------------------------------%
c        | Sort XREAL into decreasing order of algebraic. |
c        %------------------------------------------------%
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (xreal(j).lt.xreal(j+igap)) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c
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
  120    continue
         igap = igap / 2
         go to 100
c
      else if (which .eq. 'LI') then
c
c        %------------------------------------------------%
c        | Sort XIMAG into increasing order of magnitude. |
c        %------------------------------------------------%
c
  130    continue
         if (igap .eq. 0) go to 9000
         do 150 i = igap, n-1
            j = i-igap
  140       continue
c
            if (j.lt.0) go to 150
c
            if (abs(ximag(j)).gt.abs(ximag(j+igap))) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c
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
  150    continue
         igap = igap / 2
         go to 130
c
      else if (which .eq. 'SI') then
c
c        %------------------------------------------------%
c        | Sort XIMAG into decreasing order of magnitude. |
c        %------------------------------------------------%
c
  160    continue
         if (igap .eq. 0) go to 9000
         do 180 i = igap, n-1
            j = i-igap
  170       continue
c
            if (j.lt.0) go to 180
c
            if (abs(ximag(j)).lt.abs(ximag(j+igap))) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c
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
  180    continue
         igap = igap / 2
         go to 160
      end if
c
 9000 continue
      return
c
c     %---------------%
c     | End of dsortc |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dneigh
c
c\Description:
c  Compute the eigenvalues of the current upper Hessenberg matrix
c  and the corresponding Ritz estimates given the current residual norm.
c
c\Usage:
c  call dneigh
c     ( RNORM, N, H, LDH, RITZR, RITZI, BOUNDS, Q, LDQ, WORKL, IERR )
c
c\Arguments
c  RNORM   Double precision scalar.  (INPUT)
c          Residual norm corresponding to the current upper Hessenberg
c          matrix H.
c
c  N       Integer.  (INPUT)
c          Size of the matrix H.
c
c  H       Double precision N by N array.  (INPUT)
c          H contains the current upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RITZR,  Double precision arrays of length N.  (OUTPUT)
c  RITZI   On output, RITZR(1:N) (resp. RITZI(1:N)) contains the real
c          (respectively imaginary) parts of the eigenvalues of H.
c
c  BOUNDS  Double precision array of length N.  (OUTPUT)
c          On output, BOUNDS contains the Ritz estimates associated with
c          the eigenvalues RITZR and RITZI.  This is equal to RNORM
c          times the last components of the eigenvectors corresponding
c          to the eigenvalues in RITZR and RITZI.
c
c  Q       Double precision N by N array.  (WORKSPACE)
c          Workspace needed to store the eigenvectors of H.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length N**2 + 3*N.  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  This is needed to keep the full Schur form
c          of H and also in the calculation of the eigenvectors of H.
c
c  IERR    Integer.  (OUTPUT)
c          Error exit flag from dlaqrb or dtrevc.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dlaqrb  ARPACK routine to compute the real Schur form of an
c             upper Hessenberg matrix and last row of the Schur vectors.
c     second  ARPACK utility routine for timing.
c     dmout   ARPACK utility routine that prints matrices
c     dvout   ARPACK utility routine that prints vectors.
c     dlacpy  LAPACK matrix copy routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dtrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper quasi-triangular form
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dscal   Level 1 BLAS that scales a vector.
c
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c
c\SCCS Information: @(#)
c FILE: neigh.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dneigh (rnorm, n, h, ldh, ritzr, ritzi, bounds,
     &     q, ldq, workl, ierr)

      implicit none
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    ierr, n, ldh, ldq
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           bounds(n), h(ldh,n), q(ldq,n), ritzi(n), ritzr(n),
     &           workl(n*(n+3))
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
      integer(4)
     &           ione
      parameter (ione = 1)

c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
c     % dneightrue = .true. to adress logical(4)/(8)
c     % initialization for spam/spam64
c     % analogue for dneighone
      logical    select(1), dneightrue
      integer    i, iconj,  dneighone
      Double precision
     &           temp, vl(1)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dlacpy, dlaqrb, dtrevc
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlapy2, dnrm2
      external   dlapy2, dnrm2
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic  abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      dneightrue = .true.
      dneighone = 1
c
c
c     %-----------------------------------------------------------%
c     | 1. Compute the eigenvalues, the last components of the    |
c     |    corresponding Schur vectors and the full Schur form T  |
c     |    of the current upper Hessenberg matrix H.              |
c     | dlaqrb returns the full Schur form of H in WORKL(1:N**2)  |
c     | and the last components of the Schur vectors in BOUNDS.   |
c     %-----------------------------------------------------------%
c
cx      call dlacpy ('All', n, n, h, ldh, workl, n)
      call dlacpy ('All', n, n, h(1,1), ldh, workl, n)
      call dlaqrb (dneightrue, n, dneighone, n, workl, n, ritzr, ritzi,
     &             bounds, ierr)
      if (ierr .ne. 0) go to 9000
c
c
c     %-----------------------------------------------------------%
c     | 2. Compute the eigenvectors of the full Schur form T and  |
c     |    apply the last components of the Schur vectors to get  |
c     |    the last components of the corresponding eigenvectors. |
c     | Remember that if the i-th and (i+1)-st eigenvalues are    |
c     | complex conjugate pairs, then the real & imaginary part   |
c     | of the eigenvector components are split across adjacent   |
c     | columns of Q.                                             |
c     %-----------------------------------------------------------%
c
      call dtrevc ('R', 'A', select, n, workl, n, vl, n, q, ldq,
     &             n, n, workl(n*n+1), ierr)
c
      if (ierr .ne. 0) go to 9000
c
c     %------------------------------------------------%
c     | Scale the returning eigenvectors so that their |
c     | euclidean norms are all one. LAPACK subroutine |
c     | dtrevc returns each eigenvector normalized so  |
c     | that the element of largest magnitude has      |
c     | magnitude 1; here the magnitude of a complex   |
c     | number (x,y) is taken to be |x| + |y|.         |
c     %------------------------------------------------%
c
      iconj = 0
      do 10 i=1, n
         if ( abs( ritzi(i) ) .le. zero ) then
c
c           %----------------------%
c           | Real eigenvalue case |
c           %----------------------%
c
            temp = dnrm2( n, q(1,i), 1 )
            call dscal ( n, one / temp, q(1,i), ione )
         else
c
c           %-------------------------------------------%
c           | Complex conjugate pair case. Note that    |
c           | since the real and imaginary part of      |
c           | the eigenvector are stored in consecutive |
c           | columns, we further normalize by the      |
c           | square root of two.                       |
c           %-------------------------------------------%
c
            if (iconj .eq. 0) then
               temp = dlapy2( dnrm2( n, q(1,i), 1 ),
     &                        dnrm2( n, q(1,i+1), 1 ) )
               call dscal ( n, one / temp, q(1,i),
     &                      ione )
               call dscal ( n, one / temp, q(1,i+1),
     &                      ione )
               iconj = 1
            else
               iconj = 0
            end if
         end if
   10 continue
c
      call dgemv('T', n, n, one, q, ldq, bounds(1), 1, zero,
     &         workl(1), 1)
c
c
c     %----------------------------%
c     | Compute the Ritz estimates |
c     %----------------------------%
c
      iconj = 0
      do 20 i = 1, n
         if ( abs( ritzi(i) ) .le. zero ) then
c
c           %----------------------%
c           | Real eigenvalue case |
c           %----------------------%
c
            bounds(i) = rnorm * abs( workl(i) )
         else
c
c           %-------------------------------------------%
c           | Complex conjugate pair case. Note that    |
c           | since the real and imaginary part of      |
c           | the eigenvector are stored in consecutive |
c           | columns, we need to take the magnitude    |
c           | of the last components of the two vectors |
c           %-------------------------------------------%
c
            if (iconj .eq. 0) then
               bounds(i) = rnorm * dlapy2( workl(i), workl(i+1) )
               bounds(i+1) = bounds(i)
               iconj = 1
            else
               iconj = 0
            end if
         end if
   20 continue
c
 9000 continue
      return
c
c     %---------------%
c     | End of dneigh |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dlaqrb
c
c\Description:
c  Compute the eigenvalues and the Schur decomposition of an upper
c  Hessenberg submatrix in rows and columns ILO to IHI.  Only the
c  last component of the Schur vectors are computed.
c
c  This is mostly a modification of the LAPACK routine dlahqr.
c
c\Usage:
c  call dlaqrb
c     ( WANTT, N, ILO, IHI, H, LDH, WR, WI,  Z, INFO )
c
c\Arguments
c  WANTT   Logical variable.  (INPUT)
c          = .TRUE. : the full Schur form T is required;
c          = .FALSE.: only eigenvalues are required.
c
c  N       Integer.  (INPUT)
c          The order of the matrix H.  N >= 0.
c
c  ILO     Integer.  (INPUT)
c  IHI     Integer.  (INPUT)
c          It is assumed that H is already upper quasi-triangular in
c          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
c          ILO = 1). SLAQRB works primarily with the Hessenberg
c          submatrix in rows and columns ILO to IHI, but applies
c          transformations to all of H if WANTT is .TRUE..
c          1 <= ILO <= max(1,IHI); IHI <= N.
c
c  H       Double precision array, dimension (LDH,N).  (INPUT/OUTPUT)
c          On entry, the upper Hessenberg matrix H.
c          On exit, if WANTT is .TRUE., H is upper quasi-triangular in
c          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in
c          standard form. If WANTT is .FALSE., the contents of H are
c          unspecified on exit.
c
c  LDH     Integer.  (INPUT)
c          The leading dimension of the array H. LDH >= max(1,N).
c
c  WR      Double precision array, dimension (N).  (OUTPUT)
c  WI      Double precision array, dimension (N).  (OUTPUT)
c          The real and imaginary parts, respectively, of the computed
c          eigenvalues ILO to IHI are stored in the corresponding
c          elements of WR and WI. If two eigenvalues are computed as a
c          complex conjugate pair, they are stored in consecutive
c          elements of WR and WI, say the i-th and (i+1)th, with
c          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
c          eigenvalues are stored in the same order as on the diagonal
c          of the Schur form returned in H, with WR(i) = H(i,i), and, if
c          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
c          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
c
c  Z       Double precision array, dimension (N).  (OUTPUT)
c          On exit Z contains the last components of the Schur vectors.
c
c  INFO    Integer.  (OUPUT)
c          = 0: successful exit
c          > 0: SLAQRB failed to compute all the eigenvalues ILO to IHI
c               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
c               elements i+1:ihi of WR and WI contain those eigenvalues
c               which have been successfully computed.
c
c\Remarks
c  1. None.
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dlabad  LAPACK routine that computes machine constants.
c     dlamch  LAPACK routine that determines machine constants.
c     dlanhs  LAPACK routine that computes various norms of a matrix.
c     dlanv2  LAPACK routine that computes the Schur factorization of
c             2 by 2 nonsymmetric matrix in standard form.
c     dlarfg  LAPACK Householder reflection construction routine.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     drot    Level 1 BLAS that applies a rotation to a 2 by 2 matrix.

c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c               Modified from the LAPACK routine dlahqr so that only the
c               last component of the Schur vectors are computed.
c
c\SCCS Information: @(#)
c FILE: laqrb.F   SID: 2.2   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dlaqrb ( wantt, n, ilo, ihi, h, ldh, wr, wi,
     &                    z, info )
c
      implicit none
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      logical    wantt
      integer    ihi, ilo, info, ldh, n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           h( ldh, * ), wi( * ), wr( * ), z( * )
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           zero, one, dat1, dat2
      parameter (zero = 0.0D+0, one = 1.0D+0, dat1 = 7.5D-1,
     &           dat2 = -4.375D-1)
c
      integer(4) ione
      parameter (ione = 1)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      integer    i, i1, i2, itn, its, j, k, l, m, nh, nr
      integer    dlqthree
     &
      parameter (dlqthree = 3)
      Double precision
     &           cs, h00, h10, h11, h12, h21, h22, h33, h33s,
     &           h43h34, h44, h44s, ovfl, s, smlnum, sn, sum,
     &           t1, t2, t3, tst1, ulp, unfl, v1, v2, v3
      Double precision
     &           v( 3 ), work( 1 )
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
c
      Double precision
     &           dlamch, dlanhs
      external   dlamch, dlanhs
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dlabad, dlanv2, dlarfg, drot
c
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      info = 0
c
c     % if i2 not initialized; Warning: ‘i2’ may be used uninitialized
c     % in this function [-Wmaybe-uninitialized]
      i2 = 0
c
c     %--------------------------%
c     | Quick return if possible |
c     %--------------------------%
c
      if( n.eq.0 )
     &   return
      if( ilo.eq.ihi ) then
         wr( ilo ) = h( ilo, ilo )
         wi( ilo ) = zero
         return
      end if
c
c     %---------------------------------------------%
c     | Initialize the vector of last components of |
c     | the Schur vectors for accumulation.         |
c     %---------------------------------------------%
c
      do 5 j = 1, n-1
         z(j) = zero
  5   continue
      z(n) = one
c
      nh = ihi - ilo + 1
c
c     %-------------------------------------------------------------%
c     | Set machine-dependent constants for the stopping criterion. |
c     | If norm(H) <= sqrt(OVFL), overflow should not occur.        |
c     %-------------------------------------------------------------%
c
      unfl = dlamch( 'safe minimum' )
      ovfl = one / unfl
      call dlabad( unfl, ovfl )
      ulp = dlamch( 'precision' )
      smlnum = unfl*( nh / ulp )
c
c     %---------------------------------------------------------------%
c     | I1 and I2 are the indices of the first row and last column    |
c     | of H to which transformations must be applied. If eigenvalues |
c     | only are computed, I1 and I2 are set inside the main loop.    |
c     | Zero out H(J+2,J) = ZERO for J=1:N if WANTT = .TRUE.          |
c     | else H(J+2,J) for J=ILO:IHI-ILO-1 if WANTT = .FALSE.          |
c     %---------------------------------------------------------------%
c
      if( wantt ) then
         i1 = 1
         i2 = n
         do 8 i=1,i2-2
            h(i1+i+1,i) = zero
 8       continue
      else
         do 9 i=1, ihi-ilo-1
            h(ilo+i+1,ilo+i-1) = zero
 9       continue
      end if
c
c     %---------------------------------------------------%
c     | ITN is the total number of QR iterations allowed. |
c     %---------------------------------------------------%
c
      itn = 30*nh
c
c     ------------------------------------------------------------------
c     The main loop begins here. I is the loop index and decreases from
c     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
c     with the active submatrix in rows and columns L to I.
c     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
c     H(L,L-1) is negligible so that the matrix splits.
c     ------------------------------------------------------------------
c
      i = ihi
   10 continue
      l = ilo
      if( i.lt.ilo )
     &   go to 150

c     %--------------------------------------------------------------%
c     | Perform QR iterations on rows and columns ILO to I until a   |
c     | submatrix of order 1 or 2 splits off at the bottom because a |
c     | subdiagonal element has become negligible.                   |
c     %--------------------------------------------------------------%

      do 130 its = 0, itn
c
c        %----------------------------------------------%
c        | Look for a single small subdiagonal element. |
c        %----------------------------------------------%
c
         do 20 k = i, l + 1, -1
            tst1 = abs( h( k-1, k-1 ) ) + abs( h( k, k ) )
            if( abs(tst1) .le. zero)
     &         tst1 = dlanhs( '1', i-l+1, h( l, l ), ldh, work )
            if( abs( h( k, k-1 ) ).le.max( ulp*tst1, smlnum ) )
     &         go to 30
   20    continue
   30    continue
         l = k
         if( l.gt.ilo ) then
c
c           %------------------------%
c           | H(L,L-1) is negligible |
c           %------------------------%
c
            h( l, l-1 ) = zero
         end if
c
c        %-------------------------------------------------------------%
c        | Exit from loop if a submatrix of order 1 or 2 has split off |
c        %-------------------------------------------------------------%
c
         if( l.ge.i-1 )
     &      go to 140
c
c        %---------------------------------------------------------%
c        | Now the active submatrix is in rows and columns L to I. |
c        | If eigenvalues only are being computed, only the active |
c        | submatrix need be transformed.                          |
c        %---------------------------------------------------------%
c
         if( .not.wantt ) then
            i1 = l
            i2 = i
         end if
c
         if( its.eq.10 .or. its.eq.20 ) then
c
c           %-------------------%
c           | Exceptional shift |
c           %-------------------%
c
            s = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
            h44 = dat1*s
            h33 = h44
            h43h34 = dat2*s*s
c
         else
c
c           %-----------------------------------------%
c           | Prepare to use Wilkinson's double shift |
c           %-----------------------------------------%
c
            h44 = h( i, i )
            h33 = h( i-1, i-1 )
            h43h34 = h( i, i-1 )*h( i-1, i )
         end if
c
c        %-----------------------------------------------------%
c        | Look for two consecutive small subdiagonal elements |
c        %-----------------------------------------------------%
c
         do 40 m = i - 2, l, -1
c
c           %---------------------------------------------------------%
c           | Determine the effect of starting the double-shift QR    |
c           | iteration at row M, and see if this would make H(M,M-1) |
c           | negligible.                                             |
c           %---------------------------------------------------------%
c
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
            if( m.eq.l )
     &         go to 50
            h00 = h( m-1, m-1 )
            h10 = h( m, m-1 )
            tst1 = abs( v1 )*( abs( h00 )+abs( h11 )+abs( h22 ) )
            if( abs( h10 )*( abs( v2 )+abs( v3 ) ).le.ulp*tst1 )
     &         go to 50
   40    continue
   50    continue
c
c        %----------------------%
c        | Double-shift QR step |
c        %----------------------%
c
         do 120 k = m, i - 1
c
c           ------------------------------------------------------------
c           The first iteration of this loop determines a reflection G
c           from the vector V and applies it from left and right to H,
c           thus creating a nonzero bulge below the subdiagonal.
c
c           Each subsequent iteration determines a reflection G to
c           restore the Hessenberg form in the (K-1)th column, and thus
c           chases the bulge one step toward the bottom of the active
c           submatrix. NR is the order of G.
c           ------------------------------------------------------------
c
            nr = min( dlqthree, i-k+1 )
            if( k.gt.m )
     &         call dcopy( nr, h( k, k-1 ), 1, v(1), 1 )
            call dlarfg( nr, v( 1 ), v( 2 ), 1, t1 )
            if( k.gt.m ) then
               h( k, k-1 ) = v( 1 )
               h( k+1, k-1 ) = zero
               if( k.lt.i-1 )
     &            h( k+2, k-1 ) = zero
            else if( m.gt.l ) then
               h( k, k-1 ) = -h( k, k-1 )
            end if
            v2 = v( 2 )
            t2 = t1*v2
            if( nr.eq.3 ) then
               v3 = v( 3 )
               t3 = t1*v3
c
c              %------------------------------------------------%
c              | Apply G from the left to transform the rows of |
c              | the matrix in columns K to I2.                 |
c              %------------------------------------------------%
c
               do 60 j = k, i2
                  sum = h( k, j ) + v2*h( k+1, j ) + v3*h( k+2, j )
                  h( k, j ) = h( k, j ) - sum*t1
                  h( k+1, j ) = h( k+1, j ) - sum*t2
                  h( k+2, j ) = h( k+2, j ) - sum*t3
   60          continue
c
c              %----------------------------------------------------%
c              | Apply G from the right to transform the columns of |
c              | the matrix in rows I1 to min(K+3,I).               |
c              %----------------------------------------------------%
c
               do 70 j = i1, min( k+3, i )
                  sum = h( j, k ) + v2*h( j, k+1 ) + v3*h( j, k+2 )
                  h( j, k ) = h( j, k ) - sum*t1
                  h( j, k+1 ) = h( j, k+1 ) - sum*t2
                  h( j, k+2 ) = h( j, k+2 ) - sum*t3
   70          continue
c
c              %----------------------------------%
c              | Accumulate transformations for Z |
c              %----------------------------------%
c
               sum      = z( k ) + v2*z( k+1 ) + v3*z( k+2 )
               z( k )   = z( k ) - sum*t1
               z( k+1 ) = z( k+1 ) - sum*t2
               z( k+2 ) = z( k+2 ) - sum*t3

            else if( nr.eq.2 ) then
c
c              %------------------------------------------------%
c              | Apply G from the left to transform the rows of |
c              | the matrix in columns K to I2.                 |
c              %------------------------------------------------%
c
               do 90 j = k, i2
                  sum = h( k, j ) + v2*h( k+1, j )
                  h( k, j ) = h( k, j ) - sum*t1
                  h( k+1, j ) = h( k+1, j ) - sum*t2
   90          continue
c
c              %----------------------------------------------------%
c              | Apply G from the right to transform the columns of |
c              | the matrix in rows I1 to min(K+3,I).               |
c              %----------------------------------------------------%
c
               do 100 j = i1, i
                  sum = h( j, k ) + v2*h( j, k+1 )
                  h( j, k ) = h( j, k ) - sum*t1
                  h( j, k+1 ) = h( j, k+1 ) - sum*t2
  100          continue
c
c              %----------------------------------%
c              | Accumulate transformations for Z |
c              %----------------------------------%
c
               sum      = z( k ) + v2*z( k+1 )
               z( k )   = z( k ) - sum*t1
               z( k+1 ) = z( k+1 ) - sum*t2
            end if
  120    continue

  130 continue
c
c     %-------------------------------------------------------%
c     | Failure to converge in remaining number of iterations |
c     %-------------------------------------------------------%
c
      info = i
      return

  140 continue

      if( l.eq.i ) then
c
c        %------------------------------------------------------%
c        | H(I,I-1) is negligible: one eigenvalue has converged |
c        %------------------------------------------------------%
c
         wr( i ) = h( i, i )
         wi( i ) = zero

      else if( l.eq.i-1 ) then
c
c        %--------------------------------------------------------%
c        | H(I-1,I-2) is negligible;                              |
c        | a pair of eigenvalues have converged.                  |
c        |                                                        |
c        | Transform the 2-by-2 submatrix to standard Schur form, |
c        | and compute and store the eigenvalues.                 |
c        %--------------------------------------------------------%
c
         call dlanv2( h( i-1, i-1 ), h( i-1, i ), h( i, i-1 ),
     &                h( i, i ), wr( i-1 ), wi( i-1 ), wr( i ), wi( i ),
     &                cs, sn )

         if( wantt ) then
c
c           %-----------------------------------------------------%
c           | Apply the transformation to the rest of H and to Z, |
c           | as required.                                        |
c           %-----------------------------------------------------%
c
            if( i2.gt.i )
     &         call drot( i2-i, h( i-1, i+1 ), ldh, h( i, i+1 ), ldh,
     &                    cs, sn )
            call drot( i-i1-1, h( i1, i-1 ), ione, h( i1, i ), ione,
     &                 cs, sn )
            sum      = cs*z( i-1 ) + sn*z( i )
            z( i )   = cs*z( i )   - sn*z( i-1 )
            z( i-1 ) = sum
         end if
      end if
c
c     %---------------------------------------------------------%
c     | Decrement number of remaining iterations, and return to |
c     | start of the main loop with new value of I.             |
c     %---------------------------------------------------------%
c
      itn = itn - its
      i = l - 1
      go to 10

  150 continue
      return
c
c     %---------------%
c     | End of dlaqrb |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dngets
c
c\Description:
c  Given the eigenvalues of the upper Hessenberg matrix H,
c  computes the NP shifts AMU that are zeros of the polynomial of
c  degree NP which filters out components of the unwanted eigenvectors
c  corresponding to the AMU's based on some given criteria.
c
c  NOTE: call this even in the case of user specified shifts in order
c  to sort the eigenvalues, and error bounds of H for later use.
c
c\Usage:
c  call dngets
c     ( ISHIFT, WHICH, KEV, NP, RITZR, RITZI, BOUNDS, SHIFTR, SHIFTI )
c
c\Arguments
c  ISHIFT  Integer.  (INPUT)
c          Method for selecting the implicit shifts at each iteration.
c          ISHIFT = 0: user specified shifts
c          ISHIFT = 1: exact shift with respect to the matrix H.
c
c  WHICH   Character*2.  (INPUT)
c          Shift selection criteria.
c          'LM' -> want the KEV eigenvalues of largest magnitude.
c          'SM' -> want the KEV eigenvalues of smallest magnitude.
c          'LR' -> want the KEV eigenvalues of largest real part.
c          'SR' -> want the KEV eigenvalues of smallest real part.
c          'LI' -> want the KEV eigenvalues of largest imaginary part.
c          'SI' -> want the KEV eigenvalues of smallest imaginary part.
c
c  KEV      Integer.  (INPUT/OUTPUT)
c           INPUT: KEV+NP is the size of the matrix H.
c           OUTPUT: Possibly increases KEV by one to keep complex conjugate
c           pairs together.
c
c  NP       Integer.  (INPUT/OUTPUT)
c           Number of implicit shifts to be computed.
c           OUTPUT: Possibly decreases NP by one to keep complex conjugate
c           pairs together.
c
c  RITZR,  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c  RITZI   On INPUT, RITZR and RITZI contain the real and imaginary
c          parts of the eigenvalues of H.
c          On OUTPUT, RITZR and RITZI are sorted so that the unwanted
c          eigenvalues are in the first NP locations and the wanted
c          portion is in the last KEV locations.  When exact shifts are
c          selected, the unwanted part corresponds to the shifts to
c          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
c          are further sorted so that the ones with largest Ritz values
c          are first.
c
c  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.
c
c  SHIFTR, SHIFTI  *** USE deprecated as of version 2.1. ***
c
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dsortc  ARPACK sorting routine.
c     dcopy   Level 1 BLAS that copies one vector to another .
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c
c\SCCS Information: @(#)
c FILE: ngets.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     1. xxxx
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dngets ( ishift, which, kev, np, ritzr, ritzi, bounds)
c     &                    shiftr, shifti )
c
      implicit none
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  which*2
      integer    ishift, kev, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           bounds(kev+np), ritzr(kev+np), ritzi(kev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           zero
      parameter (zero = 0.0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      logical    dngetstrue
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dsortc
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      dngetstrue = .true.
c
c
c     %----------------------------------------------------%
c     | LM, SM, LR, SR, LI, SI case.                       |
c     | Sort the eigenvalues of H into the desired order   |
c     | and apply the resulting order to BOUNDS.           |
c     | The eigenvalues are sorted so that the wanted part |
c     | are always in the last KEV locations.              |
c     | We first do a pre-processing sort in order to keep |
c     | complex conjugate pairs together                   |
c     %----------------------------------------------------%
c
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
c
      call dsortc (which, dngetstrue, kev+np, ritzr, ritzi, bounds)
c
c     %-------------------------------------------------------%
c     | Increase KEV by one if the ( ritzr(np),ritzi(np) )    |
c     | = ( ritzr(np+1),-ritzi(np+1) ) and ritz(np) .ne. zero |
c     | Accordingly decrease NP by one. In other words keep   |
c     | complex conjugate pairs together.                     |
c     %-------------------------------------------------------%
c
      if (       (  abs(ritzr(np+1) - ritzr(np)) .le. zero)
     &     .and. (  abs(ritzi(np+1) + ritzi(np)) .le. zero) ) then

         np = np - 1
         kev = kev + 1
      end if
c
      if ( ishift .eq. 1 ) then
c
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first        |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when they shifts |
c        | are applied in subroutine dnapps.                     |
c        | Be careful and use 'SR' since we want to sort BOUNDS! |
c        %-------------------------------------------------------%
c
         call dsortc ( 'SR', dngetstrue, np, bounds, ritzr, ritzi )
      end if
c
      return
c
c     %---------------%
c     | End of dngets |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dnaitr
c
c\Description:
c  Reverse communication interface for applying NP additional steps to
c  a K step nonsymmetric Arnoldi factorization.
c
c  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
c
c          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
c
c  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
c
c          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
c
c  where OP and B are as in dnaupd.  The B-norm of r_{k+p} is also
c  computed and returned.
c
c\Usage:
c  call dnaitr
c     ( IDO, BMAT, N, K, NP, NB, RESID, RNORM, V, LDV, H, LDH,
c       IPNTR, WORKD, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c                    This is for the restart phase to force the new
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y,
c                    IPNTR(3) is the pointer into WORK for B * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c          When the routine is used in the "shift-and-invert" mode, the
c          vector B * Q is already available and do not need to be
c          recompute in forming OP * Q.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.  See dnaupd.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*M**x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  K       Integer.  (INPUT)
c          Current size of V and H.
c
c  NP      Integer.  (INPUT)
c          Number of additional Arnoldi steps to take.
c
c  NB      Integer.  (INPUT)
c          Blocksize to be used in the recurrence.
c          Only work for NB = 1 right now.  The goal is to have a
c          program that implement both the block and non-block method.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT:  RESID contains the residual vector r_{k}.
c          On OUTPUT: RESID contains the residual vector r_{k+p}.
c
c  RNORM   Double precision scalar.  (INPUT/OUTPUT)
c          B-norm of the starting residual on input.
c          B-norm of the updated residual r_{k+p} on output.
c
c  V       Double precision N by K+NP array.  (INPUT/OUTPUT)
c          On INPUT:  V contains the Arnoldi vectors in the first K
c          columns.
c          On OUTPUT: V contains the new NP Arnoldi vectors in the next
c          NP columns.  The first K columns are unchanged.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Double precision (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
c          H is used to store the generated upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORK for
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The calling program should not
c          use WORKD as temporary workspace during the iteration !!!!!!
c          On input, WORKD(1:N) = B*RESID and is used to save some
c          computation at the first step.
c
c  INFO    Integer.  (OUTPUT)
c          = 0: Normal exit.
c          > 0: Size of the spanning invariant subspace of OP found.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     dgetv0  ARPACK routine to generate the initial vector.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dmout   ARPACK utility routine that prints matrices
c     dvout   ARPACK utility routine that prints vectors.
c     dlabad  LAPACK routine that computes machine constants.
c     dlamch  LAPACK routine that determines machine constants.
c     dlascl  LAPACK routine for careful scaling of a matrix.
c     dlanhs  LAPACK routine that computes various norms of a matrix.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dscal   Level 1 BLAS that scales a vector.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     ddot    Level 1 BLAS that computes the scalar product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c
c\SCCS Information: @(#)
c FILE: naitr.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c  The algorithm implemented is:
c
c  restart = .false.
c  Given V_{k} = [v_{1}, ..., v_{k}], r_{k};
c  r_{k} contains the initial residual vector even for k = 0;
c  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already
c  computed by the calling program.
c
c  betaj = rnorm ; p_{k+1} = B*r_{k} ;
c  For  j = k+1, ..., k+np  Do
c     1) if ( betaj < tol ) stop or restart depending on j.
c        ( At present tol is zero )
c        if ( restart ) generate a new starting vector.
c     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];
c        p_{j} = p_{j}/betaj
c     3) r_{j} = OP*v_{j} where OP is defined as in dnaupd
c        For shift-invert mode p_{j} = B*v_{j} is already available.
c        wnorm = || OP*v_{j} ||
c     4) Compute the j-th step residual vector.
c        w_{j} =  V_{j}^T * B * OP * v_{j}
c        r_{j} =  OP*v_{j} - V_{j} * w_{j}
c        H(:,j) = w_{j};
c        H(j,j-1) = rnorm
c        rnorm = || r_(j) ||
c        If (rnorm > 0.717*wnorm) accept step and go back to 1)
c     5) Re-orthogonalization step:
c        s = V_{j}'*B*r_{j}
c        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
c        alphaj = alphaj + s_{j};
c     6) Iterative refinement step:
c        If (rnorm1 > 0.717*rnorm) then
c           rnorm = rnorm1
c           accept step and go back to 1)
c        Else
c           rnorm = rnorm1
c           If this is the first time in step 6), go to 5)
c           Else r_{j} lies in the span of V_{j} numerically.
c              Set r_{j} = 0 and rnorm = 0; go to 1)
c        EndIf
c  End Do
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnaitr
     &   (ido, bmat, n, k, np, nb, resid, rnorm, v, ldv, h, ldh,
     &     ipntr, workd, info)

      implicit none
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, nb, np
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Double precision
     &           h(ldh,k+np), resid(n), v(ldv,k+np), workd(3*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
      logical    fls
      parameter (fls = .false. )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      logical    first, orth1, orth2, rstart, step3, step4
      integer    ierr, i, infol, ipj, irj, ivj, iter, itry, j,
     &            jj
      integer    dnaitrone
      parameter (dnaitrone = 1)
c
      integer(4)
     &           ione
      parameter (ione = 1)
c
      Double precision
     &           betaj, ovfl, temp1, rnorm1, smlnum, tst1, ulp, unfl,
     &           wnorm
      save       first, orth1, orth2, rstart, step3, step4,
     &           ierr, ipj, irj, ivj, iter, itry, j,  ovfl,
     &           betaj, rnorm1, smlnum, ulp, unfl, wnorm
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   daxpy, dcopy, dscal, dgemv, dgetv0, dlabad
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2, dlanhs, dlamch
      external   ddot, dnrm2, dlanhs, dlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, sqrt
c
c     %-----------------%
c     | Data statements |
c     %-----------------%
c
      data      first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (nb .gt. 1000) then
        goto 9000
      end if
c
      if (first) then
c
c        %-----------------------------------------%
c        | Set machine-dependent constants for the |
c        | the splitting and deflation criterion.  |
c        | If norm(H) <= sqrt(OVFL),               |
c        | overflow should not occur.              |
c        | REFERENCE: LAPACK subroutine dlahqr     |
c        %-----------------------------------------%
c
         unfl = dlamch( 'safe minimum' )
         ovfl = one / unfl
         call dlabad( unfl, ovfl )
         ulp = dlamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
c
      if (ido .eq. 0) then
c
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
c         call second (t0)
cm         msglvl = mnaitr
c
c        %------------------------------%
c        | Initial call to this routine |
c        %------------------------------%
c
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
c
c     %-------------------------------------------------%
c     | When in reverse communication mode one of:      |
c     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
c     | will be .true. when ....                        |
c     | STEP3: return from computing OP*v_{j}.          |
c     | STEP4: return from computing B-norm of OP*v_{j} |
c     | ORTH1: return from computing B-norm of r_{j+1}  |
c     | ORTH2: return from computing B-norm of          |
c     |        correction to the residual vector.       |
c     | RSTART: return from OP computations needed by   |
c     |         dgetv0.                                 |
c     %-------------------------------------------------%
c
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
c
c     %-----------------------------%
c     | Else this is the first step |
c     %-----------------------------%
c
c     %--------------------------------------------------------------%
c     |                                                              |
c     |        A R N O L D I     I T E R A T I O N     L O O P       |
c     |                                                              |
c     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
c     %--------------------------------------------------------------%

 1000 continue
c
c        %---------------------------------------------------%
c        | STEP 1: Check if the B norm of j-th residual      |
c        | vector is zero. Equivalent to determing whether   |
c        | an exact j-step Arnoldi factorization is present. |
c        %---------------------------------------------------%
c
         betaj = rnorm
         if (rnorm .gt. zero) go to 40
c
c           %---------------------------------------------------%
c           | Invariant subspace found, generate a new starting |
c           | vector which is orthogonal to the current Arnoldi |
c           | basis and continue the iteration.                 |
c           %---------------------------------------------------%
c           %---------------------------------------------%
c           | ITRY is the loop variable that controls the |
c           | maximum amount of times that a restart is   |
c           | attempted. NRSTRT is used by stat.h         |
c           %---------------------------------------------%
c
            betaj  = zero
cp            nrstrt = nrstrt + 1
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
c
c           %--------------------------------------%
c           | If in reverse communication mode and |
c           | RSTART = .true. flow returns here.   |
c           %--------------------------------------%
c
            call dgetv0 (ido, bmat, fls, n, j, v, ldv,
     &                   resid, rnorm, ipntr, workd, ierr)
            if (ido .ne. 99) go to 9000
            if (ierr .lt. 0) then
               itry = itry + 1
               if (itry .le. 3) go to 20
c
c              %------------------------------------------------%
c              | Give up after several restart attempts.        |
c              | Set INFO to the size of the invariant subspace |
c              | which spans OP and exit.                       |
c              %------------------------------------------------%
c
               info = j - 1
               ido = 99
               go to 9000
            end if
c
   40    continue
c
c        %---------------------------------------------------------%
c        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
c        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
c        | when reciprocating a small RNORM, test against lower    |
c        | machine bound.                                          |
c        %---------------------------------------------------------%
c
         call dcopy (n, resid, 1, v(1,j), 1)
         if (rnorm .ge. unfl) then
             temp1 = one / rnorm
             call dscal (n, temp1, v(1,j), ione)
             call dscal (n, temp1, workd(ipj), ione)
         else
c
c            %-----------------------------------------%
c            | To scale both v_{j} and p_{j} carefully |
c            | use LAPACK routine SLASCL               |
c            %-----------------------------------------%
c
             call dlascl ('General', i, i, rnorm, one, n, 1,
     &                    v(1,j), n, infol)
             call dlascl ('General', i, i, rnorm, one, n, 1,
     &                    workd(ipj), n, infol)
         end if
c
c        %------------------------------------------------------%
c        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
c        | Note that this is not quite yet r_{j}. See STEP 4    |
c        %------------------------------------------------------%
c
         step3 = .true.
         call dcopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
c
c        %-----------------------------------%
c        | Exit in order to compute OP*v_{j} |
c        %-----------------------------------%
c
         go to 9000
   50    continue
c
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   |
c        | if step3 = .true.                |
c        %----------------------------------%
c
         step3 = .false.
c
c        %------------------------------------------%
c        | Put another copy of OP*v_{j} into RESID. |
c        %------------------------------------------%
c
         call dcopy (n, workd(irj), 1, resid(1), 1)
c
c        %---------------------------------------%
c        | STEP 4:  Finish extending the Arnoldi |
c        |          factorization to length j.   |
c        %---------------------------------------%
c
         if (bmat .eq. 'G') then
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c
c           %-------------------------------------%
c           | Exit in order to compute B*OP*v_{j} |
c           %-------------------------------------%
c
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   60    continue
c
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} |
c        | if step4 = .true.                |
c        %----------------------------------%
c
         step4 = .false.
c
c        %-------------------------------------%
c        | The following is needed for STEP 5. |
c        | Compute the B-norm of OP*v_{j}.     |
c        %-------------------------------------%
c
         if (bmat .eq. 'G') then
             wnorm = ddot (n, resid, 1, workd(ipj), 1)
             wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'I') then
            wnorm = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------%
c        | Compute the j-th residual corresponding |
c        | to the j step factorization.            |
c        | Use Classical Gram Schmidt and compute: |
c        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
c        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
c        %-----------------------------------------%
c
c
c        %------------------------------------------%
c        | Compute the j Fourier coefficients w_{j} |
c        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
c        %------------------------------------------%
c
         call dgemv ('T', n, j, one, v, ldv, workd(ipj), 1,
     &               zero, h(1,j), 1)
c
c        %--------------------------------------%
c        | Orthogonalize r_{j} against V_{j}.   |
c        | RESID contains OP*v_{j}. See STEP 3. |
c        %--------------------------------------%
c
         call dgemv ('N', n, j, -one, v, ldv, h(1,j), 1,
     &               one, resid(1), 1)
c
         if (j .gt. 1) h(j,j-1) = betaj
c
         orth1 = .true.
c
         if (bmat .eq. 'G') then
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c
c           %----------------------------------%
c           | Exit in order to compute B*r_{j} |
c           %----------------------------------%
c
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   70    continue
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH1 = .true. |
c        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
c        %---------------------------------------------------%
c
         orth1 = .false.
c
c        %------------------------------%
c        | Compute the B-norm of r_{j}. |
c        %------------------------------%
c
         if (bmat .eq. 'G') then
            rnorm = ddot (n, resid, 1, workd(ipj), 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------------------------%
c        | STEP 5: Re-orthogonalization / Iterative refinement phase |
c        | Maximum NITER_ITREF tries.                                |
c        |                                                           |
c        |          s      = V_{j}^T * B * r_{j}                     |
c        |          r_{j}  = r_{j} - V_{j}*s                         |
c        |          alphaj = alphaj + s_{j}                          |
c        |                                                           |
c        | The stopping criteria used for iterative refinement is    |
c        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
c        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
c        | Determine if we need to correct the residual. The goal is |
c        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
c        | The following test determines whether the sine of the     |
c        | angle between  OP*x and the computed residual is less     |
c        | than or equal to 0.717.                                   |
c        %-----------------------------------------------------------%
c
         if (rnorm .gt. 0.717*wnorm) go to 100
         iter  = 0
c
c        %---------------------------------------------------%
c        | Enter the Iterative refinement phase. If further  |
c        | refinement is necessary, loop back here. The loop |
c        | variable is ITER. Perform a step of Classical     |
c        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
c        %---------------------------------------------------%
c
   80    continue
c
c        %----------------------------------------------------%
c        | Compute V_{j}^T * B * r_{j}.                       |
c        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
c        %----------------------------------------------------%
c
         call dgemv ('T', n, j, one, v, ldv, workd(ipj), 1,
     &               zero, workd(irj), 1)
c
c        %---------------------------------------------%
c        | Compute the correction to the residual:     |
c        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
c        | The correction to H is v(:,1:J)*H(1:J,1:J)  |
c        | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
c        %---------------------------------------------%
c
         call dgemv ('N', n, j, -one, v, ldv, workd(irj), 1,
     &               one, resid(1), 1)
         call daxpy (j, one, workd(irj), 1, h(1,j), 1)
c
         orth2 = .true.
         if (bmat .eq. 'G') then
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c
c           %-----------------------------------%
c           | Exit in order to compute B*r_{j}. |
c           | r_{j} is the corrected residual.  |
c           %-----------------------------------%
c
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   90    continue
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH2 = .true. |
c        %---------------------------------------------------%
c
c        %-----------------------------------------------------%
c        | Compute the B-norm of the corrected residual r_{j}. |
c        %-----------------------------------------------------%
c
         if (bmat .eq. 'G') then
             rnorm1 = ddot (n, resid, 1, workd(ipj), 1)
             rnorm1 = sqrt(abs(rnorm1))
         else if (bmat .eq. 'I') then
             rnorm1 = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------%
c        | Determine if we need to perform another |
c        | step of re-orthogonalization.           |
c        %-----------------------------------------%
c
         if (rnorm1 .gt. 0.717*rnorm) then
c
c           %---------------------------------------%
c           | No need for further refinement.       |
c           | The cosine of the angle between the   |
c           | corrected residual vector and the old |
c           | residual vector is greater than 0.717 |
c           | In other words the corrected residual |
c           | and the old residual vector share an  |
c           | angle of less than arcCOS(0.717)      |
c           %---------------------------------------%
c
            rnorm = rnorm1
c
         else
c
c           %-------------------------------------------%
c           | Another step of iterative refinement step |
c           | is required. NITREF is used by stat.h     |
c           %-------------------------------------------%
c
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter .le. 1) go to 80
c
c           %-------------------------------------------------%
c           | Otherwise RESID is numerically in the span of V |
c           %-------------------------------------------------%
c
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue
            rnorm = zero
         end if
c
c        %----------------------------------------------%
c        | Branch here directly if iterative refinement |
c        | wasn't necessary or after at most NITER_REF  |
c        | steps of iterative refinement.               |
c        %----------------------------------------------%
c
  100    continue
c
         rstart = .false.
         orth2  = .false.
c
c        %------------------------------------%
c        | STEP 6: Update  j = j+1;  Continue |
c        %------------------------------------%
c
         j = j + 1
         if (j .gt. k+np) then
            ido = 99
            do 110 i = max(dnaitrone,k), k+np-1
c
c              %--------------------------------------------%
c              | Check for splitting and deflation.         |
c              | Use a standard test as in the QR algorithm |
c              | REFERENCE: LAPACK subroutine dlahqr        |
c              %--------------------------------------------%
c
               tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
               if( abs(tst1) .le. zero )
     &              tst1 = dlanhs( '1', k+np, h, ldh, workd(n+1) )
               if( abs( h( i+1,i ) ).le.max( ulp*tst1, smlnum ) )
     &              h(i+1,i) = zero
 110        continue
c
            go to 9000
         end if
c
c        %--------------------------------------------------------%
c        | Loop back to extend the factorization by another step. |
c        %--------------------------------------------------------%
c
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 9000 continue
      return
c
c     %---------------%
c     | End of dnaitr |
c     %---------------%
c
      end
c
