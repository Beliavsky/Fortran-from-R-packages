c
      subroutine d_ope ( n, x, y, a, ja, ia )
c
c ope computes A * x for a sparse matrix A.
c
      implicit none
c
      integer n
c
      double precision  a(*)
      integer i
      integer ia(n+1)
      integer ja(*)
      integer k1
      integer k2
      double precision  x(*)
      double precision  y(*)
c
c spasrse matrix * vector multiplication
c
      do i=1,n
         k1 = ia(i)
         k2 = ia(i+1) -1
         y(i) = dot_product ( a(k1:k2), x(ja(k1:k2)) )
      end do
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine ds_eigen_f (maxnev, ncv, maxitr,
     &                       n, iwhich,
     &                       na, a, ja, ia,
     &                       v, d,
     &                       iparam)
c
      implicit none
c
c     %--------------------%
c     | Input Declarations |
c     %--------------------%
c
      integer          maxnev, ncv, maxitr, n, na, ja(*), ia(na+1),
     &                 iparam(8), iwhich
c
      Double precision
     &                 a(*),
     &                 v(n, ncv), d(maxnev),
     &                 workl(ncv*(ncv+8)), workd(3*n), resid(n)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      Double precision
     &                 tol, sigma
      integer          ipntr(11)
c
      character        bmat*1, which*2
      integer          ido, info, lworkl,
     &                 ishfts, mode1, ierr
      logical          rvec, select(ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                 zero
      parameter        (zero = 0.0D+0)
c
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------------------------%
c     | The following include statement and assignments |
c     | initiate trace output from the internal         |
c     | actions of ARPACK.  See debug.doc in the        |
c     | DOCUMENTS directory for usage.  Initially, the  |
c     | most useful information will be a breakdown of  |
c     | time spent in the various stages of computation |
c     | given by setting msaupd = 1.                    |
c     %-------------------------------------------------%
c
c     %-------------------------------------------------%
c     | The following sets dimensions for this problem. |
c     %-------------------------------------------------%
c
      bmat  = 'I'
c
      if (iwhich .eq. 1) then
            which = 'LM'
      else if (iwhich .eq. 2) then
            which = 'SM'
      else if (iwhich .eq. 7) then
            which = 'LA'
      else if (iwhich .eq. 8) then
            which = 'SA'
      else if (iwhich .eq. 9) then
            which = 'BE'
      else
c
        goto 9000
      end if
c
c     %-----------------------------------------------------%
c     |                                                     |
c     | Specification of stopping rules and initial         |
c     | conditions before calling DSAUPD                    |
c     |                                                     |
c     | TOL  determines the stopping criterion.             |
c     |                                                     |
c     |      Expect                                         |
c     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |           (machine precision) is used.              |
c     |                                                     |
c     | IDO  is the REVERSE COMMUNICATION parameter         |
c     |      used to specify actions to be taken on return  |
c     |      from DSAUPD. (See usage below.)                |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to DSAUPD.                                |
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     |
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID.)  |
c     |                                                     |
c     | The work array WORKL is used in DSAUPD as           |
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     |                                                     |
c     %-----------------------------------------------------%
c
      lworkl = ncv * (ncv + 8)
      tol = zero
      info = 0
      ido = 0
c
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting PARAM(1) = 1).              |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of DSAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | DSAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      mode1 = 1
c
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode1
c
c     %------------------------------------------------%
c     | M A I N   L O O P (Reverse communication loop) |
c     %------------------------------------------------%
c
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DSAUPD and take |
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
c
         call dsaupd ( ido, bmat, n, which, maxnev, tol, resid,
     &                 ncv, v, n, iparam, ipntr, workd, workl,
     &                 lworkl, info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %--------------------------------------%
c           | Perform matrix vector multiplication |
c           |              y <--- OP*x             |
c           | The user should supply his/her own   |
c           | matrix vector multiplication routine |
c           | here that takes workd(ipntr(1)) as   |
c           | the input, and return the result to  |
c           | workd(ipntr(2)).                     |
c           %--------------------------------------%
c
         call d_ope (na, workd(ipntr(1)), workd(ipntr(2)),
     &               a, ja, ia)
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DSAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         end if
c
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message. Check the |
c        | documentation in DSAUPD. |
c        %--------------------------%
c
         goto 9000
c
      else
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using DSEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |
c        |                                           |
c        | Eigenvectors may be also computed now if  |
c        | desired.  (indicated by rvec = .true.)    |
c        |                                           |
c        | The routine DSEUPD now called to do this  |
c        | post processing (Other modes may require  |
c        | more complicated post processing than     |
c        | mode1.)                                   |
c        |                                           |
c        %-------------------------------------------%
c
          rvec = .true.
c
          call dseupd ( rvec, 'A', select, d, v, n, sigma,
     &         bmat, n, which, maxnev, tol, resid, ncv, v, n,
     &         iparam, ipntr, workd, workl, lworkl, ierr )
c
c         %----------------------------------------------%
c         | Eigenvalues are returned in the first column |
c         | of the two dimensional array D and the       |
c         | corresponding eigenvectors are returned in   |
c         | the first NCONV (=IPARAM(5)) columns of the  |
c         | two dimensional array V if requested.        |
c         | Otherwise, an orthogonal basis for the       |
c         | invariant subspace corresponding to the      |
c         | eigenvalues in D is returned in V.           |
c         %----------------------------------------------%
c
          if ( ierr .ne. 0) then
c
c            %------------------------------------%
c            | Error condition:                   |
c            | Check the documentation of DSEUPD. |
c            %------------------------------------%
c
             goto 9000
c
         end if
c
      end if
c
 9000 continue

      end
c
c
