c-----------------------------------------------------------------------

      subroutine lirslv(R,ldr,n,cndtol, stepadj,
     *                  qtf,dn,ierr,rcond, rcdwrk,icdwrk)

      integer ldr,n,ierr
      double precision  cndtol,R(ldr,*)
      double precision  dn(*),qtf(*)
      double precision  rcdwrk(*)
      integer           icdwrk(*)
      double precision  rcond
      logical           stepadj

c-----------------------------------------------------------------------
c
c     Solve R * dn = -qtf safely with quard against ill-conditioning
c
c     Arguments
c
c     Inout    R       Real(ldr,*)     upper triangular matrix from QR
c                                      if ill conditioned result from liqrev
c                                      Levenberg-Marquardt correction
c                                      Warning: lower triangular part of R destroyed
c     In       ldr     Integer         leading dimension of R
c     In       n       Integer         dimension of problem
c     In       cndtol  Real            tolerance of test for ill conditioning
c     In       stepadj Logical         allow adjusting step for singular/illconditioned jacobian
c     In       qtf     Real(*)         trans(Q)*f()
c     Out      dn      Real(*)         Newton direction
c     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
c                                      1 indicating Jacobian ill-conditioned
c                                      2 indicating Jacobian completely singular
c                                      3 indicating almost zero LM correction
c     Out      rcond   Real            inverse condition of upper triangular R of QR
c     Wk       rcdwrk  Real(*)         workspace
c     Wk       icdwrk  Integer(*)      workspace
c
c-----------------------------------------------------------------------

      integer k

      double precision Rone
      parameter(Rone=1.0d0)
      double precision mu

c     check for singularity or ill conditioning

      call cndjac(n,R,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)

      if( ierr .eq. 0 ) then
c         Normal Newton step
c         solve Jacobian*dn  =  -fn
c         ==> R*dn = - qtf

          call dcopy(n,qtf,1,dn,1)
          call dtrsv('U','N','N',n,r,ldr,dn,1)
          call dscal(n, -Rone, dn, 1)

      elseif( stepadj ) then
c         Adjusted Newton step
c         approximately from pseudoinverse(Jac+)
c         use mu to solve (trans(R)*R + mu*I*mu*I) * x = - trans(R) * fn
c         directly from the QR decomposition of R stacked with mu*I
c         a la Levenberg-Marquardt
          call compmu(R,ldr,n,mu,rcdwrk,ierr)
          if( ierr .eq. 0 ) then
             call liqrev(n,R,ldr,mu,qtf,dn,
     *                   rcdwrk(1+n),rcdwrk(2*n+1))
             call dscal(n, -Rone, dn, 1)

c            copy lower triangular R to upper triangular
             do k=1,n
                call dcopy (n-k+1,R(k,k),1,R(k,k),ldr)
                R(k,k) = rcdwrk(1+n+k-1)
             enddo
          endif
      endif

      return
      end
