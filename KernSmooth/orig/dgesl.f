c     A LINPACK routine updated for Fortran 90.

      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(*),job
      double precision a(lda,*),b(*)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c           end do
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) then
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .ge. 1) then
         do k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .ne. k) then
               b(l) = b(k)
               b(k) = t
            end if
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
         end do
      end if
c
c        now solve  u*x = y
c
         do kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
         end do
         return
      end if
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
         end do
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) return
         do kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .ne. k) then
               t = b(l)
               b(l) = b(k)
               b(k) = t
            end if
         end do
      return
      end
