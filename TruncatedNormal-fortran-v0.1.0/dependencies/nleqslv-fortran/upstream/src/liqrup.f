      subroutine liqrup(Q,ldq,n,R,ldr,u,v,wk)
      integer ldq,n,ldr
      double precision Q(ldq,*),R(ldr,*),u(*),v(*),wk(*)

c-----------------------------------------------------------------------------
c
c     Arguments
c
c     Inout  Q       Real(ldq,n)      orthogonal matrix from QR
c     In     ldq     Integer          leading dimension of Q
c     In     n       Integer          order of Q and R
c     Inout  R       Real(ldr,n)      upper triangular matrix R from QR
c     In     ldr     Integer          leading dimension of R
c     In     u       Real(*)          vector u of size n
c     In     v       Real(*)          vector v of size n
c     Out    wk      Real(*)          workspace of size n
c
c     on return
c
c        Q       Q is the matrix with orthonormal columns in a QR
c                decomposition of the matrix B = A + u*v'
c
c        R       R is the upper triangular matrix in a QR
c                decomposition of the matrix B = A + u*v'
c
c     Description
c
c     The matrices Q and R are a QR decomposition of a square matrix
c     A = Q*R.
c     Given Q and R, qrupdt computes a QR decomposition of the rank one
c     modification B = A + u*trans(v) of A. Here u and v are vectors.
c
c     Code based on Reichel and Cragg and simplified and modernized
c           Reichel, L. and W.B. Gragg (1990), Algorithm 686: FORTRAN subroutines for updating the QR decomposition,
c           ACM Trans. Math. Softw., 16, 4, pp. 369--377.
c
c-----------------------------------------------------------------------------

c     Local variables and functions

      integer k,i
      double precision  ddot, c, s

c     calculate wk = trans(Q)*u

      do i=1,n
         wk(i) = ddot(n,Q(1,i),1,u,1)
      enddo

c     nuvgiv uses Lapack dlartg and
c     sets its first argument x to value of the last argument returned by dlartg.

c     zero components wk(n),wk(n-1)...wk(2)
c     and apply rotators to R and Q.

      do k=n-1,1,-1
         call nuvgiv(wk(k),wk(k+1),c,s)
         call drot(n-k+1,R(k,k),ldr,R(k+1,k),ldr,c,s)
         call drot(n    ,Q(1,k),1  ,Q(1,k+1),1  ,c,s)
      enddo

c     r(1,1:n) += wk(1)*v(1:n)
      call daxpy(n,wk(1),v,1,R(1,1),ldr)

c     R is of upper hessenberg form. Triangularize R.

      do k=1,n-1
         call nuvgiv(R(k,k),R(k+1,k),c,s)
         call drot(n-k,R(k,k+1),ldr,R(k+1,k+1),ldr,c,s)
         call drot(n  ,Q(1,k)  ,1  ,Q(1,k+1)  ,1  ,c,s)
      enddo

      return
      end
