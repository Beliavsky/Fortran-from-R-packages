C--- Common to "Original" (Old)         matrexpO() -- ./matrexpO.f
C--- and to new potentially BLAS-based  matrexp()  -- ./matrexp.f

      subroutine subtract(m,x,y,z)

c	subtracting a matrix y from a matrix x

      integer m
      double precision x(m,m),y(m,m),z(m,m)
      integer i,j

      do i=1,m
        do j=1,m
           z(i,j)=x(i,j)-y(i,j)
        enddo
      enddo

      return
      end

      subroutine addtodiag(m,x,s)

c	add a scalar s to the main diagonal elements of a matrix x

      integer m
      double precision x(m,m), s

      integer i
      do i=1,m
         x(i,i)=x(i,i)+s
      enddo
      return
      end

      subroutine minus(m,x,y)

c	the minus of a matrix

      integer m
      double precision x(m,m),y(m,m)

      integer i,j
      do i=1,m
         do j=1,m
            y(i,j)=-x(i,j)
         enddo
      enddo
      return
      end


      double precision function dl1norm(m,x)

c	L_1 norm of a matrix :=  max_j  sum_i |x_{ij}|

      integer m
      double precision x(m,m)

      double precision sum
      integer j,i

      dl1norm=0.d0
      do i=1,m
        sum=0.d0
        do j=1,m
           sum=sum+dabs(x(j,i))
        enddo
        if(sum.gt.dl1norm) dl1norm=sum
      enddo
      return
      end

      double precision function dlinfnorm(m,x)

c	L_infty norm of a matrix :=  max_i  sum_j |x_{ij}|

      integer m
      double precision x(m,m)

      double precision sum
      integer i,j

      dlinfnorm=0.d0
      do i=1,m
         sum=0.d0
         do j=1,m
            sum=sum+dabs(x(i,j))
         enddo
         if(sum.gt.dlinfnorm) dlinfnorm=sum
      enddo
      return
      end


      subroutine zero(m,x)

c	zeroing a vector

      integer m
      double precision x(m)

      integer i
      do i=1,m
         x(i)=0.d0
      enddo
      end

      subroutine add(m,x,y,z)

c	adding 2 vectors			z[] := x[] + y[]

      integer m
      double precision x(m),y(m),z(m)

      integer i
      do i=1,m
         z(i)=x(i)+y(i)
      enddo
      return
      end

      subroutine sub(m,x,y,z)

c	subtracting a vector y from a vector x
      integer m
      double precision x(m),y(m),z(m)

      integer i
      do i=1,m
         z(i)=x(i)-y(i)
      enddo
      return
      end

      subroutine comb(m,x,a,y,z)

c	linear combination of 2 vectors		z[] := x[] + a* y[]

      integer m
      double precision a, x(m),y(m),z(m)

      integer i
      do i=1,m
         z(i)=x(i)+a*y(i)
      enddo
      return
      end


      subroutine multiplyvector(m,a,x,y)

c	multiplying matrix times vector  y := A . x

      integer m
      double precision a(m,m),x(m),y(m)
      integer i,j

      do i=1,m
         y(i) = 0d0
         do j=1,m
            y(i)=y(i)+a(i,j)*x(j)
         enddo
      enddo
      return
      end

      integer function log2(x)

c       the least integer larger than log_2(x)

      double precision x

      log2=0

 8    log2 = log2+1
      if(dble(2**log2).lt.x) goto 8
      return
      end

      integer function nfact(n)

c	factorial function
      integer n

      integer i

      nfact=1
      do  i=1,n
         nfact=nfact*i
      enddo
      return
      end

      double precision function c(n,k)

c	kth coefficient in the nth Pade polynom

      integer n,k
      double precision padenom, padedenom

      integer nfact

      padenom=dble(nfact(2*n-k)*nfact(n))
      padedenom=dble(nfact(2*n)*nfact(k)*nfact(n-k))
      c=padenom/padedenom
      return
      end
