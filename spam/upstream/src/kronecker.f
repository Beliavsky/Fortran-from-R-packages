      subroutine kroneckermult(xnrow,xent,xcol,xrow,
     &     ynrow,yncol,yent,ycol,yrow,
     &     ent, col, row)

      implicit none
      integer xnrow,ynrow,yncol, xcol(*), xrow(*)
      integer   ycol(*), yrow(*), col(*), row(*)

      double precision  xent(*), yent(*), ent(*)

      integer i,k,j,l,n,nr,xdiffi,ydiffk


      n = 1
      nr = 2
      row(1) = 1
      do i = 1,xnrow
         xdiffi = xrow(i+1)-xrow(i)
         do k = 1,ynrow
            ydiffk = yrow(k+1)-yrow(k)
            do j = 1,xdiffi
               do l = 1,ydiffk
                  ent(n) = xent(j+xrow(i)-1)*yent(l+yrow(k)-1)
                  col(n) = ycol(l+yrow(k)-1)+
     &                    (xcol(j+xrow(i)-1)-1)*yncol

                  n=n+1
               enddo
            enddo
            row(nr) = n
            nr = nr+1
         enddo
      enddo

      return
      end


      subroutine kroneckerf(xnrow,xent,xcol,xrow,
     &     ynrow,yncol,yent,ycol,yrow,
     &     ent1, ent2, col, row)

      implicit none
      integer xnrow,ynrow,yncol, xcol(*), xrow(*)
      integer   ycol(*), yrow(*), col(*), row(*)

      double precision  xent(*), yent(*), ent1(*), ent2(*)

      integer i,k,j,l,n,nr,xdiffi,ydiffk


      n = 1
      nr = 2
      row(1) = 1
      do i = 1,xnrow
         xdiffi = xrow(i+1)-xrow(i)-1
         do k = 1,ynrow
            ydiffk = yrow(k+1)-yrow(k)-1
            do j = 0,xdiffi
               do l = 0,ydiffk
                  ent1(n) = xent(j+xrow(i))
                  ent2(n)= yent(l+yrow(k))
                  col(n) = ycol(l+yrow(k))+
     &                    (xcol(j+xrow(i))-1)*yncol

                  n=n+1
               enddo
            enddo
            row(nr) = n
            nr = nr+1
         enddo
      enddo

      return
      end

