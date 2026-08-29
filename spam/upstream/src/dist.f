C     closestdistXY,  distance between x and x or between x and y
C
C
C     We have four distances implemented:
C       c("euclidean", "maximum", "minkowski", "greatcircle")
C
c     In case we need the distance matrix between x and x, then the
c     following parameters are used as well:
c     if part=-1, lower tri, part=0 the entire matrix
c        part= 1, upper tri only.
c     only values smaller than eta are considered.
c     p power for minkowski



      double precision function euclid(x,y,p)
      implicit none
      double precision x,y,p
      p=p ! to avoid: Warning: Unused dummy argument 'p' at (1)
      euclid=(x-y)**2
      return
      end



      double precision function minkowski(x,y,p)
      implicit none
      double precision x,y,p
      minkowski=abs(x-y)**p
      return
      end



      subroutine closestdist( ncol, x,nrowx, y, nrowy,
     &    part, p, method,
     &    eta, colindices, rowpointers, entries, nnz, iflag)

      implicit none


      double precision euclid, minkowski
      external euclid, minkowski

      integer ncol,nrowx, nrowy, nnz, method,  part,  iflag
      integer colindices(nnz), rowpointers(nrowx+1)
      double precision p, x(nrowx,ncol),y(nrowy,ncol)
      double precision eta, entries(nnz)




      if (method.eq.1) then
         p=2.0
         call closestEdistXY( ncol, x,nrowx, y, nrowy,
     &        part, p, euclid,
     &        eta, colindices, rowpointers, entries, nnz, iflag)
      endif
      if (method.eq.2) then
         p=1.0
         call closestMAXdistXY( ncol, x,nrowx, y, nrowy,
     &        part,
     &        eta, colindices, rowpointers, entries, nnz, iflag)
      endif
      if (method.eq.3) then
         call closestEdistXY( ncol, x,nrowx, y, nrowy,
     &        part, p, minkowski,
     &        eta, colindices, rowpointers, entries, nnz, iflag)
      endif
      if (method.eq.4) then
         call closestGCdistXY( x,nrowx, y, nrowy,
     &        part, p,
     &        eta, colindices, rowpointers, entries, nnz, iflag)
      endif
      return
      end



      subroutine closestEdistXY( ncol, x,xnrow, y, ynrow,
     &    part, p, distfcn,
     &    eta, colindices, rowpointers, entries, nnz, iflag)

      implicit none

      double precision distfcn
      external distfcn

      integer ncol,xnrow, ynrow, nnz,  part,iflag
      integer colindices(nnz),rowpointers(xnrow+1)
      double precision p,x(xnrow,ncol), y(ynrow,ncol)
      double precision eta, entries(nnz)

c     local variables
      integer jja, i,j,k, jfrom, jto
      double precision etap, tmp,pinv


      etap=eta**p
      pinv=1/p

      jja=1

      rowpointers(1)=1
      jfrom = 1
      jto = ynrow


c cycle over all rows of x (independent of part)

      do i= 1,xnrow

         if (part .lt. 0) then
            jto = i
         endif
         if (part .gt. 0) then
            jfrom = i
         endif

         do 10 j = jfrom,jto

c Start calculating the distance (until delta is exceeded)
            tmp = 0.0
            do  k = 1, ncol
               tmp = tmp + distfcn(x(i,k),y(j,k),p)
               if( tmp.gt.etap) goto 10
            enddo
c Delta is not exceeded.

c     in case nnz was too small, recall line to get a better estimate
            if( jja .gt. nnz) then
               iflag = i
               goto 20
            endif

            colindices(jja) = j
            if (abs(p-2) .le. 0.D0) then
               entries(jja) = sqrt(tmp)
            else
               if (abs(p-1) .le. 0.D0) then
                  entries(jja) = tmp
               else
                  entries(jja) = tmp**pinv
               endif
            endif
            jja = jja + 1


 10      continue
         rowpointers(i+1)=jja
      enddo


      if (part.gt.0) then
         rowpointers(xnrow+1)=jja
      endif
      nnz=jja-1
 20   continue

      return
      end

      subroutine closestMAXdistXY( ncol, x,xnrow, y, ynrow,
     &    part,
     &    eta, colindices, rowpointers, entries, nnz, iflag)

      implicit none

      integer ncol,xnrow, ynrow, nnz,  part,iflag
      integer colindices(nnz),rowpointers(xnrow+1)
      double precision x(xnrow,ncol), y(ynrow,ncol)
      double precision eta, entries(nnz)

c     local variables
      integer jja, i,j,k, jfrom, jto
      double precision  tmp



      jja=1

      rowpointers(1)=1
      jfrom = 1
      jto = ynrow


      do i= 1,xnrow

         if (part .lt. 0) then
            jto = i
         endif
         if (part .gt. 0) then
            jfrom = i
         endif

         do 10 j = jfrom,jto

c Start calculating the distance
            tmp = 0.0
            do  k = 1, ncol
               tmp = max(tmp, abs(x(i,k)-y(j,k)))
               if( tmp.gt.eta) goto 10
            enddo

c Delta is not exceeded.
c     (i,j) has a distance smaller than eta.



c     in case nnz was too small, recall line to get a better estimate
            if( jja .gt. nnz) then
               iflag = i
               goto 20
            endif

            colindices(jja) = j
            entries(jja) = tmp
            jja = jja + 1


 10      continue
         rowpointers(i+1)=jja
      enddo


      if (part.gt.0) then
         rowpointers(xnrow+1)=jja
      endif
      nnz=jja-1
 20   continue

      return
      end



      subroutine closestGCdistXY( x,nx, y, ny,
     &    part,p,
     &    eta, colindices, rowpointers, entries, nnz, iflag)

      implicit none

      integer nx, ny, nnz, colindices(nnz),rowpointers(nx+1)
      integer part, iflag
      double precision x(nx,2), y(ny,2), p, eta, entries(nnz)

c     local variables
      logical equi
      integer jja, i,j, jfrom, jto
      double precision etap, tmp, tmp1, tmp2
      double precision rad, thres
      double precision scy12(ny), ccy12(ny), sy2(ny)
      double precision scx12,     ccx12,     sx2

      parameter (rad = 0.01745329251994329)
      parameter (thres = 0.99999999999)


c     Great savings if we know that x=y.
c     Changes for archaic d... to ... trigonometric fcn
      if (p .lt. 0) then
         equi=.TRUE.
         p=-p
      else
         equi= .FALSE.
      endif

      jja=1

      etap=cos(eta*rad)
      rowpointers(1)=1
      jfrom = 1
      jto = ny


      DO j=1,ny
         tmp1=y(j,1)*rad
         tmp2=y(j,2)*rad
         ccy12(j)=cos(tmp1)*cos(tmp2)
         scy12(j)=sin(tmp1)*cos(tmp2)
         sy2(j)=sin(tmp2)
      ENDDO



      do i= 1,nx

c     x2 is missing if equi=.TRUE. and we reuse the y stuff
         if (equi .eqv. .TRUE.) then
            ccx12=ccy12(i)
            scx12=scy12(i)
            sx2=sy2(i)
         else
            tmp1=x(i,1)*rad
            tmp2=x(i,2)*rad
            ccx12=cos(tmp1)*cos(tmp2)
            scx12=sin(tmp1)*cos(tmp2)
            sx2=sin(tmp2)
         endif

         if (part .lt. 0) then
            jto = i
         endif
         if (part .gt. 0) then
            jfrom = i
         endif

         do 10 j = jfrom,jto


c     Start calculating the distance
            tmp = ccx12 * ccy12(j) + scx12 * scy12(j)  + sx2*sy2(j)


            if (tmp .lt. etap) goto 10
c     Delta is not exceeded.

c     Due to numerical instabilities, we need the following... 0.15-2:
c     Patch suggested at code clinics.
            if  (tmp .ge. thres) then
               tmp = 0.0
            else
               tmp = acos( tmp)
            endif

c     (i,j) has a distance smaller than eta.

c     In case nnz was too small, recall line to get a better estimate
            if( jja .gt. nnz) then
               iflag = i
               goto 20
            endif

            colindices(jja) = j
            entries(jja) = tmp*p
            jja = jja + 1


 10      continue
         rowpointers(i+1)=jja
      enddo


      if (part.gt.0) then
         rowpointers(nx+1)=jja
      endif
      nnz=jja-1
 20   continue

      return
      end




