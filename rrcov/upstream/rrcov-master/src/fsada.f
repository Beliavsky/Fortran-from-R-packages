      SUBROUTINE FSADA(X,NE,NV,NG,IGRP,XM,XC,XND,
     1          IHALF,NIT,IDSAV,MRAND,IERR,DETC,ITRACE)

        IMPLICIT DOUBLEPRECISION (A-H,O-Z)

C...Default number of iterations to be performed(if NIT=0)
        PARAMETER (ITK0=500)


        DIMENSION X(NV,NE), IGRP(NE)
        DIMENSION XM(NV,NG), XC(NV,NV), XND(NG)

        DIMENSION ave(NV,NG),cov(NV,NV),
     1 datimg(NV,NE),ustu(NE),ibasis(NE),
     3 isin(NE),wtv(NG),cinv(NV,NV),
     4 IDSAV(NE)

C...Set the default number of iterations if not supplied
        IF(NIT.EQ.0) NIT=ITK0

      if(ITRACE .ge. 2) then
        CALL INTPR('Entering FSADA - NIT: ',-1,NIT,1)
      endif

C...Set the coverage for max breakdown point, if not supplied
      IF(IHALF.LE.0) IHALF = (NE+NV+1)/2

      call reader(X,NE,NV,NG,ave,cov,cinv,datimg,wtv,ustu,DETC,
     2 IGRP,ibasis,isin,IHALF,
     3 XC,XM,XND,IDSAV,NIT,MRAND,ITRACE)

C...the cross-product matrix is returned - scale it
        DO I=1,NV
           DO J=1,NV
              XC(I,J) = XC(I,J)/(IHALF-NG)
           ENDDO
        ENDDO


        END


      subroutine reader(data,ncas,nord,npop,
     1ave,cov,cinv,datimg,wtv,ustu,deter,
     2igp,ibasis,isin,ncover,
     3bescov,besave,beswtv,ibsbes,nsim,iseed,ITRACE)

      implicit doubleprecision (a-h,o-z)
      dimension data(nord,ncas),ave(nord,npop),cov(nord,nord),
     1 datimg(nord,ncas),
     2 ustu(ncas),igp(ncas),ibasis(ncas),
     3 isin(ncas),wtv(npop),cinv(nord,nord),
     4 bescov(nord,nord),besave(nord,npop),ibsbes(ncas),
     5 beswtv(npop),
     6 index(5000),devi(100),
     7 detrtb(100),icount(100),isgntb(100)

      data xix, xia, xic /.31415925d0,17119.d0,0.1221d0/,
     1 one /1.d0/

      xix = xix * iseed

CCC     WRITE(*,*) ncas,' cases'
c       ncover = cover * ncas + 0.5
      if(ITRACE .ge. 2) then
        CALL INTPR('Entering READER - ncas: ',-1,ncas,1)
        CALL INTPR('Entering READER - ncover: ',-1,ncover,1)
      endif

        do j = 1, ncas
           index(j) = j
        enddo

        do j = 1, npop
           wtv(j) = 0
           do i = 1, nord
              ave(i,j) = 0
           enddo
        enddo

      do j = 1, nord
        do i = 1, nord
           cov(i,j) = 0
        enddo
      enddo

      do icas = 1, ncas
        ngp = igp(icas)
        wold = wtv(ngp)
        wnew = wold + one
        wtv(ngp) = wnew
        rati = wold / wnew
        do i = 1, nord
           devi(i) =  data(i,icas) - ave(i,ngp)
           ave(i,ngp) = ave(i,ngp) + devi(i) / wnew
           do j = 1, i
              cov(i,j) = cov(i,j) + devi(i) * devi(j) * rati
           enddo
        enddo
      enddo

      edf = ncas - npop
      do i = 1, nord
        do j = 1, i
           cov(i,j) = cov(i,j) / edf
           cov(j,i) = cov(i,j)
        enddo
      enddo

CCC      WRITE(*,'('' Full data set''/'' Grp  size    averages'')')
CCC      WRITE(8,'('' Full data set''/'' Grp  size    averages'')')
CCC      do j = 1, npop
CCC		WRITE(*,101) j,int(wtv(j)),(ave(i,j),i=1,nord)
CCC		WRITE(8,101) j,int(wtv(j)),(ave(i,j),i=1,nord)
CCC 101		format(/i4,i5,(t12,6g11.4))
CCC	enddo

CCC      WRITE(*,*) ' Covariance matrix'
CCC      WRITE(8,'(/'' Covariance matrix'')')
CCC      do j = 1, nord
CCC		WRITE(*,102) j,(cov(i,j),i=1,nord)
CCC		WRITE(8,102) j,(cov(i,j),i=1,nord)
CCC 102		format(/i5,(t12,6g11.4))
CCC      enddo

      deter = 1
      ixlo = 1
      do i = 1, nord
CCC     WRITE(*,'('' Sweeping pivot, deter'',2g15.7)') cov(i,i),deter
        call zsweep(cov,nord,i,deter)
      enddo

CCC      WRITE(*,'(/'' Log10 determinant of overall covariance matrix'',
CCC     1 f10.4)') log10(deter)
CCC      WRITE(8,'(/'' Log10 determinant of overall covariance matrix'',
CCC     1 f10.4)') log10(deter)

      if(ITRACE .ge. 2) then
        xdet = log10(deter)
        CALL DBLEPR('Initialization ready - log det: ',-1,xdet,1)
      endif

      verbes = 1.d30
      if (ncover .ge. ncas) return


      cover = 1.0*ncover/ncas

CCC      WRITE(*,'(//'' Start search for high breakdown estimators''
CCC     1 '' of coverage'',f7.3)') cover
CCC      WRITE(8,'(//'' Start search for high breakdown estimators''
CCC     1 '' of coverage'',f7.3)') cover

      edf = ncover - npop
      corter = nord * log10(edf)
      nsol = 0
      do 80 loopo = 1, nsim
        do i = ncas, 1, -1
           fi = i
           xix = mod(xix * xia + xic, one)
           inx = xix * fi + one
           if(inx .ne. i) then
              itemp = index (i)
              index (i) = index (inx)
              index (inx) = itemp
           endif
           ibasis(i) = index(i)
        enddo

      if(ITRACE .ge. 2) then
        CALL INTPR('Entering iteration: ',-1,loopo,1)
      endif

        call itera(data,ave,cov,cinv,datimg,wtv,ustu,deter,
     1 igp,ibasis,isin,nord,ncas,npop,ncover)
        isgn = isigna(ibasis,ncover)

        do i = 1, nsol
           if(isgn.eq.isgntb(i).and.abs(deter/detrtb(i)-one)
     1    .lt. .001) then
              icount(i) = icount(i) + 1
              go to 135
           endif
        enddo

      do i = 1, ncover
        do j = 1, i
           if (ibasis(j) .gt. ibasis(i)) then
              itemp = ibasis(i)
              ibasis(i) = ibasis(j)
              ibasis(j) = itemp
           endif
        enddo
      enddo

CCC   WRITE(*,'(/''Loop'',i6,'' New feasible solution.  Retained '',
CCC     1 ''cases''/(20i4))') loopo,(ibasis(i),i=1,ncover)
CCC      detlog = log10(deter) - corter
CCC      WRITE(*,'('' log10 determinant'',
CCC     1 g15.7/'' Grp  size    averages'')') detlog

CCC      do 145 j = 1, npop
CCC 145  WRITE(*,'(i4,i5,(t12,6g11.4))') j,int(wtv(j)),(ave(i,j),i=1,nord)

CCC     WRITE(*,*) ' Covariance matrix'
CCC     do j = 1, nord
CCC             WRITE(*,'(t12,6g11.4)') (cov(i,j)/edf,i=1,nord)
CCC     enddo

        nsol = nsol + 1
        isgntb(nsol) = isgn
        detrtb(nsol) = deter
        icount(nsol) = 1
  135  continue

        if (deter .lt. 0.999999d0 * verbes) then
CCC        WRITE(*,'(//'' New optimum'')')
CCC        WRITE(8,'(/''Loop'',i6,'' New feasible solution.  Retained '',
CCC     1   ''cases''/(20i4))') loopo,(ibasis(i),i=1,ncover)
CCC        WRITE(8,'('' log10 determinant'',g15.7)') detlog
        verbes = deter
        do i = 1, nord
                do j = 1, npop
                beswtv(j) = wtv(j)
                besave(i,j) = ave(i,j)
                enddo
                do j = 1, nord
                bescov(i,j) = cov(i,j)
                enddo
        enddo
        do i = 1, ncover
                ibsbes(i) = ibasis(i)
        enddo
        endif
 80   continue

CCC      WRITE(*,'('' Criterion and count of different feasible '',
CCC     1 ''solutions'')')
CCC      WRITE(*,'(5(f7.3,i5))') (log10(detrtb(i))-corter,
CCC     1  icount(i), i=1,nsol)
CCC      WRITE(8,'('' Criterion and count of different feasible '',
CCC     1 ''solutions'')')
CCC      WRITE(8,'(5(f7.3,i5))') (log10(detrtb(i))-corter,
CCC     1  icount(i), i=1,nsol)
CCC      write(*,103) (ibsbes(i),i=1,ncover)
CCC      write(8,103) (ibsbes(i),i=1,ncover)
CCC 103  format(/' Best feasible solution.  Cases covered are'/(15i5))
CCC      WRITE(*,'(/'' Grp  size    averages'')')
CCC      WRITE(8,'(/'' Grp  size    averages'')')

CCC      do j = 1, npop
CCC        WRITE(*,101) j,int(beswtv(j)),(besave(i,j),i=1,nord)
CCC        WRITE(8,101) j,int(beswtv(j)),(besave(i,j),i=1,nord)
CCC      enddo

CCC      WRITE(*,*) ' Covariance matrix'
CCC      WRITE(8,'(/'' Covariance matrix'')')

CCC        do j = 1, nord
CCC        WRITE(*,102) j,(bescov(i,j)/edf,i=1,nord)
CCC        WRITE(8,102) j,(bescov(i,j)/edf,i=1,nord)
CCC        enddo

      end

      subroutine itera(data,ave,cov,cinv,datimg,wtv,ustu,deter,
     1 igp,ibasis,isin,nord,ncas,npop,ncover)

      implicit double precision (a-h,o-z)
      dimension data(nord,ncas),ave(nord,npop),cov(nord,nord),
     1 cinv(nord,nord),datimg(nord,ncas),ustu(ncas),
     2 ibasis(ncover),isin(ncas),igp(ncas),devi(200),devj(200),
     3 upfac(100), dnfac(100),wtv(npop)
      data one /1.d0/, big /1.d10/

C       Initialize to avoid warnings
        iout = 0
        icasot = 0
        jin = 0

        do j = 1, npop
           wtv(j) = 0
           do i = 1, nord
              ave(i,j) = 0
           enddo
        enddo
        do j = 1, nord
           do i = 1, nord
              cov(i,j) = 0
           enddo
        enddo
        do i = 1, ncas
           isin(i) = 0
        enddo

c     Get initial covariance matrix
      do inx = 1, ncover
        icas = ibasis(inx)
        isin(icas) = 1
        ngp = igp(icas)
        wold = wtv(ngp)
        wnew = wold + one
        wtv(ngp) = wnew
        rati = wold / wnew
        do i = 1, nord
           devi(i) =  data(i,icas) - ave(i,ngp)
           ave(i,ngp) = ave(i,ngp) + devi(i) / wnew
           do j = 1, i
              cov(i,j) = cov(i,j) + devi(i) * devi(j) * rati
           enddo
        enddo
      enddo

      do i = 1, nord
        do j = 1, i
           cov(j,i) = cov(i,j)
           cinv(i,j) = cov(i,j)
           cinv(j,i) = cov(i,j)
        enddo
      enddo

      deter = 1
      ixlo = 1
      do i = 1, nord
        call zsweep(cinv,nord,i,deter)
      enddo
      do i = 1, nord
        do j = 1, nord
           cinv(i,j) = -cinv(i,j)
        enddo
      enddo
      call verify1(cov,cinv,nord)

c     Major loop point
 30   continue
      do  i = 1, npop
        wtvi = wtv(i)
        upfac(i) = sqrt(wtvi / (wtvi + one))
        dnfac(i) = big
        if (wtvi .gt. one) dnfac(i) = sqrt(wtvi / (wtvi - one))
      enddo

c     Get images of cases
      do 41 j = 1, ncas
        ngp = igp(j)
        ustu(j) = 0
        do 40 i = 1, nord
           sum = 0
           do 45 k = 1, nord
              sum = sum + cinv(i,k) * (data(k,j) - ave(k,ngp))
45         continue
              datimg(i,j) = sum
              ustu(j) = ustu(j) + sum * (data(i,j) - ave(i,ngp))
40      continue
41    continue

c     Consider possible case swaps
        best = one
        do 50 i = 1, ncover
           icas = ibasis(i)
           ngp = igp(icas)
           if(wtv(ngp) .eq. one) go to 50
c
c     dont remove the only case in a group
c
           firfac = one - dnfac(ngp) ** 2 * ustu(icas)
           if (firfac .gt. best) go to 50
           do 55 j = 1, ncas
              if (isin(j) .eq. 1) go to 55
c
c             do pretest
c
              jgp = igp(j)
              if(jgp .ne. ngp) then
c
c       (we need special codes when the two are in the same group)
c
                 factor=firfac*(one+upfac(jgp)**2*ustu(j))
                 if (factor .ge. best) go to 55
c       (cannot beat what we have already)
c
                 sum = 0
                 do 60 k = 1, nord
                    sum=sum+(data(k,icas)-ave(k,ngp))*datimg(k,j)
60               continue
                    factor=factor+(sum*upfac(jgp)*dnfac(ngp))**2
                    if(factor.lt.0) then
CCC          WRITE(*,*) 'Impossible factor. dnfac,ustu(icas),firfac',
CCC     1       dnfac(ngp),ustu(icas),firfac,' upfac,ustu(j)',
CCC     2       upfac(jgp),ustu(j),' sum', sum
CCC          WRITE(*,*) ' wtv(ngp)', wtv(ngp)
                       do 155 ik = 1, ncover
                          ikk = ibasis(ik)
CCC                          if(igp(ikk).eq.ngp) WRITE(*,*) 'hit',ik,ikk
155                    continue
                    endif
                    if(factor .lt. best) then
                        best = factor
                        iout = i
                        icasot = icas
                        jin = j
                    endif
                 else
                    sum1 = 0
                    sum2 = 0
                    divis = wtv(ngp) - one
                    do 61 k = 1, nord
                       ui = data(k,icas) - ave(k,ngp)
                       vi = data(k,j) - ave(k,ngp) + ui / divis
                       vimg = datimg(k,j) + datimg(k,icas) / divis
                       sum1 = sum1 + vi * vimg
                       sum2 = sum2 + ui * vimg
61                  continue
                    factor = firfac * (one + divis * sum1 / wtv(ngp)) +
     1                  sum2 ** 2
                    if (factor .lt. best) then
                    best = factor
                    iout = i
                    icasot = icas
                    jin = j
                endif
             endif
55         continue
50      continue

        if(best .eq. one) go to 90
c
c     There is a swap that improves things.  Make it
c
           isin(icasot) = 0
           ibasis(iout) = jin
           isin(jin) = 1
           deter = deter * best
c
c     Do the downdate, removing case icasot
c
        ngp = igp(icasot)
        jgp = igp(jin)
        wold = wtv(ngp)
        wnew = wold - one
        rati = wold / wnew
        wtv(ngp) = wnew
        fact = rati / (one - rati * ustu(icasot))
        do 70 i = 1, nord
            devi(i) = data(i,icasot) - ave(i,ngp)
            ave(i,ngp) = ave(i,ngp) - devi(i) / wnew
            devj(i) = data(i,jin) - ave(i,jgp)
            do 71 j = 1, i
                cov(i,j) = cov(i,j) - devi(i) * devi(j) * rati
                cov(j,i) = cov(i,j)
                cinv(i,j) = cinv(i,j) + datimg(i,icasot) *
     1              datimg(j,icasot) * fact
                cinv(j,i) = cinv(i,j)
71          continue
70      continue
        call verify1(cov,cinv,nord)
c
c     Now do the update, adding case jin
c
        wold = wtv(jgp)
        wnew = wold + one
        wtv(jgp) = wnew
        rati = wold / wnew
        sum2 = 0
        do 80 i = 1, nord
            ave(i,jgp) = ave(i,jgp) + devj(i) / wnew
            sum = 0
            do 81 j = 1, nord
                cov(i,j) = cov(i,j) + rati * devj(i) * devj(j)
                sum = sum + cinv(i,j) * devj(j)
81          continue
            devi(i) = sum
            sum2 = sum2 + devi(i) * devj(i)
80      continue
        factor = rati / (one + rati * sum2)
        do 85 i = 1, nord
            do 86 j = 1, nord
                cinv(i,j) = cinv(i,j) - devi(i) * devi(j) * factor
86          continue
85      continue
        call verify1(cov,cinv,nord)
      go to 30

 90   return
      end

      subroutine verify1(cov,cinv,nord)
      implicit double precision (a-h,o-z)
      dimension cov(nord,nord),cinv(nord,nord)
      data one /1.d0/
      biger = 0
      do 5 i = 1, nord
      do 10 j = 1, nord
      sum = 0
      do 15 k = 1, nord
      sum = sum + cov(i,k) * cinv(k,j)
 15   continue
      if (i .eq. j) then
        biger = max(biger,abs(sum-one))
        else
        biger = max(biger,abs(sum))
        endif
 10   continue
 5    continue
      if (biger .gt. 0.001) then
CCC        WRITE(*,*) 'Inversion error, departure from I is',biger
        return
        endif
      return
      end

      SUBROUTINE ZSWEEP (COV,NORD,NEL,DETER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COV(NORD,NORD)
      DATA ONE/1.D0/
      TEMP=COV(NEL,NEL)
      DETER=DETER*TEMP
      IF (NORD.GT.1) GO TO 10
      COV(1,1)=ONE/TEMP
      RETURN
10    DO 30 I=1,NORD
        IF (I.EQ.NEL) GO TO 30
        DO 20 J=1,I
          IF (J.EQ.NEL) GO TO 20
          COV(J,I)=COV(I,J)-COV(I,NEL)*COV(NEL,J)/TEMP
          COV(I,J)=COV(J,I)
20      CONTINUE
30    CONTINUE
      COV(NEL,NEL)=ONE
      DO 40 I=1,NORD
        COV(NEL,I)=-COV(I,NEL)/TEMP
        COV(I,NEL)=COV(NEL,I)
40    CONTINUE
      RETURN
      END

      FUNCTION ISIGNA (LIST,NROW)
      DIMENSION LIST(NROW)
      DATA NPRIM1/30931/, NPRIM2/59473/
      ISIG1 = 43
      ISIG2 = 23
      DO 10 I=1,NROW
        ITEMP=LIST(I)+1000
        ISIG1 = MOD(ISIG1*ITEMP,NPRIM1)
        ISIG2 = MOD(ISIG2*ITEMP,NPRIM2)
10    CONTINUE
      ISIGNA = ISIG1 * ISIG2
      RETURN
      END
