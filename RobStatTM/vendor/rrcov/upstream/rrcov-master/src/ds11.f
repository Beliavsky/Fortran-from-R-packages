      subroutine rlds(n,np,nresamp,x,tune,wk,locat,cov,maxres,
     1     nresper,w,z,icent,iwork)
c
c
c  INPUT:  'n' = number  observations (integer);
c         'np' = number of indepent  variables (integer);
c    'nresamp' = mumber of  resamples required (integer), may not be reached
c                if too many of the subsamples of size 'np', chosen out of
c                the observed vectors, are in a hyperplane;
c                if nresamp=0, ALL subsamples are taken;
c          'x' = real (n*np) matrix of observed values;
c      'tune'  = tuning constant used in calculation of weights;
c                by default should be equal to the square root of the
c                95-percentile of the chi-square distribution with 'np'
c                 degrees of freedom (real).
c       'seed' = seed for calculating random numbers (real);
c         'wk' = real work vector of length (4*n+np).  On output,
c                wk(1),..,wk(n) contain the weights assigned to
c                each observation;
c     'maxres' = maximum nunber of  resamples to be performed (integer),
c                including those which are discarded due to linearly
c                dependent subsamples.
c        'icent'= if 0 the observations are centered with 0
c
c  OUTPUT:  'locat' = real vector of location parameters, of length
c                     'np';
c             'cov' =  real (np*np) Donoho-Stahel covariance matrix;
c         'nresper' = number of  valid resamples performed (integer).
c               'w' = weights
c               'z' = outlyigness
c
c  NOTE:  all real variables are double precision.
c
      implicit double precision (a-h,o-z)
      double precision locat(np)
      dimension x(n,np),wk(4*n+np),cov(np,np),w(n),z(n), iwork(4*n+np)
      nr=1
      nind=nr+n
      naux=nind+np

C      CALL INTPR('ENTER rlds. nresamp=',-1,nresamp,1)
C      CALL INTPR('maxres=',-1,maxres,1)
C      CALL INTPR('icent=',-1,icent,1)

      call rndstart()
      call rlweights(n,np,nresamp,x,tune,w,z,locat,wk(nr),iwork(nind),
     1             cov,wk(naux),maxres,nresper,icent)

      call rldonostah(n,np,x,w,locat,cov,icent)

      call rndend()
      return
      end


      subroutine rlweights(n,np,nresamp,x,c,w,z,a,b,ind,wk,
     1                   u,maxres,nresper,icent)
      implicit double precision (a-h,o-z)
      dimension x(n,np),z(n),a(np),b(n),w(n),ind(np),u(n)
      dimension wk(np,np)

C      CALL INTPR('ENTER rlweights',-1,0,0)

      k1=(np-1)+(n+1)/2
      k2=(np-1)+(n+2)/2
      z1=dble(k1)
      zn=dble(n)
      z3=(1+(z1/zn))/2

        call rlquntbi(z3, cc)

        do i=1,n
           z(i)=-1.
        enddo

      nresper=0
      if(np.eq.1) then
         call rlprocess(n,np,nresper,x,a,b,w,z,ind,wk,u,k1,
     +        k2,cc,icent)
      elseif (nresamp.eq.0) then
         call rlall(n,np,nresper,x,a,b,w,z,ind,wk,u,k1,k2,cc,icent)
      else
         k=0
         do while (k.lt.maxres.and.nresper.lt.nresamp)
            k=k+1
            call rlsubsamp(n,np,ind)
            call rlprocess(n,np,nresper,x,a,b,w,z,ind,wk,u,k1,
     +           k2,cc,icent)
          enddo
      endif

C      CALL DBLEPR('zi',-1,z,n)

      do i=1,n
         call rlrwetml(z(i)/c,w(i))
      enddo

C      CALL DBLEPR('EXIT rlweights: wi=',-1,w,n)
      return
      end

      subroutine rlall(n,np,nresper,x,a,b,w,z,ind,wk,u,k1,k2,cc,icent)
      implicit double precision (a-h,o-z)
      dimension x(n,np),z(n),a(np),b(n),w(n),ind(np),u(n)
      dimension wk(np,np)
      do j=1,np
         ind(j)=j
         enddo
      call rlprocess(n,np,nresper,x,a,b,w,z,ind,wk,u,k1,k2,cc,icent)
      j=0
      do while (np-j.ge.1)
         if (ind(np-j).eq.n-j) then
            j=j+1
            else
            ind(np-j)=ind(np-j)+1
            do k=np-j+1,np
               ind(k)=ind(k-1)+1
            enddo
            call rlprocess(n,np,nresper,x,a,b,w,z,ind,
     +           wk,u,k1,k2,cc,icent)
            j=0
         endif
      enddo
      return
      end


      subroutine rlprocess(n,np,nresper,x,a,b,w,z,ind,wk,
     +     u,k1,k2,cc,icent)
      implicit double precision (a-h,o-z)
      dimension x(n,np),z(n),a(np),b(n),w(n),ind(np),u(n)
      dimension wk(np,np)
      data tola,tolr,big1,big2 /1.d-15, 1.d-8,1.d+2,1.d+15/

C      CALL INTPR('ENTER rlprocess',-1,0,0)

      CALL RCHKUSR()

      ierr=0
      if(np.gt.1) then
         call rlvectora(n,np,x,a,ind,wk,icent,ierr)
      endif

C      CALL INTPR('IERR',-1,ierr,1)
C      CALL DBLEPR('A',-1,a,np)

      if (ierr.eq.0) then
         nresper=nresper+1

C       VT::19.07.2010
C       Handle the univariate case
         if(np.eq.1) then
            do i=1,n
               b(i)=x(i,1)
            enddo
         else
            do i=1,n
               b(i)=0.
               do j=1,np
                  b(i)=b(i)+x(i,j)*a(j)
               enddo
            enddo
         endif

         bmed=0.0d0
         if(icent.ne.0) bmed=rlamed(b,n,u)
         do i=1,n
            w(i)=abs(b(i)-bmed)
         enddo
         ww=0
         do i=1,n
            ww=ww+w(i)
         enddo
         ww=ww/n
         if(ww.ge.tola) then
            call rlsort(w,n,1)
            bmad=(w(k1)+w(k2))/2
            bmad=bmad/cc
            if(bmad.ge.tolr *ww) then
               do i=1,n
                   aux=abs(b(i)-bmed)/bmad
                   if (aux.gt.z(i)) z(i)=aux
               enddo
            else
               do i=1,n
                  if(abs(b(i)-bmed).gt. big1*bmad) z(i)=big2
               enddo
            endif
         endif
      endif
      return
      end

      subroutine rlvectora(n,np,x,a,ind,wk,icent,ierr)
      implicit double precision (a-h,o-z)
      dimension x(n,np),a(np),ind(np),wk(np,np)
      do k=1,np
         do j=1,np
            wk(j,k)=x(ind(k),j)
         enddo
      enddo
      call rldirec(wk,np,np,icent,ierr,a)
      return
      end


      subroutine rldonostah(n,np,x,w,locat,cov,icent)
      implicit double precision (a-h,o-z)
      double precision locat(np)
      dimension x(n,np),w(n),cov(np,np)
      sumw=0.
      sumw2=0.
      do i=1,n
         sumw=sumw+w(i)
         sumw2=sumw2+w(i)*w(i)
      enddo
      do j=1,np
         locat(j)=.0
      enddo
      if(icent.eq.1)then
         do j=1,np
            locat(j)=0.
            do i=1,n
               locat(j)=locat(j)+w(i)*x(i,j)
            enddo
            locat(j)=locat(j)/sumw
         enddo
      endif
      do j=1,np
         do k=1,np
            cov(j,k)=0.
            do i=1,n
               cov(j,k)=cov(j,k)+w(i)*(x(i,j)-locat(j))*
     1         w(i)*(x(i,k)-locat(k))
            enddo
            cov(j,k)=cov(j,k)/sumw2
         enddo
      enddo
      return
      end


      subroutine rlsubsamp(n,np,ind)
      implicit double precision (a-h,o-z)
      dimension ind(np)
      en=dble(n)
C      call roblibrunif(RND)
      RND = unifrnd()
      ind(1)=int(en*RND+1.)
      if (np.eq.1) return
      k=2
c   10 call roblibrunif(RND)
10    RND = unifrnd()
      ind(k)=int(en*RND+1.)
      do i=1,k-1
         if (ind(i).eq.ind(k)) go to 10
      enddo
      if (k.eq.np) return
      k=k+1
      go to 10
      end


      double precision function rlamed(z,n,aux)
      implicit double precision (a-h,o-z)
      DIMENSION Z(n),aux(n)
      DO 100 I=1,N
      AUX(I)=Z(I)
  100 CONTINUE
      CALL rlSORT (AUX,N,1)
      I=N/2
      K=I*2
      rlamed=AUX(I+1)
      IF (k.GE.N) rlamed=(rlamed+AUX(I))/2.
      RETURN
      END

      SUBROUTINE RLSORT (A,N,SWITCH)
      implicit double precision (a-h,o-z)
      DIMENSION A(n)
      INTEGER SWITCH
      IF (N.LE.1) GO TO 999
      M=1
106   M=M+M
      IF(M.LE.N) GO TO 106
      M=M-1
994    M=M/2
      IF (M.EQ.0) GO TO 999
      KK=N-M
      J=1
992   I=J
996   IM=I+M

CCCC    VT::18.07.2022
CCCC      IF(SWITCH) 810,810,800
      IF(SWITCH.LE.0) GO TO 810
800    IF (A(I).GT.A(IM)) GO TO 110
      GO TO 995
810    IF(A(I).LT.A(IM)) GO TO 110
995   J=J+1
      IF(J.GT.KK) GO TO 994
      GO TO 992
110   TEMP=A(I)
      A(I)=A(IM)
      A(IM)=TEMP
       I=I-M
      IF (I.LT.1) GO TO 995
      GO TO 996
999    RETURN
      END


      double precision function rldprodd(x,y,nn)
      implicit double precision (a-h,o-z)
      dimension x(nn), y(nn)
      rldprodd=0.
      do  i=1,nn
         rldprodd=rldprodd+x(i)*y(i)
      enddo
      return
      end


        double precision function rlrobustdnorm(x,nn)
        implicit double precision (a-h,o-z)
        dimension x(nn)

        rlrobustdnorm=rldprodd(x,x,nn)
        rlrobustdnorm=dsqrt(rlrobustdnorm)
        return
        end


        subroutine rlxnorma(x,nn,ierr,tol)
        implicit double precision (a-h,o-z)
        dimension x(nn)
        ierr=1
        dn=rlrobustdnorm(x,nn)
        if (dn.le.tol) then
           ierr=1
           return
        else
           ierr=0
        endif
        do  i=1,nn
           x(i)=x(i)/dn
        enddo
        return
        end

        subroutine rlorthog(xx,nn,mm,nmain,ierr)
        implicit double precision (a-h,o-z)
        dimension xx(nmain,mm)
        data tola,tolr /1.d-15, 1.d-8/
C In original code tolb was never initialized (was 0 on Solaris, random on HP)
        tolb = tola
        do  j=1,mm
           call  rlxnorma(xx(1,j),nn,ierr,tola)
           if (ierr.gt.0) return
        enddo
        mm1=mm-1
        do  j=1,mm1
           call  rlxnorma(xx(1,j),nn,ierr,tolr)
           if (ierr.ne.0) return
           j1=j+1
           do k=j1,mm
              dp=rldprodd(xx(1,j),xx(1,k),nn)
              do  i=1,nn
                 xx(i,k)=xx(i,k)-xx(i,j)*dp
              enddo
           enddo
        enddo
        call  rlxnorma(xx(1,mm),nn,ierr,tolb)
C       if (ierr .ne. 0) write(*,*) 'rlxnorma(...,tolb) failed!'
        return
        end


        subroutine rlortdir(xx,mm,nmain,dire)
        implicit double precision (a-h,o-z)
        dimension xx(nmain,mm), dire(mm)
        tol=1./dsqrt(dble(mm))
        mm1=mm-1
        do  k=1,mm
           do  i=1,mm
              dire(i)=0
              do  j=1,mm1
                 dire(i)=dire(i)-xx(i,j)*xx(k,j)
              enddo
           enddo
           dire(k)=dire(k)+1
           dn=rlrobustdnorm(dire,mm)
           if (dn.ge.tol) goto 40
        enddo
40      do  i=1,mm
           dire(i)=dire(i)/dn
        enddo
        return
        end


        subroutine rldirec(xx,mm,nmain,icent,ierr,dire)
        implicit double precision (a-h,o-z)
        dimension xx(nmain,mm), dire(mm)
        mm1=mm
        if (icent.ne.0)then
           mm1=mm-1
           do  k=1,mm1
              do  i=1,mm
                 xx(i,k)=xx(i,k)-xx(i,mm)
              enddo
           enddo
        endif
        call rlorthog(xx,mm,mm1,nmain,ierr)
        if (ierr.eq.0) call rlortdir(xx,mm,nmain,dire)
        return
        end

      SUBROUTINE RLRWETML(X,P)
C.......................................................................
      DOUBLE PRECISION X,AX,P,COEF,ZERO
      DIMENSION COEF(4)
      DATA ZERO/0.D0/
      DATA COEF/-19.7187928669416D0,
     +           82.3045267489739D0,
     +         -105.4526748971229D0,
     +           42.8669410150906D0/
C-----------------------------------------------------------------------
C     CALCULATE WEIGHT FUNCTION FOR REWEIGHTING
C     NOTE: X >= 0
C-----------------------------------------------------------------------
      AX = DABS(X)
      IF (AX .GE. 1.D0) THEN
         P = ZERO
      ELSE IF (AX .LE. 0.8D0) THEN
         P = 1.D0
      ELSE
         P = COEF(1)+COEF(2)*AX**2+COEF(3)*AX**4+COEF(4)*AX**6
      ENDIF

C      CALL DBLEPR('IN RLRWETML',-1,p,1)

      RETURN
      END
C=======================================================================
      SUBROUTINE RLQUNTBI(P,X)
C.......................................................................
      DOUBLE PRECISION C(6),P,P1,T,X,XN,XZ
      DATA C(1),C(2),C(3),C(4),C(5),C(6)/
     +     2.515517D0,0.802853D0,0.010328D0,
     +     1.432788D0,0.189269D0,0.001308D0/
C-----------------------------------------------------------------------
C     INVERSE OF GAUSSIAN DISTRIBUTION FUNCTION
C-----------------------------------------------------------------------
C     P: I, PROBABILITY,
C     X: O, QUANTILE.
C-----------------------------------------------------------------------
      P1=P
      IF (P .GT. 0.5D0) P1=1.D0-P
      T=DSQRT(-2.D0*DLOG(P1))
      XZ=(C(3)*T+C(2))*T+C(1)
      XN=((C(6)*T+C(5))*T+C(4))*T+1.D0
      X=T-XZ/XN
      IF (P .LT. 0.5D0) X=-X
      RETURN
      END
