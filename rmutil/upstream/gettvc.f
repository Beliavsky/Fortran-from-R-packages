c
c  rmutil : A Library of Special Functions for Repeated Measurements
c  Copyright (C) 1998, 1999, 2000, 2001 J.K. Lindsey
c
c  This program is free software; you can redistribute it and/or modify
c  it under the terms of the GNU General Public License as published by
c  the Free Software Foundation; either version 2 of the License, or
c  (at your option) any later version.
c
c  This program is distributed in the hope that it will be useful,
c  but WITHOUT ANY WARRANTY; without even the implied warranty of
c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c  GNU General Public License for more details.
c
c  You should have received a copy of the GNU General Public License
c  along with this program; if not, write to the Free Software
c  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
c
c  SYNOPSIS
c
c    subroutine gettvc_f(x,y,xtvc,tvcov,nobs,nind,nknt,ties,
c   +     xu,ndelta,tvcov2,nu,wu,nld,tvcov3)
c
c  DESCRIPTION
c
c    Function to find the most recent value of a time-varying
c covariate not recorded at the same time as the response.
c
c
      subroutine gettvc_f(x,y,xtvc,tvcov,nobs,nind,nknt,ties,
     +     xu,ndelta,tvcov2,nu,wu,nld,tvcov3)
c
c Merge and sort x and xtvc, putting the result, of length nu,
c in xu. Identify members of xtvc with 1 in the binary ndelta,
c and determine the relevant covariate values for those members,
c putting that data in tvcov2.
c                      Daniel F. Heitjan, 27 June 1990
c                      Jim Lindsey, revised 6 October, 1991
c
      implicit none
      integer i,j,n,indx,indu,indk,nm,nm1,nld,nind
      integer nobs(nind),nknt(nind),nu(nind)
      logical ndelta(nind,2*nld),ldone,ties
      double precision x(1),y(1),xtvc(1),wu(1),tvcov(1),
     +     tvcov3(1),xu(nind,2*nld),tvcov2(nind,2*nld),
     +     recx,reck
      nm=0
      nm1=0
      do 1 n=1,nind
         do 9 i=1,nobs(n)
            xu(n,i)=x(nm+i)
 9       continue
         do 19 i=1,nld*2
            ndelta(n,i)=.false.
            tvcov2(n,i)=-1e30
 19      continue
         indx=1
         indk=1
         indu=1
         recx=x(nm+indx)
         reck=xtvc(nm1+indk)
         ldone=.false.
 10      continue
         if(.not.ldone)then
            if(recx.lt.reck)then
               xu(n,indu)=recx
               indx=indx+1
            else
               xu(n,indu)=reck
               ndelta(n,indu)=.true.
               tvcov2(n,indu)=tvcov(nm1+indk)
               if(reck.ne.recx)indk=indk+1
               if(reck.eq.recx)indx=indx+1
            endif
            indu=indu+1
            if(indx.le.nobs(n))then
               recx=x(nm+indx)
            else
               recx=1e30
            endif
            if(indk.le.nknt(n))then
               reck=xtvc(nm1+indk)
            else
               reck=1e30
            endif
            if((recx.ge.1e30).and.(reck.ge.1e30)) ldone=.true.
            go to 10
         endif
         nu(n)=indu-1
         nm=nm+nobs(n)
         nm1=nm1+nknt(n)
 1    continue
      nm=0
      nm1=0
      do 18 n=1,nind
c
c     find time-varying covariate values
c
         if((x(nm+1).ge.xtvc(nm1+1).and.ties).or.
     +        (x(nm+1).gt.xtvc(nm1+1).and..not.ties))then
            reck=tvcov2(n,1)
         else
            reck=0
         endif
         do 49 i=1,nu(n)
            if(i.gt.1.and.ndelta(n,i-1))reck=tvcov2(n,i-1)
            wu(i)=reck
 49      continue
c
c pick the components of wu that correspond to the sampled
c x times, and place them in tvcov
c
         do 2 j=1,nobs(n)
            do 29 i=j,nu(n)
               if(x(nm+j).eq.xu(n,i))then
                  if(x(nm+j).eq.xu(n,i+1).and.ties)then
                     tvcov3(nm+j)=wu(i+1)
                  else
                     tvcov3(nm+j)=wu(i)
                  endif
                  goto 2
               endif
 29         continue
 2       continue
         nm1=nm1+nknt(n)
         nm=nm+nobs(n)
 18   continue
      return
      end
