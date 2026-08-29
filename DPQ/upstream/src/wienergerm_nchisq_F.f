      double precision function h(y)
      double precision y, one,half, dlog
      data one,half/1.d0, 0.5d0/

      h=((one-y)* dlog(one-y) + y- half*y*y)/y/y
      return
      end

      subroutine noncechi(variant,argument,noncentr,df,p,ifault)

      integer variant ! 1="f"(irst)  2="s"(econd)  [Fortan char N.A. in call from R]
      double precision argument,noncentr,df,p
      integer ifault

      double precision mu2,s,y,help1,help2,eta,corre
      external h
      double precision h, erfc, sqrt, dlog

      double precision one,two,three,four,nine
      data one,two,three,four,nine/1.d0,2.d0,3.d0,4.d0,9.d0/
c
c check of initial values
c
      if (noncentr.lt. 0.) then
         ifault=1
         return
      endif

      if (df.le. 0.) then
         ifault=2
         return
      endif

      if (argument.le. 0.) then
         p= 0.
         return
      endif
c
c start calculation c

      mu2=noncentr/df
      s=(-one+ sqrt(one+four*argument*mu2/df))/(two*mu2)
      if (s.eq.one) then
         p=one/two
         return
      endif
      y=one-(one/s)
      help1= -h(y)
      help2= df*(s-one)*(s-one)*(one/s/two+mu2-(one/s)* help1) -
     +     dlog(one/s - two*help1/s/(one+two*mu2*s))

      if (variant.eq. 1) then
c                    "f"
c
c       calculating the improved _f_irst order approximation
c
         help2=help2+two*(one+three*mu2)*(one+three*mu2)/
     +        (one+two*mu2)**three/nine/df

      else ! variant = "s"
c
c        calculating the _s_econd order approximation c

         corre=-1.5*(one+four*mu2*s)/(one+two*mu2*s)**2+
     +        5.0*(one+three*mu2*s)**two/(one+two*mu2*s)**three/three

         eta=(one+two*mu2*s-two*help1-s-two*mu2*s*s)/
     +        (one+two*mu2*s-two*help1)

         corre=two*(corre+ two*(one+three*mu2*s)/(s-1)/(one +
     1        two*mu2*s)**two + three*eta/(s-one)**two/(one +
     2        two*mu2*s)-0.5*(one+two*h(eta))*eta*eta/
     3        (one+two*mu2*s)/(s-1)/(s-1))/df

         help2=help2+corre

      endif

c 10
      help2 = sqrt(help2)

      if (s.lt.one) help2 = -help2

      p=0.5*erfc(-help2/sqrt(two))

      end

