module suppdists_johnson
   use suppdists_kinds, only : dp
   use suppdists_special, only : normal_pdf, normal_cdf, normal_quantile, randn
   use suppdists_stats, only : dist_stats
   implicit none
   private
   integer, parameter, public :: johnson_sn=1, johnson_sl=2, johnson_su=3, johnson_sb=4
   type, public :: johnson_parms
      real(dp) :: gamma=0.0_dp, delta=1.0_dp, xi=0.0_dp, lambda=1.0_dp
      integer :: family=johnson_sn
   end type johnson_parms
   public :: djohnson, pjohnson, qjohnson, rjohnson, sjohnson, johnson_fit_quantiles, johnson_fit_moments
contains
   pure real(dp) function zjohnson(x,p) result(z)
      real(dp), intent(in) :: x
      type(johnson_parms), intent(in) :: p
      real(dp) :: u
      u=(x-p%xi)/p%lambda
      select case(p%family)
      case(johnson_sn)
      case(johnson_sl)
         u=log(u)
      case(johnson_su)
         u=asinh(u)
      case(johnson_sb)
         u=log(u/(1.0_dp-u))
      end select
      z=p%gamma+p%delta*u
   end function zjohnson

   pure real(dp) function pjohnson(x,p) result(prob)
      real(dp), intent(in) :: x
      type(johnson_parms), intent(in) :: p
      real(dp) :: z,u
      u=(x-p%xi)/p%lambda
      if(p%family==johnson_sb)then
         if(u<=0.0_dp)then;prob=0.0_dp;return;end if
         if(u>=1.0_dp)then;prob=1.0_dp;return;end if
      else if(p%family==johnson_sl .and. u<=0.0_dp)then
         prob=0.0_dp;return
      end if
      z=zjohnson(x,p)
      prob=normal_cdf(z)
   end function pjohnson

   pure real(dp) function qjohnson(prob,p) result(x)
      real(dp), intent(in) :: prob
      type(johnson_parms), intent(in) :: p
      real(dp) :: z,u,e
      z=normal_quantile(prob); u=(z-p%gamma)/p%delta
      select case(p%family)
      case(johnson_sn)
      case(johnson_sl)
         u=exp(u)
      case(johnson_su)
         u=sinh(u)
      case(johnson_sb)
         if(u>=0.0_dp)then;e=exp(-u);u=1.0_dp/(1.0_dp+e)
         else;e=exp(u);u=e/(1.0_dp+e);end if
      end select
      x=p%xi+p%lambda*u
   end function qjohnson

   pure real(dp) function djohnson(x,p) result(f)
      real(dp), intent(in) :: x
      type(johnson_parms), intent(in) :: p
      real(dp) :: z,u,jac
      u=(x-p%xi)/p%lambda
      if(p%family==johnson_sl .and. u<=0.0_dp)then;f=0.0_dp;return;end if
      if(p%family==johnson_sb .and. (u<=0.0_dp .or. u>=1.0_dp))then;f=0.0_dp;return;end if
      z=zjohnson(x,p)
      select case(p%family)
      case(johnson_sn); jac=p%delta/p%lambda
      case(johnson_sl); jac=p%delta/(p%lambda*u)
      case(johnson_su); jac=p%delta/(p%lambda*sqrt(1.0_dp+u*u))
      case(johnson_sb); jac=p%delta/(p%lambda*u*(1.0_dp-u))
      end select
      f=normal_pdf(z)*abs(jac)
   end function djohnson

   real(dp) function rjohnson(p) result(x)
      type(johnson_parms), intent(in) :: p
      real(dp) :: z,u,e
      z=randn(); u=(z-p%gamma)/p%delta
      select case(p%family)
      case(johnson_sn)
      case(johnson_sl); u=exp(u)
      case(johnson_su); u=sinh(u)
      case(johnson_sb)
         if(u>=0.0_dp)then;e=exp(-u);u=1.0_dp/(1.0_dp+e)
         else;e=exp(u);u=e/(1.0_dp+e);end if
      end select
      x=p%xi+p%lambda*u
   end function rjohnson

   function sjohnson(p) result(s)
      type(johnson_parms), intent(in) :: p
      type(dist_stats) :: s
      integer, parameter :: n=4001
      real(dp) :: z,x,w,sumw,m1,m2,m3,m4,dz,bestf,f
      integer :: i
      dz=16.0_dp/real(n-1,dp); sumw=0.0_dp; m1=0.0_dp
      do i=1,n
         z=-8.0_dp+real(i-1,dp)*dz
         w=normal_pdf(z); if(i==1 .or. i==n)w=0.5_dp*w
         x=x_from_z(z,p); sumw=sumw+w; m1=m1+w*x
      end do
      m1=m1/sumw; m2=0.0_dp;m3=0.0_dp;m4=0.0_dp
      do i=1,n
         z=-8.0_dp+real(i-1,dp)*dz
         w=normal_pdf(z); if(i==1 .or. i==n)w=0.5_dp*w
         x=x_from_z(z,p); m2=m2+w*(x-m1)**2; m3=m3+w*(x-m1)**3; m4=m4+w*(x-m1)**4
      end do
      s%mean=m1;s%variance=m2/sumw;s%third_central=m3/sumw;s%fourth_central=m4/sumw
      s%median=qjohnson(0.5_dp,p)
      s%mode=s%median;bestf=djohnson(s%mode,p)
      do i=0,500
         x=qjohnson((real(i,dp)+0.5_dp)/501.0_dp,p); f=djohnson(x,p)
         if(f>bestf)then;bestf=f;s%mode=x;end if
      end do
   contains
      pure real(dp) function x_from_z(z,p) result(x)
         real(dp),intent(in)::z; type(johnson_parms),intent(in)::p
         real(dp)::u,e
         u=(z-p%gamma)/p%delta
         select case(p%family)
         case(johnson_sn)
         case(johnson_sl);u=exp(u)
         case(johnson_su);u=sinh(u)
         case(johnson_sb)
            if(u>=0)then;e=exp(-u);u=1/(1+e);else;e=exp(u);u=e/(1+e);end if
         end select
         x=p%xi+p%lambda*u
      end function x_from_z
   end function sjohnson

   function johnson_fit_quantiles(x05,x206,x50,x794,x95) result(p)
      real(dp), intent(in) :: x05,x206,x50,x794,x95
      type(johnson_parms) :: p
      real(dp), parameter :: zn=1.64485363_dp
      real(dp) :: t,tu,tb,tbu,b,a,uu(5),yy(5),xm,ym,sxx,sxy
      integer :: i
      t=(x95-x50)/(x50-x05)
      tu=(x95-x05)/(x794-x206)
      tb=0.5_dp*((x794-x50)*(x95-x05)/((x95-x794)*(x50-x05))+ &
         (x206-x50)*(x05-x95)/((x05-x206)*(x50-x95)))
      tbu=tb/tu
      if(abs(abs(tbu)-1.0_dp)<0.1_dp .and. abs(abs(t-1.0_dp))<0.1_dp)then
         p%family=johnson_sn;p%delta=1.0_dp;p%gamma=0.0_dp
      else if(abs(abs(tbu)-1.0_dp)<0.1_dp)then
         p%family=johnson_sl;p%delta=zn/log(t);p%gamma=0.0_dp
      else if(tbu>1.0_dp)then
         p%family=johnson_sb;tb=0.5_dp*tb;b=tb+sqrt(max(0.0_dp,tb*tb-1.0_dp))
         p%delta=zn/(2.0_dp*log(b));b=b*b;a=(t-b)/(1.0_dp-t*b);p%gamma=-p%delta*log(a)
      else
         p%family=johnson_su;tu=0.5_dp*tu;b=tu+sqrt(max(0.0_dp,tu*tu-1.0_dp))
         p%delta=zn/(2.0_dp*log(b));b=b*b;a=(1.0_dp-t*b)/(t-b);p%gamma=-0.5_dp*p%delta*log(a)
      end if
      uu=[-zn,-zn/2.0_dp,0.0_dp,zn/2.0_dp,zn]
      yy=[x05,x206,x50,x794,x95]
      do i=1,5
         select case(p%family)
         case(johnson_sn)
         case(johnson_sl); uu(i)=exp(uu(i)/p%delta)
         case(johnson_su); uu(i)=sinh((uu(i)-p%gamma)/p%delta)
         case(johnson_sb)
            a=exp((uu(i)-p%gamma)/p%delta);uu(i)=a/(1.0_dp+a)
         end select
      end do
      xm=sum(uu)/5.0_dp;ym=sum(yy)/5.0_dp
      sxx=sum((uu-xm)**2);sxy=sum((uu-xm)*(yy-ym))
      p%lambda=sxy/sxx;p%xi=ym-p%lambda*xm
   end function johnson_fit_quantiles

   function johnson_fit_moments(mean,sd,skew,excess) result(p)
      real(dp), intent(in) :: mean,sd,skew,excess
      type(johnson_parms) :: p
      real(dp) :: b1,b2,x,y,u,w,b2log,test,best,err,g,d,stepg,stepd
      real(dp) :: mu0,var0,sk0,ku0
      integer :: fam,ig,id,it
      b1=skew*skew; b2=excess+3.0_dp
      if(abs(skew)<=0.1_dp .and. abs(excess)<=0.1_dp)then
         p%family=johnson_sn;p%gamma=0.0_dp;p%delta=1.0_dp;p%xi=mean;p%lambda=sd
         return
      end if
      x=0.5_dp*b1+1.0_dp;y=sqrt(max(0.0_dp,b1+0.25_dp*b1*b1));u=(x+y)**(1.0_dp/3.0_dp)
      w=u+1.0_dp/u-1.0_dp;b2log=w*w*(3.0_dp+w*(2.0_dp+w))-3.0_dp;test=b2log-b2
      if(abs(test)<0.1_dp .and. w>1.0_dp)then
         p%family=johnson_sl;p%lambda=1.0_dp;p%delta=1.0_dp/sqrt(log(w))
         p%gamma=0.5_dp*p%delta*log(w*(w-1.0_dp)/(sd*sd))
         p%xi=mean-sd/sqrt(w-1.0_dp)
         return
      end if
      if(test<=0.0_dp)then;fam=johnson_su;else;fam=johnson_sb;end if
      best=huge(1.0_dp);g=0.0_dp;d=1.0_dp
      do ig=-10,10
         do id=0,20
            call shape_moments(fam,0.5_dp*real(ig,dp),exp(-1.5_dp+0.15_dp*real(id,dp)),mu0,var0,sk0,ku0)
            err=(sk0-skew)**2+(ku0-excess)**2
            if(err<best)then
               best=err;g=0.5_dp*real(ig,dp);d=exp(-1.5_dp+0.15_dp*real(id,dp))
            end if
         end do
      end do
      stepg=0.5_dp;stepd=0.25_dp*d
      do it=1,12
         best=huge(1.0_dp)
         do ig=-2,2
            do id=-2,2
               call shape_moments(fam,g+stepg*real(ig,dp),max(0.05_dp,d+stepd*real(id,dp)),mu0,var0,sk0,ku0)
               err=(sk0-skew)**2+(ku0-excess)**2
               if(err<best)then
                  best=err;x=g+stepg*real(ig,dp);y=max(0.05_dp,d+stepd*real(id,dp))
               end if
            end do
         end do
         g=x;d=y;stepg=0.5_dp*stepg;stepd=0.5_dp*stepd
      end do
      call shape_moments(fam,g,d,mu0,var0,sk0,ku0)
      p%family=fam;p%gamma=g;p%delta=d
      p%lambda=sd/sqrt(max(var0,tiny(1.0_dp)))
      p%xi=mean-p%lambda*mu0
   contains
      subroutine shape_moments(family,gamma,delta,mu,var,sk,ku)
         integer,intent(in)::family
         real(dp),intent(in)::gamma,delta
         real(dp),intent(out)::mu,var,sk,ku
         integer,parameter::ng=1601
         real(dp)::z,uu,e,wgt,sw,m1,m2,m3,m4,dz
         integer::j
         dz=14.0_dp/real(ng-1,dp);sw=0.0_dp;m1=0.0_dp;m2=0.0_dp;m3=0.0_dp;m4=0.0_dp
         do j=1,ng
            z=-7.0_dp+real(j-1,dp)*dz;wgt=normal_pdf(z);if(j==1 .or. j==ng)wgt=0.5_dp*wgt
            uu=(z-gamma)/delta
            select case(family)
            case(johnson_su);uu=sinh(uu)
            case(johnson_sb)
               if(uu>=0.0_dp)then;e=exp(-uu);uu=1.0_dp/(1.0_dp+e)
               else;e=exp(uu);uu=e/(1.0_dp+e);end if
            end select
            sw=sw+wgt;m1=m1+wgt*uu;m2=m2+wgt*uu*uu;m3=m3+wgt*uu**3;m4=m4+wgt*uu**4
         end do
         m1=m1/sw;m2=m2/sw;m3=m3/sw;m4=m4/sw
         var=max(0.0_dp,m2-m1*m1);mu=m1
         if(var<=tiny(1.0_dp))then;sk=0.0_dp;ku=0.0_dp;return;end if
         sk=(m3-3.0_dp*m1*m2+2.0_dp*m1**3)/var**1.5_dp
         ku=(m4-4.0_dp*m1*m3+6.0_dp*m1*m1*m2-3.0_dp*m1**4)/(var*var)-3.0_dp
      end subroutine shape_moments
   end function johnson_fit_moments

end module suppdists_johnson
