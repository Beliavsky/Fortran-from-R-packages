! Derived from locfit src/family.c, GPL-2-or-later.
module locfit_families
  use locfit_kinds, only : dp
  use locfit_constants
  use locfit_math, only : lf_exp, expit, logit, normal_cdf, gamma_p, beta_i, ptail
  implicit none
  private
  public :: default_link, valid_link, link_value, inverse_link, family_terms
  public :: cumulant_b2, cumulant_b3, cumulant_b4

contains

  pure integer function default_link(link, family) result(ans)
    integer, intent(in) :: link, family
    integer :: f
    f=iand(family,63)
    ans=link
    if(link==ldefau)then
      select case(f)
      case(tden,trat,thaz,tgamm,tgeom,tprob,tpois); ans=llog
      case(tcirc,tgaus,tcauc,trobt); ans=lident
      case(trbin,tlogt); ans=llogit
      end select
    else if(link==lcanon)then
      select case(f)
      case(tden,trat,thaz,tprob,tpois); ans=llog
      case(tgeom,tgamm); ans=linver
      case(tcirc,tgaus,tcauc,trobt); ans=lident
      case(trbin,tlogt); ans=llogit
      end select
    end if
  end function default_link

  pure logical function valid_link(link,family) result(ok)
    integer,intent(in)::link,family
    select case(iand(family,63))
    case(tden,trat,thaz); ok=(link==llog .or. link==lident)
    case(tgaus); ok=any(link==[lident,llog,llogit])
    case(trobt,tcauc,tcirc); ok=(link==lident)
    case(tlogt); ok=any(link==[llogit,lident,lasin])
    case(trbin); ok=(link==llogit)
    case(tgamm); ok=any(link==[llog,linver,lident])
    case(tgeom); ok=any(link==[llog,lident])
    case(tpois,tprob); ok=any(link==[llog,lsqrt,lident])
    case default; ok=.false.
    end select
  end function valid_link

  pure real(dp) function link_value(y,link) result(th)
    real(dp),intent(in)::y
    integer,intent(in)::link
    select case(link)
    case(lident); th=y
    case(llog); th=log(y)
    case(llogit); th=logit(y)
    case(linver); th=1.0_dp/y
    case(lsqrt); th=sqrt(abs(y))
    case(lasin); th=asin(sqrt(y))
    case default; th=0.0_dp
    end select
  end function link_value

  pure real(dp) function inverse_link(th,link) result(y)
    real(dp),intent(in)::th
    integer,intent(in)::link
    select case(link)
    case(lident); y=th
    case(llog); y=lf_exp(th)
    case(llogit); y=expit(th)
    case(linver); y=1.0_dp/th
    case(lsqrt); y=th*abs(th)
    case(lasin); y=sin(th)**2
    case(linit); y=0.0_dp
    case default; y=0.0_dp
    end select
  end function inverse_link

  pure subroutine robustify_terms(res,rs)
    real(dp),intent(inout)::res(llen)
    real(dp),intent(in)::rs
    real(dp)::sc,z
    sc=rs*huberc
    if(sc<=0.0_dp)return
    if(res(zlik)>-sc*sc/2.0_dp)then
      res(zlik)=res(zlik)/(sc*sc)
      res(zdll)=res(zdll)/(sc*sc)
      res(zddll)=res(zddll)/(sc*sc)
      return
    end if
    z=sqrt(max(0.0_dp,-2.0_dp*res(zlik)))
    if(z==0.0_dp)return
    res(zddll)=(-sc*res(zdll)**2/z**3+sc*res(zddll)/z)/(sc*sc)
    res(zdll)=res(zdll)/(z*sc)
    res(zlik)=0.5_dp-z/sc
  end subroutine robustify_terms

  pure subroutine family_terms(th,y,family,link,res,status,censored,prior_weight,robust_scale)
    real(dp),intent(in)::th,y
    integer,intent(in)::family,link
    real(dp),intent(out)::res(llen)
    integer,intent(out)::status
    logical,intent(in),optional::censored
    real(dp),intent(in),optional::prior_weight,robust_scale
    logical::cens
    real(dp)::w,rs,mean,z,pz,dpv,wmu,pt,dq,dg,p,wp,s2y,yy,sw
    integer::f
    res=0.0_dp
    cens=.false.; if(present(censored))cens=censored
    w=1.0_dp; if(present(prior_weight))w=prior_weight
    rs=1.0_dp; if(present(robust_scale))rs=robust_scale
    f=iand(family,63)
    status=lf_ok
    mean=inverse_link(th,link)
    res(zmean)=mean

    select case(f)
    case(tden,trat,thaz)
      if(.not.cens)then
        res(zlik)=w*th; res(zdll)=w; res(zddll)=w
      end if

    case(tgaus)
      if(link==linit)then
        res(zdll)=w*y; return
      end if
      z=y-mean
      if(cens)then
        if(link/=lident)then; status=lf_lnk; return; end if
        pz=normal_cdf(-z)
        if(pz<=tiny(1.0_dp))then; status=lf_badp; return; end if
        if(z>6.0_dp)then
          dpv=ptail(-z)/s2pi
        else
          dpv=exp(-z*z/2.0_dp)/(pz*s2pi)
        end if
        res(zlik)=w*log(pz); res(zdll)=w*dpv; res(zddll)=w*dpv*(dpv-z)
      else
        res(zlik)=-w*z*z/2.0_dp
        select case(link)
        case(lident); res(zdll)=w*z; res(zddll)=w
        case(llog); res(zdll)=w*z*mean; res(zddll)=w*mean*mean
        case(llogit); res(zdll)=w*z*mean*(1.0_dp-mean); res(zddll)=w*mean**2*(1.0_dp-mean)**2
        case default; status=lf_lnk
        end select
      end if

    case(trobt)
      if(link==linit)then; res(zdll)=w*y; return; end if
      if(rs<=0.0_dp)then; status=lf_badp; return; end if
      sw=merge(1.0_dp,sqrt(w),w==1.0_dp)
      z=sw*(y-mean)/rs
      if(abs(z)<huberc)then
        res(zlik)=-z*z/2.0_dp
      else
        res(zlik)=huberc*(huberc/2.0_dp-abs(z))
      end if
      if(z < -huberc)then
        res(zdll)=-sw*huberc/rs; res(zddll)=0.0_dp
      else if(z > huberc)then
        res(zdll)=sw*huberc/rs; res(zddll)=0.0_dp
      else
        res(zdll)=sw*z/rs; res(zddll)=w/(rs*rs)
      end if

    case(tcauc)
      if(link/=lident .or. rs<=0.0_dp)then; status=lf_lnk; return; end if
      z=w*(y-th)/rs
      res(zlik)=-log(1.0_dp+z*z)
      res(zdll)=2.0_dp*w*z/(rs*(1.0_dp+z*z))
      res(zddll)=2.0_dp*w*w*(1.0_dp-z*z)/(rs*rs*(1.0_dp+z*z)**2)

    case(trbin)
      if(link==linit)then; res(zdll)=y; return; end if
      if(y<0.0_dp .or. y>w)return
      p=mean
      if(th<0.0_dp)then
        res(zlik)=th*y-w*log(1.0_dp+exp(th))
      else
        res(zlik)=th*(y-w)-w*log(1.0_dp+exp(-th))
      end if
      if(y>0.0_dp)res(zlik)=res(zlik)-y*log(y/w)
      if(y<w)res(zlik)=res(zlik)-(w-y)*log(1.0_dp-y/w)
      res(zdll)=y-w*p; res(zddll)=w*p*(1.0_dp-p)
      if(-res(zlik)>huberc*huberc/2.0_dp)then
        s2y=sqrt(-2.0_dp*res(zlik)); res(zlik)=huberc*(huberc/2.0_dp-s2y)
        res(zdll)=res(zdll)*huberc/s2y
        res(zddll)=huberc/s2y*(res(zddll)-w*p*(1.0_dp-p)/(s2y*s2y))
      end if

    case(tlogt)
      if(link==linit)then
        res(zdll)=min(w,max(0.0_dp,y)); return
      end if
      p=mean; wp=w*p
      select case(link)
      case(lident)
        if((p<=0.0_dp .and. y>0.0_dp) .or. (p>=1.0_dp .and. y<w))then; status=lf_badp; return; end if
        if(y>0.0_dp)then; res(zlik)=res(zlik)+y*log(wp/y); res(zdll)=res(zdll)+y/p; res(zddll)=res(zddll)+y/(p*p); end if
        if(y<w)then; res(zlik)=res(zlik)+(w-y)*log((w-wp)/(w-y)); res(zdll)=res(zdll)-(w-y)/(1.0_dp-p); &
          res(zddll)=res(zddll)+(w-y)/(1.0_dp-p)**2; end if
      case(llogit)
        if(y<0.0_dp .or. y>w)return
        if(th<0.0_dp)then
          res(zlik)=th*y-w*log(1.0_dp+exp(th))
        else
          res(zlik)=th*(y-w)-w*log(1.0_dp+exp(-th))
        end if
        if(y>0.0_dp)res(zlik)=res(zlik)-y*log(y/w)
        if(y<w)res(zlik)=res(zlik)-(w-y)*log(1.0_dp-y/w)
        res(zdll)=y-wp; res(zddll)=wp*(1.0_dp-p)
      case(lasin)
        if ((p<=0.0_dp .and. y>0.0_dp) .or. (p>=1.0_dp .and. y<w) .or. &
            th<0.0_dp .or. th>pi/2.0_dp) then
          status=lf_badp
          return
        end if
        if(y>0.0_dp)then; res(zdll)=res(zdll)+2.0_dp*y*sqrt((1.0_dp-p)/p); res(zlik)=res(zlik)+y*log(wp/y); end if
        if(y<w)then; res(zdll)=res(zdll)-2.0_dp*(w-y)*sqrt(p/(1.0_dp-p)); res(zlik)=res(zlik)+(w-y)*log((w-wp)/(w-y)); end if
        res(zddll)=4.0_dp*w
      case default; status=lf_lnk
      end select

    case(tpois,tprob)
      if(link==linit)then; res(zdll)=max(y,0.0_dp); return; end if
      wmu=w*mean
      if(cens)then
        if(y<=0.0_dp)return
        pt=gamma_p(wmu,y)
        if(pt<=tiny(1.0_dp) .or. wmu<=0.0_dp)then; status=lf_badp; return; end if
        dpv=exp((y-1.0_dp)*log(wmu)-wmu-log_gamma(y))/pt
        dq=dpv*((y-1.0_dp)/wmu-1.0_dp)
        res(zlik)=log(pt)
        select case(link)
        case(llog); res(zdll)=dpv*wmu; res(zddll)=-(dq-dpv*dpv)*wmu*wmu-dpv*wmu
        case(lident); res(zdll)=dpv*w; res(zddll)=-(dq-dpv*dpv)*w*w
        case(lsqrt); res(zdll)=dpv*2.0_dp*w*th; res(zddll)=-(dq-dpv*dpv)*(4.0_dp*w*w*mean)-2.0_dp*dpv*w
        case default; status=lf_lnk
        end select
      else
        select case(link)
        case(llog)
          if(y<0.0_dp)return
          res(zlik)=y-wmu; res(zdll)=y-wmu; if(y>0.0_dp)res(zlik)=res(zlik)+y*(th-log(y/w)); res(zddll)=wmu
        case(lident)
          if(mean<=0.0_dp .and. y>0.0_dp)then; status=lf_badp; return; end if
          res(zlik)=y-wmu; res(zdll)=-w
          if(y>0.0_dp)then; res(zlik)=res(zlik)+y*log(wmu/y); res(zdll)=res(zdll)+y/mean; res(zddll)=y/(mean*mean); end if
        case(lsqrt)
          if(mean<=0.0_dp .and. y>0.0_dp)then; status=lf_badp; return; end if
          res(zlik)=y-wmu; res(zdll)=-2.0_dp*w*th; res(zddll)=2.0_dp*w
          if (y>0.0_dp) then
            res(zlik)=res(zlik)+y*log(wmu/y)
            res(zdll)=res(zdll)+2.0_dp*y/th
            res(zddll)=res(zddll)+2.0_dp*y/mean
          end if
        case default; status=lf_lnk
        end select
      end if

    case(tgamm)
      if(link==linit)then; res(zdll)=max(y,0.0_dp); return; end if
      if(mean<=0.0_dp .and. y>0.0_dp)then; status=lf_badp; return; end if
      if(cens)then
        if(y<=0.0_dp)return
        select case(link)
        case(llog)
          pt=1.0_dp-gamma_p(y/mean,w); if(pt<=tiny(1.0_dp))then;status=lf_badp;return;end if
          dg=exp((w-1.0_dp)*log(y/mean)-y/mean-log_gamma(w))
          res(zlik)=log(pt); res(zdll)=y*dg/(mean*pt); res(zddll)=dg*(w*y/mean-y*y/(mean*mean))/pt+res(zdll)**2
        case(linver)
          pt=1.0_dp-gamma_p(th*y,w); if(pt<=tiny(1.0_dp))then;status=lf_badp;return;end if
          dg=exp((w-1.0_dp)*log(th*y)-th*y-log_gamma(w))
          res(zlik)=log(pt); res(zdll)=-y*dg/pt; res(zddll)=dg*y*((w-1.0_dp)*mean-y)/pt+res(zdll)**2
        case default; status=lf_lnk
        end select
      else
        select case(link)
        case(llog)
          res(zlik)=-y/mean+w*(1.0_dp-th); if(y>0.0_dp)res(zlik)=res(zlik)+w*log(y/w)
          res(zdll)=y/mean-w; res(zddll)=y/mean
        case(linver)
          res(zlik)=-y/mean+w-w*log(mean); if(y>0.0_dp)res(zlik)=res(zlik)+w*log(y/w)
          res(zdll)=-y+w*mean; res(zddll)=w*mean*mean
        case(lident)
          res(zlik)=-y/mean+w-w*log(mean); if(y>0.0_dp)res(zlik)=res(zlik)+w*log(y/w)
          res(zdll)=(y-mean)/(mean*mean); res(zddll)=w/(mean*mean)
        case default; status=lf_lnk
        end select
      end if

    case(tgeom)
      if(link==linit)then; res(zdll)=max(y,0.0_dp); return; end if
      p=1.0_dp/(1.0_dp+mean)
      if(cens)then
        if(y<=0.0_dp)return
        pt=1.0_dp-beta_i(p,w,y); if(pt<=tiny(1.0_dp))then;status=lf_badp;return;end if
        dpv=-exp(log_gamma(w+y)-log_gamma(w)-log_gamma(y)+(y-1.0_dp)*th+(w+y-2.0_dp)*log(p))/pt
        dq=((w-1.0_dp)/p-(y-1.0_dp)/(1.0_dp-p))*dpv
        res(zlik)=log(pt); res(zdll)=-dpv*p*(1.0_dp-p)
        res(zddll)=-((dq-dpv*dpv)*p*p*(1.0_dp-p)**2+dpv*(1.0_dp-2.0_dp*p)*p*(1.0_dp-p))
      else
        res(zlik)=(y+w)*log((y/w+1.0_dp)/(mean+1.0_dp)); if(y>0.0_dp)res(zlik)=res(zlik)+y*log(w*mean/y)
        select case(link)
        case(llog); res(zdll)=(y-w*mean)*p; res(zddll)=(y+w)*p*(1.0_dp-p)
        case(lident); res(zdll)=(y-w*mean)/(mean*(1.0_dp+mean)); res(zddll)=w/(mean*(1.0_dp+mean))
        case default; status=lf_lnk
        end select
      end if

    case(tweib)
      yy=y**w
      if(link==linit)then; res(zdll)=max(yy,0.0_dp); return; end if
      if(mean<=0.0_dp)then; status=lf_badp; return; end if
      if(cens)then
        res(zlik)=-yy/mean; res(zdll)=yy/mean; res(zddll)=yy/mean
      else
        res(zlik)=1.0_dp-yy/mean-th; if(yy>0.0_dp)res(zlik)=res(zlik)+log(w*yy)
        res(zdll)=-1.0_dp+yy/mean; res(zddll)=yy/mean
      end if

    case(tcirc)
      if(link==linit)then; res(zdll)=w*sin(y); res(zlik)=w*cos(y); return; end if
      res(zdll)=w*sin(y-mean); res(zddll)=w*cos(y-mean); res(zlik)=res(zddll)-w

    case default
      status=lf_fam
    end select
    if(status==lf_ok .and. link/=linit .and. iand(family,128)==128)call robustify_terms(res,rs)
  end subroutine family_terms

  pure real(dp) function cumulant_b2(th,family,w) result(v)
    real(dp),intent(in)::th,w; integer,intent(in)::family; real(dp)::y
    select case(iand(family,63))
    case(tgaus); v=w
    case(tpois); v=w*lf_exp(th)
    case(tlogt); y=expit(th); v=w*y*(1.0_dp-y)
    case default; v=0.0_dp
    end select
  end function cumulant_b2

  pure real(dp) function cumulant_b3(th,family,w) result(v)
    real(dp),intent(in)::th,w; integer,intent(in)::family; real(dp)::y
    select case(iand(family,63))
    case(tgaus); v=0.0_dp
    case(tpois); v=w*lf_exp(th)
    case(tlogt); y=expit(th); v=w*y*(1.0_dp-y)*(1.0_dp-2.0_dp*y)
    case default; v=0.0_dp
    end select
  end function cumulant_b3

  pure real(dp) function cumulant_b4(th,family,w) result(v)
    real(dp),intent(in)::th,w; integer,intent(in)::family; real(dp)::y
    select case(iand(family,63))
    case(tgaus); v=0.0_dp
    case(tpois); v=w*lf_exp(th)
    case(tlogt); y=expit(th); y=y*(1.0_dp-y); v=w*y*(1.0_dp-6.0_dp*y)
    case default; v=0.0_dp
    end select
  end function cumulant_b4

end module locfit_families
