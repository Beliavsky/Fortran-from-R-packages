module discretedists_family
   use discretedists_kinds, only : dp
   use discretedists_numerics, only : logistic, logit, normal_cdf, normal_quantile, mean_sample, variance_sample
   use discretedists_distributions
   use discretedists_estimators
   use compoissonreg_distributions, only : cmp_stats, vcmp
   implicit none
   private

   integer, parameter, public :: DD_FAM_BERG=1,DD_FAM_COMPO=2,DD_FAM_COMPO2=3,DD_FAM_DBH=4
   integer, parameter, public :: DD_FAM_DGEII=5,DD_FAM_DIKUM=6,DD_FAM_DLD=7,DD_FAM_DMOLBE=8
   integer, parameter, public :: DD_FAM_DPERKS=9,DD_FAM_DSPA=10,DD_FAM_GGEO=11,DD_FAM_HYPERPO=12
   integer, parameter, public :: DD_FAM_HYPERPO2=13,DD_FAM_POISXL=14

   type, public :: discrete_family_t
      integer :: id=0
      integer :: npar=0
      character(len=12) :: short_name=''
      character(len=12) :: mu_link='identity'
      character(len=12) :: sigma_link='identity'
   contains
      procedure :: density => family_density
      procedure :: cdf => family_cdf
      procedure :: score => family_score
      procedure :: curvature => family_curvature
      procedure :: deviance_increment => family_deviance
      procedure :: mean_variance => family_mean_variance
      procedure :: initial => family_initial
      procedure :: valid => family_valid
      procedure :: linkfun => family_linkfun
      procedure :: linkinv => family_linkinv
      procedure :: mu_eta => family_mu_eta
      procedure :: rqres => family_rqres
   end type discrete_family_t

   public :: berg,compo,compo2,dbh,dgeii,dikum,dld,dmolbe,dperks,dspa,ggeo,hyperpo,hyperpo2,poisxl

contains

   function make_family(id,name,npar,mu_link,sigma_link) result(f)
      integer,intent(in)::id,npar
      character(*),intent(in)::name,mu_link,sigma_link
      type(discrete_family_t)::f
      f%id=id;f%npar=npar;f%short_name=name;f%mu_link=mu_link;f%sigma_link=sigma_link
   end function make_family

   function berg(mu_link,sigma_link) result(f)
      character(*),intent(in),optional::mu_link,sigma_link
      type(discrete_family_t)::f
      character(len=12)::m,s
      m='sqrt';s='log';if(present(mu_link))m=mu_link;if(present(sigma_link))s=sigma_link
      f=make_family(DD_FAM_BERG,'BerG',2,m,s)
   end function berg
   function compo(mu_link,sigma_link) result(f)
      character(*),intent(in),optional::mu_link,sigma_link
      type(discrete_family_t)::f;character(len=12)::m,s
      m='log';s='log';if(present(mu_link))m=mu_link;if(present(sigma_link))s=sigma_link
      f=make_family(DD_FAM_COMPO,'COMPO',2,m,s)
   end function compo
   function compo2(mu_link,sigma_link) result(f)
      character(*),intent(in),optional::mu_link,sigma_link
      type(discrete_family_t)::f;character(len=12)::m,s
      m='log';s='identity';if(present(mu_link))m=mu_link;if(present(sigma_link))s=sigma_link
      f=make_family(DD_FAM_COMPO2,'COMPO2',2,m,s)
   end function compo2
   function dbh(mu_link) result(f)
      character(*),intent(in),optional::mu_link
      type(discrete_family_t)::f;character(len=12)::m
      m='logit';if(present(mu_link))m=mu_link
      f=make_family(DD_FAM_DBH,'DBH',1,m,'identity')
   end function dbh
   function dgeii(mu_link,sigma_link) result(f)
      character(*),intent(in),optional::mu_link,sigma_link
      type(discrete_family_t)::f;character(len=12)::m,s
      m='logit';s='log';if(present(mu_link))m=mu_link;if(present(sigma_link))s=sigma_link
      f=make_family(DD_FAM_DGEII,'DGEII',2,m,s)
   end function dgeii
   function dikum(mu_link,sigma_link) result(f)
      character(*),intent(in),optional::mu_link,sigma_link
      type(discrete_family_t)::f;character(len=12)::m,s
      m='log';s='log';if(present(mu_link))m=mu_link;if(present(sigma_link))s=sigma_link
      f=make_family(DD_FAM_DIKUM,'DIKUM',2,m,s)
   end function dikum
   function dld(mu_link) result(f)
      character(*),intent(in),optional::mu_link
      type(discrete_family_t)::f;character(len=12)::m
      m='log';if(present(mu_link))m=mu_link
      f=make_family(DD_FAM_DLD,'DLD',1,m,'identity')
   end function dld
   function dmolbe(mu_link,sigma_link) result(f)
      character(*),intent(in),optional::mu_link,sigma_link
      type(discrete_family_t)::f;character(len=12)::m,s
      m='log';s='log';if(present(mu_link))m=mu_link;if(present(sigma_link))s=sigma_link
      f=make_family(DD_FAM_DMOLBE,'DMOLBE',2,m,s)
   end function dmolbe
   function dperks(mu_link,sigma_link) result(f)
      character(*),intent(in),optional::mu_link,sigma_link
      type(discrete_family_t)::f;character(len=12)::m,s
      m='log';s='log';if(present(mu_link))m=mu_link;if(present(sigma_link))s=sigma_link
      f=make_family(DD_FAM_DPERKS,'DPERKS',2,m,s)
   end function dperks
   function dspa(mu_link,sigma_link) result(f)
      character(*),intent(in),optional::mu_link,sigma_link
      type(discrete_family_t)::f;character(len=12)::m,s
      m='log';s='logit';if(present(mu_link))m=mu_link;if(present(sigma_link))s=sigma_link
      f=make_family(DD_FAM_DSPA,'DsPA',2,m,s)
   end function dspa
   function ggeo(mu_link,sigma_link) result(f)
      character(*),intent(in),optional::mu_link,sigma_link
      type(discrete_family_t)::f;character(len=12)::m,s
      m='logit';s='log';if(present(mu_link))m=mu_link;if(present(sigma_link))s=sigma_link
      f=make_family(DD_FAM_GGEO,'GGEO',2,m,s)
   end function ggeo
   function hyperpo(mu_link,sigma_link) result(f)
      character(*),intent(in),optional::mu_link,sigma_link
      type(discrete_family_t)::f;character(len=12)::m,s
      m='log';s='log';if(present(mu_link))m=mu_link;if(present(sigma_link))s=sigma_link
      f=make_family(DD_FAM_HYPERPO,'HYPERPO',2,m,s)
   end function hyperpo
   function hyperpo2(mu_link,sigma_link) result(f)
      character(*),intent(in),optional::mu_link,sigma_link
      type(discrete_family_t)::f;character(len=12)::m,s
      m='log';s='log';if(present(mu_link))m=mu_link;if(present(sigma_link))s=sigma_link
      f=make_family(DD_FAM_HYPERPO2,'HYPERPO2',2,m,s)
   end function hyperpo2
   function poisxl(mu_link) result(f)
      character(*),intent(in),optional::mu_link
      type(discrete_family_t)::f;character(len=12)::m
      m='log';if(present(mu_link))m=mu_link
      f=make_family(DD_FAM_POISXL,'POISXL',1,m,'identity')
   end function poisxl

   real(dp) function family_density(self,y,mu,sigma,log_p) result(v)
      class(discrete_family_t),intent(in)::self
      real(dp),intent(in)::y,mu
      real(dp),intent(in),optional::sigma
      logical,intent(in),optional::log_p
      real(dp)::s;logical::lg
      s=1.0_dp;if(present(sigma))s=sigma;lg=.false.;if(present(log_p))lg=log_p
      select case(self%id)
      case(DD_FAM_BERG);v=dberg(y,mu,s,lg)
      case(DD_FAM_COMPO);v=dcompo(y,mu,s,lg)
      case(DD_FAM_COMPO2);v=dcompo2(y,mu,s,lg)
      case(DD_FAM_DBH);v=ddbh(y,mu,lg)
      case(DD_FAM_DGEII);v=ddgeii(y,mu,s,lg)
      case(DD_FAM_DIKUM);v=ddikum(y,mu,s,lg)
      case(DD_FAM_DLD);v=ddld(y,mu,lg)
      case(DD_FAM_DMOLBE);v=ddmolbe(y,mu,s,lg)
      case(DD_FAM_DPERKS);v=ddperks(y,mu,s,lg)
      case(DD_FAM_DSPA);v=ddspa(y,mu,s,lg)
      case(DD_FAM_GGEO);v=dggeo(y,mu,s,lg)
      case(DD_FAM_HYPERPO);v=dhyperpo(y,mu,s,lg)
      case(DD_FAM_HYPERPO2);v=dhyperpo2(y,mu,s,lg)
      case(DD_FAM_POISXL);v=dpoisxl(y,mu,lg)
      case default;v=merge(-huge(1.0_dp),0.0_dp,lg)
      end select
   end function family_density

   real(dp) function family_cdf(self,y,mu,sigma,lower_tail,log_p) result(v)
      class(discrete_family_t),intent(in)::self
      real(dp),intent(in)::y,mu
      real(dp),intent(in),optional::sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::s;logical::lo,lg
      s=1.0_dp;if(present(sigma))s=sigma;lo=.true.;if(present(lower_tail))lo=lower_tail
      lg=.false.;if(present(log_p))lg=log_p
      select case(self%id)
      case(DD_FAM_BERG);v=pberg(y,mu,s,lo,lg)
      case(DD_FAM_COMPO);v=pcompo(y,mu,s,lo,lg)
      case(DD_FAM_COMPO2);v=pcompo2(y,mu,s,lo,lg)
      case(DD_FAM_DBH);v=pdbh(y,mu,lo,lg)
      case(DD_FAM_DGEII);v=pdgeii(y,mu,s,lo,lg)
      case(DD_FAM_DIKUM);v=pdikum(y,mu,s,lo,lg)
      case(DD_FAM_DLD);v=pdld(y,mu,lo,lg)
      case(DD_FAM_DMOLBE);v=pdmolbe(y,mu,s,lo,lg)
      case(DD_FAM_DPERKS);v=pdperks(y,mu,s,lo,lg)
      case(DD_FAM_DSPA);v=pdspa(y,mu,s,lo,lg)
      case(DD_FAM_GGEO);v=pggeo(y,mu,s,lo,lg)
      case(DD_FAM_HYPERPO);v=phyperpo(y,mu,s,lo,lg)
      case(DD_FAM_HYPERPO2);v=phyperpo2(y,mu,s,lo,lg)
      case(DD_FAM_POISXL);v=ppoisxl(y,mu,lo,lg)
      case default;v=0.0_dp
      end select
   end function family_cdf

   subroutine family_score(self,y,mu,sigma,dmu,dsigma)
      class(discrete_family_t),intent(in)::self
      real(dp),intent(in)::y,mu
      real(dp),intent(in),optional::sigma
      real(dp),intent(out)::dmu
      real(dp),intent(out),optional::dsigma
      real(dp)::s,h1,h2,lm,lp,t1,t2,a,b1,b2,da,db1,db2,logz,mn,var,elf
      s=1.0_dp;if(present(sigma))s=sigma
      select case(self%id)
      case(DD_FAM_BERG)
         if(nint(y)==0)then
            dmu=2.0_dp*(s+1.0_dp)/((mu-s-1.0_dp)*(mu+s+1.0_dp))
            if(present(dsigma))dsigma=2.0_dp*mu/((s-mu+1.0_dp)*(s+mu+1.0_dp))
         else
            dmu=1.0_dp/mu+(y-1.0_dp)/(mu+s-1.0_dp)-(y+1.0_dp)/(mu+s+1.0_dp)
            if(present(dsigma))dsigma=(y-1.0_dp)/(mu+s-1.0_dp)-(y+1.0_dp)/(mu+s+1.0_dp)
         end if
      case(DD_FAM_COMPO)
         call cmp_stats(mu,s,logz=logz,mean=mn,var=var,elogfact=elf)
         dmu=(y-mn)/mu;if(present(dsigma))dsigma=-log_gamma(y+1.0_dp)+elf
      case(DD_FAM_DBH)
         dmu=(-y*y*mu+y*y-2.0_dp*y*mu+2.0_dp*y-mu)/(mu*(y+2.0_dp-mu*(y+1.0_dp)))
      case(DD_FAM_DLD)
         dmu=-y-1.0_dp/(1.0_dp+mu)+(exp(mu)-1.0_dp+2.0_dp*mu+y*exp(mu)-y+mu*y)/ &
            (mu*exp(mu)-2.0_dp*mu+exp(mu)+y*mu*exp(mu)-1.0_dp-mu*y)
      case(DD_FAM_DPERKS)
         dmu=1.0_dp/mu+1.0_dp/(mu+1.0_dp)-exp(s*y)/(1.0_dp+mu*exp(s*y)) &
            -exp(s*(y+1.0_dp))/(1.0_dp+mu*exp(s*(y+1.0_dp)))
         if(present(dsigma))dsigma=exp(s)/(exp(s)-1.0_dp)+y-mu*y*exp(s*y)/(1.0_dp+mu*exp(s*y)) &
            -mu*(y+1.0_dp)*exp(s*(y+1.0_dp))/(1.0_dp+mu*exp(s*(y+1.0_dp)))
      case(DD_FAM_DMOLBE)
         t1=(1.0_dp+y/mu)*exp(-y/mu);t2=(1.0_dp+(y+1.0_dp)/mu)*exp(-(y+1.0_dp)/mu)
         a=t1-t2;b1=1.0_dp-(1.0_dp-s)*t1;b2=1.0_dp-(1.0_dp-s)*t2
         da=y*y/mu**3*exp(-y/mu)-(y+1.0_dp)**2/mu**3*exp(-(y+1.0_dp)/mu)
         db1=-(1.0_dp-s)*y*y/mu**3*exp(-y/mu)
         db2=-(1.0_dp-s)*(y+1.0_dp)**2/mu**3*exp(-(y+1.0_dp)/mu)
         dmu=da/a-db1/b1-db2/b2
         if(present(dsigma))dsigma=1.0_dp/s-t1/b1-t2/b2
      case(DD_FAM_GGEO)
         dmu=y/mu-1.0_dp/(1.0_dp-mu) &
            -(y+1.0_dp)*(s-1.0_dp)*mu**y/((s-1.0_dp)*mu**(y+1.0_dp)+1.0_dp) &
            -y*(s-1.0_dp)*mu**(y-1.0_dp)/((s-1.0_dp)*mu**y+1.0_dp)
         if(present(dsigma))dsigma=1.0_dp/s-mu**(y+1.0_dp)/(1.0_dp-mu**(y+1.0_dp)*(1.0_dp-s)) &
            -mu**y/(1.0_dp-mu**y*(1.0_dp-s))
      case default
         h1=parameter_step(self,mu,.true.);lm=self%density(y,mu-h1,s,.true.);lp=self%density(y,mu+h1,s,.true.)
         dmu=(lp-lm)/(2.0_dp*h1)
         if(present(dsigma))then
            if(self%npar==2)then
               h2=parameter_step(self,s,.false.)
               lm=self%density(y,mu,s-h2,.true.);lp=self%density(y,mu,s+h2,.true.)
               dsigma=(lp-lm)/(2.0_dp*h2)
            else
               dsigma=0.0_dp
            end if
         end if
      end select
      if(present(dsigma).and.self%npar==1)dsigma=0.0_dp
   end subroutine family_score

   real(dp) function parameter_step(self,x,is_mu) result(h)
      class(discrete_family_t),intent(in)::self
      real(dp),intent(in)::x
      logical,intent(in)::is_mu
      logical::bounded
      bounded=(is_mu.and.(self%id==DD_FAM_DBH.or.self%id==DD_FAM_DGEII.or.self%id==DD_FAM_GGEO)) .or. &
         ((.not.is_mu).and.self%id==DD_FAM_DSPA)
      if(bounded)then
         h=max(1.0e-7_dp,min(1.0e-5_dp,0.2_dp*min(x,1.0_dp-x)))
      else
         h=max(1.0e-7_dp,min(1.0e-4_dp,0.2_dp*max(x,1.0e-6_dp)))
      end if
   end function parameter_step

   subroutine family_curvature(self,y,mu,sigma,d2mu,d2musigma,d2sigma)
      class(discrete_family_t),intent(in)::self
      real(dp),intent(in)::y,mu
      real(dp),intent(in),optional::sigma
      real(dp),intent(out)::d2mu
      real(dp),intent(out),optional::d2musigma,d2sigma
      real(dp)::s,a,b
      s=1.0_dp;if(present(sigma))s=sigma
      call self%score(y,mu,s,a,b)
      if(self%id==DD_FAM_BERG)then
         if(nint(y)==0)then
            d2mu=-4.0_dp*(s+1.0_dp)*mu/((mu-s-1.0_dp)**2*(mu+s+1.0_dp)**2)
            if(present(d2sigma))d2sigma=-4.0_dp*mu*(s+1.0_dp)/((s-mu+1.0_dp)**2*(s+mu+1.0_dp)**2)
            if(present(d2musigma))d2musigma=2.0_dp*((s+1.0_dp)**2+mu**2)/ &
               ((s-mu+1.0_dp)**2*(mu+s+1.0_dp)**2)
         else
            d2mu=(y+1.0_dp)/(mu+s+1.0_dp)**2+(1.0_dp-y)/(mu+s-1.0_dp)**2-1.0_dp/mu**2
            if(present(d2sigma))d2sigma=(y+1.0_dp)/(mu+s+1.0_dp)**2+(1.0_dp-y)/(mu+s-1.0_dp)**2
            if(present(d2musigma))d2musigma=2.0_dp*((mu+s)*(mu+s-2.0_dp*y)+1.0_dp)/ &
               ((mu+s+1.0_dp)**2*(mu+s-1.0_dp)**2)
         end if
      else
         d2mu=min(-1.0e-15_dp,-a*a)
         if(present(d2sigma))d2sigma=merge(min(-1.0e-15_dp,-b*b),0.0_dp,self%npar==2)
         if(present(d2musigma))d2musigma=merge(-a*b,0.0_dp,self%npar==2)
      end if
   end subroutine family_curvature

   real(dp) function family_deviance(self,y,mu,sigma,weight) result(v)
      class(discrete_family_t),intent(in)::self
      real(dp),intent(in)::y,mu
      real(dp),intent(in),optional::sigma,weight
      real(dp)::s,w
      s=1.0_dp;if(present(sigma))s=sigma;w=1.0_dp;if(present(weight))w=weight
      v=-2.0_dp*w*self%density(y,mu,s,.true.)
   end function family_deviance

   subroutine family_mean_variance(self,mu,sigma,mean,var)
      class(discrete_family_t),intent(in)::self
      real(dp),intent(in)::mu
      real(dp),intent(in),optional::sigma
      real(dp),intent(out)::mean,var
      real(dp)::s,p,m1,m2,x,lambda,nu
      integer::k
      s=1.0_dp;if(present(sigma))s=sigma
      select case(self%id)
      case(DD_FAM_BERG);mean=mu;var=mu*s
      case(DD_FAM_COMPO)
         mean=compo_mean_exact(mu,s);var=compo_variance_exact(mu,s)
      case(DD_FAM_COMPO2)
         call mu_phi_2_lambda_nu_compo2(mu,s,lambda,nu);mean=mu;var=vcmp(lambda,nu)
      case(DD_FAM_DBH)
         mean=-log(1.0_dp-mu)/mu-1.0_dp
         var=mu/2.0_dp*((2.0_dp*(mu-3.0_dp))/(mu*(mu-1.0_dp))+6.0_dp*log(1.0_dp-mu)/mu**2)-mean**2
      case(DD_FAM_DLD)
         mean=exp(-mu)*(2.0_dp*mu-mu*exp(-mu)+1.0_dp-exp(-mu))/((1.0_dp+mu)*(1.0_dp-exp(-mu))**2)
         var=exp(-mu)*(mu*exp(-2.0_dp*mu)-3.0_dp*mu**2*exp(-mu)+3.0_dp*mu+2.0_dp*mu**2 &
            -2.0_dp*exp(-mu)+1.0_dp+exp(-2.0_dp*mu)-4.0_dp*mu*exp(-mu))/ &
            ((1.0_dp+mu)**2*(1.0_dp-exp(-mu))**4)
      case(DD_FAM_HYPERPO);call mean_var_hp(mu,s,mean,var)
      case(DD_FAM_HYPERPO2);call mean_var_hp2(mu,s,mean,var)
      case default
         m1=0.0_dp;m2=0.0_dp
         do k=0,100000
            x=real(k,dp);p=self%density(x,mu,s,.false.)
            m1=m1+x*p;m2=m2+x*x*p
            if(k>100.and.p<1.0e-13_dp)exit
         end do
         mean=m1;var=max(0.0_dp,m2-m1*m1)
      end select
   end subroutine family_mean_variance

   subroutine family_initial(self,y,mu,sigma,status)
      class(discrete_family_t),intent(in)::self
      real(dp),intent(in)::y(:)
      real(dp),intent(out)::mu
      real(dp),intent(out),optional::sigma
      integer,intent(out),optional::status
      real(dp)::p2(2);integer::st
      st=0
      select case(self%id)
      case(DD_FAM_BERG)
         mu=max(mean_sample(y),1.0e-6_dp);p2(2)=max(variance_sample(y)/mu,1.0e-6_dp)
      case(DD_FAM_COMPO);p2=estim_mu_sigma_compo(y);mu=p2(1)
      case(DD_FAM_COMPO2);mu=max(mean_sample(y),1.0e-6_dp);p2(2)=0.0_dp
      case(DD_FAM_DBH);mu=estim_mu_dbh(y,st)
      case(DD_FAM_DGEII);p2=estim_mu_sigma_dgeii(y,st);mu=p2(1)
      case(DD_FAM_DIKUM);p2=estim_mu_sigma_dikum(y,st);mu=p2(1)
      case(DD_FAM_DLD);mu=estim_mu_dld(y,st)
      case(DD_FAM_DMOLBE);p2=estim_mu_sigma_dmolbe(y,st);mu=p2(1)
      case(DD_FAM_DPERKS);mu=1.0_dp;p2(2)=1.0_dp
      case(DD_FAM_DSPA);p2=estim_mu_sigma_dspa(y,st);mu=p2(1)
      case(DD_FAM_GGEO);p2=estim_mu_sigma_ggeo(y,st);mu=p2(1)
      case(DD_FAM_HYPERPO);p2=estim_mu_sigma_hyperpo(y,st);mu=p2(1)
      case(DD_FAM_HYPERPO2);mu=max(mean_sample(y),1.0e-6_dp);p2(2)=1.0_dp
      case(DD_FAM_POISXL);mu=estim_mu_poisxl(y,st)
      case default;mu=max(mean_sample(y),1.0e-6_dp);p2(2)=1.0_dp;st=1
      end select
      if(self%id==DD_FAM_BERG)mu=max(mean_sample(y),1.0e-6_dp)
      if(present(sigma))then
         if(self%npar==2)sigma=p2(2)
         if(self%npar==1)sigma=0.0_dp
      end if
      if(present(status))status=st
   end subroutine family_initial

   logical function family_valid(self,mu,sigma) result(ok)
      class(discrete_family_t),intent(in)::self
      real(dp),intent(in)::mu
      real(dp),intent(in),optional::sigma
      real(dp)::s
      s=1.0_dp;if(present(sigma))s=sigma
      select case(self%id)
      case(DD_FAM_DBH,DD_FAM_DGEII,DD_FAM_GGEO);ok=mu>0.0_dp.and.mu<1.0_dp
      case default;ok=mu>0.0_dp
      end select
      if(self%npar==2)then
         if(self%id==DD_FAM_COMPO2)then;ok=ok
         else if(self%id==DD_FAM_DSPA)then;ok=ok.and.s>0.0_dp.and.s<1.0_dp
         else;ok=ok.and.s>0.0_dp
         end if
      end if
   end function family_valid

   real(dp) function family_linkfun(self,x,which) result(v)
      class(discrete_family_t),intent(in)::self
      real(dp),intent(in)::x
      integer,intent(in),optional::which
      character(len=12)::link
      link=self%mu_link;if(present(which))then;if(which==2)link=self%sigma_link;end if
      select case(trim(link))
      case('identity');v=x
      case('log');v=log(x)
      case('sqrt');v=sqrt(x)
      case('logit');v=logit(x)
      case('probit');v=normal_quantile(x)
      case('cloglog');v=log(-log(1.0_dp-x))
      case('cauchit');v=tan(acos(-1.0_dp)*(x-0.5_dp))
      case('inverse');v=1.0_dp/x
      case default;v=x
      end select
   end function family_linkfun

   real(dp) function family_linkinv(self,eta,which) result(v)
      class(discrete_family_t),intent(in)::self
      real(dp),intent(in)::eta
      integer,intent(in),optional::which
      character(len=12)::link
      link=self%mu_link;if(present(which))then;if(which==2)link=self%sigma_link;end if
      select case(trim(link))
      case('identity');v=eta
      case('log');v=exp(eta)
      case('sqrt');v=eta*eta
      case('logit');v=logistic(eta)
      case('probit');v=normal_cdf(eta)
      case('cloglog');v=1.0_dp-exp(-exp(eta))
      case('cauchit');v=0.5_dp+atan(eta)/acos(-1.0_dp)
      case('inverse');v=1.0_dp/eta
      case default;v=eta
      end select
   end function family_linkinv

   real(dp) function family_mu_eta(self,eta,which) result(v)
      class(discrete_family_t),intent(in)::self
      real(dp),intent(in)::eta
      integer,intent(in),optional::which
      real(dp)::x
      character(len=12)::link
      link=self%mu_link;if(present(which))then;if(which==2)link=self%sigma_link;end if
      x=self%linkinv(eta,which)
      select case(trim(link))
      case('identity');v=1.0_dp
      case('log');v=x
      case('sqrt');v=2.0_dp*eta
      case('logit');v=x*(1.0_dp-x)
      case('probit');v=exp(-0.5_dp*eta*eta)/sqrt(2.0_dp*acos(-1.0_dp))
      case('cloglog');v=exp(eta-exp(eta))
      case('cauchit');v=1.0_dp/(acos(-1.0_dp)*(1.0_dp+eta*eta))
      case('inverse');v=-1.0_dp/(eta*eta)
      case default;v=1.0_dp
      end select
   end function family_mu_eta

   real(dp) function family_rqres(self,y,mu,sigma,u) result(r)
      class(discrete_family_t),intent(in)::self
      integer,intent(in)::y
      real(dp),intent(in)::mu
      real(dp),intent(in),optional::sigma,u
      real(dp)::s,a,b,z,uu
      s=1.0_dp;if(present(sigma))s=sigma
      a=self%cdf(real(y-1,dp),mu,s,.true.,.false.);b=self%cdf(real(y,dp),mu,s,.true.,.false.)
      if(present(u))then;uu=max(0.0_dp,min(1.0_dp,u));else;call random_number(uu);end if
      z=a+uu*(b-a);z=max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),z));r=normal_quantile(z)
   end function family_rqres

end module discretedists_family
