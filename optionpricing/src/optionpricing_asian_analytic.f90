! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from OptionPricing 0.1.2 by Wolfgang Hormann and Kemal Dingec.
module optionpricing_asian_analytic
   use optionpricing_kinds, only : dp
   use optionpricing_math, only : normal_cdf, normal_pdf, adaptive_integral, bisection_root
   use optionpricing_types, only : moments_result
   use optionpricing_linalg, only : identity_matrix, lower_ones
   implicit none
   private
   public :: bs_a, covariance_conditional_log_prices, conditional_average_moments
   public :: eval_ecv_a, asian_call_app_lord, eval_ecv, find_bcv, eval_lb
   public :: eval_equad, eval_eqcv
contains

   pure function bs_a(k,beta,gamma) result(price)
      real(dp), intent(in) :: k,beta,gamma
      real(dp) :: price,m
      if (gamma<=0.0_dp) then
         price=max(exp(beta)-k,0.0_dp)
      else if (k>0.0_dp) then
         m=(beta-log(k))/gamma
         price=exp(beta+0.5_dp*gamma*gamma)*normal_cdf(m+gamma)-k*normal_cdf(m)
      else
         price=exp(beta+0.5_dp*gamma*gamma)
      end if
   end function bs_a

   function covariance_conditional_log_prices(t,d,sigma) result(c)
      real(dp), intent(in) :: t,sigma
      integer, intent(in) :: d
      real(dp), allocatable :: c(:,:)
      real(dp), allocatable :: cov(:,:), cw(:)
      real(dp) :: dt, denom
      integer :: i,j
      allocate(c(d,d),cov(d,d),cw(d))
      dt=t/real(d,dp)
      do j=1,d
         do i=1,d
            cov(i,j)=sigma*sigma*dt*real(min(i,j),dp)
         end do
      end do
      cw=sum(cov,dim=2)/real(d,dp)
      denom=sum(cw)/real(d,dp)
      c=cov-spread(cw,2,d)*spread(cw,1,d)/denom
   end function covariance_conditional_log_prices

   subroutine conditional_average_moments(z,t,d,r,sigma,s0,mean_a,var_a,c_in)
      real(dp), intent(in) :: z,t,r,sigma,s0
      integer, intent(in) :: d
      real(dp), intent(out) :: mean_a,var_a
      real(dp), intent(in), optional :: c_in(:,:)
      real(dp), allocatable :: c(:,:), muhat(:), vec(:), evec(:)
      real(dp) :: dt, sigmas, covi, total
      integer :: i,j
      allocate(c(d,d))
      if (present(c_in)) then
         c=c_in
      else
         c=covariance_conditional_log_prices(t,d,sigma)
      end if
      allocate(muhat(d),vec(d),evec(d))
      dt=t/real(d,dp)
      sigmas=sigma/real(d,dp)*sqrt(dt*real(d*(d+1)*(2*d+1),dp)/6.0_dp)
      do i=1,d
         covi=sigma*sigma*dt/real(d,dp)*(real(i*(i+1),dp)/2.0_dp+real((d-i)*i,dp))
         muhat(i)=log(s0)+(r-0.5_dp*sigma*sigma)*dt*real(i,dp)+covi*z/sigmas
         vec(i)=muhat(i)+0.5_dp*c(i,i)
      end do
      evec=exp(vec)
      total=sum(evec)
      mean_a=total/real(d,dp)
      var_a=0.0_dp
      do i=1,d
         do j=1,d
            var_a=var_a+evec(i)*(exp(vec(j)+c(i,j))-evec(j))
         end do
      end do
      var_a=max(0.0_dp,var_a/real(d*d,dp))
   end subroutine conditional_average_moments

   pure function eval_ecv_a(t,d,k,r,sigma,s0) result(price)
      real(dp), intent(in) :: t,k,r,sigma,s0
      integer, intent(in) :: d
      real(dp) :: price,dt,mus,sigmas,zcut,i1,i2,covi,muls,sigma2ls
      integer :: i
      dt=t/real(d,dp)
      mus=log(s0)+(r-0.5_dp*sigma*sigma)*dt*real(d+1,dp)/2.0_dp
      sigmas=sigma/real(d,dp)*sqrt(dt*real(d*(d+1)*(2*d+1),dp)/6.0_dp)
      zcut=(log(k)-mus)/sigmas
      i2=k*normal_cdf(-zcut)
      i1=0.0_dp
      do i=1,d
         muls=log(s0)+(r-0.5_dp*sigma*sigma)*dt*real(i,dp)
         sigma2ls=sigma*sigma*dt*real(i,dp)
         covi=sigma*sigma*dt/real(d,dp)*(real(i*(i+1),dp)/2.0_dp+real((d-i)*i,dp))
         i1=i1+exp(muls+0.5_dp*sigma2ls)*normal_cdf(-zcut+covi/sigmas)
      end do
      price=exp(-r*t)*(i1/real(d,dp)-i2)
   end function eval_ecv_a

   function asian_call_app_lord(t,d,k,r,sigma,s0,full) result(price)
      real(dp), intent(in) :: t,k,r,sigma,s0
      integer, intent(in) :: d
      logical, intent(in), optional :: full
      real(dp) :: price,dt,mus,sigmas,zcut,lower
      logical :: include_all
      real(dp), allocatable :: c(:,:)
      include_all=.true.
      if (present(full)) include_all=full
      dt=t/real(d,dp)
      mus=log(s0)+(r-0.5_dp*sigma*sigma)*dt*real(d+1,dp)/2.0_dp
      sigmas=sigma/real(d,dp)*sqrt(dt*real(d*(d+1)*(2*d+1),dp)/6.0_dp)
      zcut=(log(k)-mus)/sigmas
      allocate(c(d,d))
      c=covariance_conditional_log_prices(t,d,sigma)
      lower=min(-12.0_dp,zcut-8.0_dp)
      price=exp(-r*t)*adaptive_integral(integrand,lower,zcut,1.0e-9_dp,24)
      if (include_all) price=price+eval_ecv_a(t,d,k,r,sigma,s0)
   contains
      function integrand(z) result(q)
         real(dp), intent(in) :: z
         real(dp) :: q,mean_a,var_a,me,gamma2,gamma,beta,shift
         call conditional_average_moments(z,t,d,r,sigma,s0,mean_a,var_a,c)
         shift=exp(mus+sigmas*z)
         me=mean_a-shift
         if (me<=tiny(1.0_dp) .or. var_a<0.0_dp) then
            q=0.0_dp
            return
         end if
         gamma2=log(1.0_dp+var_a/(me*me))
         gamma=sqrt(max(0.0_dp,gamma2))
         beta=log(me)-0.5_dp*gamma2
         q=bs_a(k-shift,beta,gamma)*normal_pdf(z)
      end function integrand
   end function asian_call_app_lord

   pure function eval_ecv(t,d,k,r,sigma,s0) result(res)
      real(dp), intent(in) :: t,k,r,sigma,s0
      integer, intent(in) :: d
      type(moments_result) :: res
      real(dp) :: dt,varx,mus,sigmas,zcut,sum1,hk,hdk,a
      integer :: i,j
      dt=t/real(d,dp)
      varx=real(d*(d+1)*(2*d+1),dp)/6.0_dp
      mus=log(s0)+(r-0.5_dp*sigma*sigma)*dt*real(d+1,dp)/2.0_dp
      sigmas=sigma/real(d,dp)*sqrt(dt*varx)
      zcut=(log(k)-mus)/sigmas
      sum1=0.0_dp; hk=0.0_dp; hdk=0.0_dp
      do i=1,d
         a=sigma*sqrt(dt)*sum([(real(d-j+1,dp)/sqrt(varx),j=1,i)])
         sum1=sum1+exp(r*real(i,dp)*dt)*normal_cdf(-zcut+a)
         hk=hk+exp(a*zcut+r*real(i,dp)*dt-0.5_dp*a*a)
         hdk=hdk+a*exp(a*zcut+r*real(i,dp)*dt-0.5_dp*a*a)
      end do
      hk=s0*hk/real(d,dp)
      hdk=s0*hdk/real(d,dp)
      res%price=exp(-r*t)*(s0*sum1/real(d,dp)-k*normal_cdf(-zcut))
      res%delta=exp(-r*t)*(sum1/real(d,dp)+(hk-k)*normal_pdf(zcut)/(s0*sigmas))
      res%gamma=exp(-r*t)*(2.0_dp*hk*sigmas-hdk-(hk-k)*(sigmas-zcut))* &
         normal_pdf(zcut)/(s0*sigmas)**2
   end function eval_ecv

   function find_bcv(t,d,k,r,sigma,s0,status) result(root)
      real(dp), intent(in) :: t,k,r,sigma,s0
      integer, intent(in) :: d
      integer, intent(out), optional :: status
      real(dp) :: root,dt,mus,sigmas,left,right,fl,fr
      integer :: istat,expand
      dt=t/real(d,dp)
      mus=log(s0)+(r-0.5_dp*sigma*sigma)*dt*real(d+1,dp)/2.0_dp
      sigmas=sigma/real(d,dp)*sqrt(dt*real(d*(d+1)*(2*d+1),dp)/6.0_dp)
      left=-10.0_dp; right=10.0_dp
      fl=f(left); fr=f(right)
      do expand=1,12
         if (fl*fr<=0.0_dp) exit
         left=left-5.0_dp; right=right+5.0_dp
         fl=f(left); fr=f(right)
      end do
      root=bisection_root(f,left,right,1.0e-12_dp,250,istat)
      if (present(status)) status=istat
   contains
      function f(z) result(value)
         real(dp), intent(in) :: z
         real(dp) :: value,x,muls,sigma2ls,covi
         integer :: i
         x=mus+sigmas*z
         value=0.0_dp
         do i=1,d
            muls=log(s0)+(r-0.5_dp*sigma*sigma)*dt*real(i,dp)
            sigma2ls=sigma*sigma*dt*real(i,dp)
            covi=sigma*sigma*dt/real(d,dp)*(real(i*(i+1),dp)/2.0_dp+real((d-i)*i,dp))
            value=value+exp(muls+covi/(sigmas*sigmas)*(x-mus)+ &
               0.5_dp*(sigma2ls-covi*covi/(sigmas*sigmas)))
         end do
         value=value/real(d,dp)-k
      end function f
   end function find_bcv

   function eval_lb(t,d,k,r,sigma,s0,full) result(res)
      real(dp), intent(in) :: t,k,r,sigma,s0
      integer, intent(in) :: d
      logical, intent(in), optional :: full
      type(moments_result) :: res,ecv
      real(dp) :: dt,varx,mus,sigmas,zcut,bcv,sum1,hk,hdk,a
      real(dp) :: k0,h0bcv,hpbcv,bcv0,c1,c2
      integer :: i,j
      logical :: include_all
      include_all=.false.; if (present(full)) include_all=full
      dt=t/real(d,dp); varx=real(d*(d+1)*(2*d+1),dp)/6.0_dp
      mus=log(s0)+(r-0.5_dp*sigma*sigma)*dt*real(d+1,dp)/2.0_dp
      sigmas=sigma/real(d,dp)*sqrt(dt*varx)
      zcut=(log(k)-mus)/sigmas
      bcv=find_bcv(t,d,k,r,sigma,s0)
      sum1=0.0_dp; hk=0.0_dp; hdk=0.0_dp; h0bcv=0.0_dp; hpbcv=0.0_dp
      do i=1,d
         a=sigma*sqrt(dt)*sum([(real(d-j+1,dp)/sqrt(varx),j=1,i)])
         sum1=sum1+exp(r*real(i,dp)*dt)*(normal_cdf(zcut-a)-normal_cdf(bcv-a))
         hk=hk+exp(a*zcut+r*real(i,dp)*dt-0.5_dp*a*a)
         hdk=hdk+a*exp(a*zcut+r*real(i,dp)*dt-0.5_dp*a*a)
         h0bcv=h0bcv+exp(a*bcv+r*real(i,dp)*dt-0.5_dp*a*a)
         hpbcv=hpbcv+a*exp(a*bcv+r*real(i,dp)*dt-0.5_dp*a*a)
      end do
      hk=s0*hk/real(d,dp); hdk=s0*hdk/real(d,dp)
      h0bcv=h0bcv/real(d,dp); hpbcv=s0*hpbcv/real(d,dp)
      res%price=exp(-r*t)*(s0*sum1/real(d,dp)-k*(normal_cdf(zcut)-normal_cdf(bcv)))
      bcv0=-h0bcv/hpbcv
      res%delta=exp(-r*t)*(sum1/real(d,dp)-(hk-k)*normal_pdf(zcut)/(s0*sigmas))
      k0=-1.0_dp/(s0*sigmas)
      c1=(hk*normal_pdf(zcut)*k0-k*normal_pdf(bcv)*bcv0)/s0
      c2=normal_pdf(zcut)/(s0*sigmas)**2*(hdk-hk*sigmas+(hk-k)*(sigmas-zcut))
      res%gamma=exp(-r*t)*(c1+c2)
      if (include_all) then
         ecv=eval_ecv(t,d,k,r,sigma,s0)
         res%price=res%price+ecv%price
         res%delta=res%delta+ecv%delta
         res%gamma=res%gamma+ecv%gamma
      end if
   end function eval_lb

   function eval_equad(t,d,k,r,sigma,s0) result(value)
      real(dp), intent(in) :: t,k,r,sigma,s0
      integer, intent(in) :: d
      real(dp) :: value
      type(moments_result) :: temp
      temp=eval_eqcv(t,d,k,r,sigma,s0)
      value=temp%price
   end function eval_equad

   function eval_eqcv(t,d,k,r,sigma,s0) result(res)
      real(dp), intent(in) :: t,k,r,sigma,s0
      integer, intent(in) :: d
      type(moments_result) :: res
      real(dp), allocatable :: avec(:), gamma_v(:), gamma0(:), gamma00(:), mu(:), vr(:), h(:)
      real(dp), allocatable :: varcov(:,:), mmat(:,:), sigmas2(:,:), gmat(:,:), g0mat(:,:), g00mat(:,:)
      real(dp) :: dt,varx,mus,sigmas,zcut,bcv,eta,term1,term2,term3
      real(dp) :: k0,k00,h0bcv,hpbcv,hppbcv,bcv0,bcv00,eta0,eta00
      real(dp) :: delta_v,gamma_vv,mat1,mat2,vec1,vec2
      integer :: i,j
      dt=t/real(d,dp); varx=real(d*(d+1)*(2*d+1),dp)/6.0_dp
      mus=log(s0)+(r-0.5_dp*sigma*sigma)*dt*real(d+1,dp)/2.0_dp
      sigmas=sigma/real(d,dp)*sqrt(dt*varx); zcut=(log(k)-mus)/sigmas
      bcv=find_bcv(t,d,k,r,sigma,s0)
      allocate(avec(d),gamma_v(d),gamma0(d),gamma00(d),mu(d),vr(d),h(d))
      allocate(varcov(d,d),mmat(d,d),sigmas2(d,d),gmat(d,d),g0mat(d,d),g00mat(d,d))
      do i=1,d
         avec(i)=sigma*sqrt(dt)*sum([(real(d-j+1,dp)/sqrt(varx),j=1,i)])
         gamma_v(i)=exp(0.5_dp*avec(i)**2)*(normal_cdf(zcut-avec(i))-normal_cdf(bcv-avec(i)))/real(d,dp)
         mu(i)=(r-0.5_dp*sigma*sigma)*real(i,dp)*dt
      end do
      eta=normal_cdf(zcut)-normal_cdf(bcv)
      varcov=matmul(lower_ones(d),matmul(identity_matrix(d)-outer(v_vector(d,varx),v_vector(d,varx)), &
         transpose(lower_ones(d))))
      do i=1,d
         vr(i)=varcov(i,i)
         h(i)=exp(mu(i)+0.5_dp*sigma*sigma*dt*vr(i))
      end do
      term1=0.0_dp
      do i=1,d
         do j=1,d
            sigmas2(i,j)=sigma*sigma*dt*(vr(i)+vr(j)+2.0_dp*varcov(i,j))
            mmat(i,j)=exp(mu(i)+mu(j)+0.5_dp*sigmas2(i,j))
            term1=term1+s0*s0*gamma_v(i)*gamma_v(j)*mmat(i,j)
         end do
      end do
      term2=-2.0_dp*s0*sum(gamma_v*h)*k*eta
      term3=(k*eta)**2
      res%price=term1+term2+term3

      k0=-1.0_dp/(s0*sigmas); k00=1.0_dp/(s0*s0*sigmas)
      h0bcv=0.0_dp; hpbcv=0.0_dp; hppbcv=0.0_dp
      do i=1,d
         h0bcv=h0bcv+exp(avec(i)*bcv+r*real(i,dp)*dt-0.5_dp*avec(i)**2)
         hpbcv=hpbcv+avec(i)*exp(avec(i)*bcv+r*real(i,dp)*dt-0.5_dp*avec(i)**2)
         hppbcv=hppbcv+avec(i)**2*exp(avec(i)*bcv+r*real(i,dp)*dt-0.5_dp*avec(i)**2)
      end do
      h0bcv=h0bcv/real(d,dp); hpbcv=s0*hpbcv/real(d,dp); hppbcv=s0*hppbcv/real(d,dp)
      bcv0=-h0bcv/hpbcv
      bcv00=-bcv0*(2.0_dp/s0+hppbcv*bcv0/hpbcv)
      do i=1,d
         gamma0(i)=exp(0.5_dp*avec(i)**2)*(normal_pdf(zcut-avec(i))*k0- &
            normal_pdf(bcv-avec(i))*bcv0)/real(d,dp)
         gamma00(i)=exp(0.5_dp*avec(i)**2)*(normal_pdf(zcut-avec(i))* &
            (k00-k0*k0*(zcut-avec(i)))-normal_pdf(bcv-avec(i))* &
            (bcv00-bcv0*bcv0*(bcv-avec(i))))/real(d,dp)
      end do
      eta0=normal_pdf(zcut)*k0-normal_pdf(bcv)*bcv0
      eta00=normal_pdf(zcut)*(k00-k0*k0*zcut)-normal_pdf(bcv)*(bcv00-bcv0*bcv0*bcv)
      gmat=spread(gamma_v,2,d); g0mat=spread(gamma0,2,d); g00mat=spread(gamma00,2,d)
      delta_v=0.0_dp; gamma_vv=0.0_dp
      do i=1,d
         do j=1,d
            mat1=2.0_dp*s0*gamma_v(i)*gamma_v(j)+s0*s0*(gamma0(i)*gamma_v(j)+gamma_v(i)*gamma0(j))
            mat2=2.0_dp*gamma_v(i)*gamma_v(j)+4.0_dp*s0*(gamma0(i)*gamma_v(j)+gamma_v(i)*gamma0(j))+ &
               s0*s0*(gamma00(i)*gamma_v(j)+2.0_dp*gamma0(i)*gamma0(j)+gamma_v(i)*gamma00(j))
            delta_v=delta_v+mat1*mmat(i,j)
            gamma_vv=gamma_vv+mat2*mmat(i,j)
         end do
      end do
      vec1=sum((gamma_v*eta+s0*(gamma0*eta+gamma_v*eta0))*h)
      vec2=sum((2.0_dp*(gamma0*eta+gamma_v*eta0)+s0*(gamma00*eta+2.0_dp*gamma0*eta0+gamma_v*eta00))*h)
      res%delta=delta_v-2.0_dp*k*vec1+2.0_dp*k*k*eta*eta0
      res%gamma=gamma_vv-2.0_dp*k*vec2+2.0_dp*k*k*(eta0*eta0+eta*eta00)
   contains
      pure function v_vector(n,varx_local) result(v)
         integer, intent(in) :: n
         real(dp), intent(in) :: varx_local
         real(dp) :: v(n)
         integer :: ii
         do ii=1,n
            v(ii)=real(n-ii+1,dp)/sqrt(varx_local)
         end do
      end function v_vector
      pure function outer(x,y) result(a)
         real(dp), intent(in) :: x(:),y(:)
         real(dp) :: a(size(x),size(y))
         a=spread(x,2,size(y))*spread(y,1,size(x))
      end function outer
   end function eval_eqcv
end module optionpricing_asian_analytic
