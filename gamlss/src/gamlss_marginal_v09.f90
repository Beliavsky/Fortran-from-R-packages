! Upstream-compatible marginal prediction for random() terms.
! Mirrors gamlss::getMarginal() methods: integrate, qfunction, random, none.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_marginal_v09
   use gamlss_kinds, only : dp
   use gamlss_fit, only : family_npar, map_parameters
   use gamlss_special, only : normal_pdf, normal_quantile
   use gamlss_random, only : random_normal
   use gamlss_v02_numerics, only : adaptive_integral
   use gamlss_random_effects, only : random_intercept_result_t
   implicit none
   private

   integer, parameter, public :: MARGINAL_INTEGRATE = 1
   integer, parameter, public :: MARGINAL_QFUNCTION = 2
   integer, parameter, public :: MARGINAL_RANDOM = 3
   integer, parameter, public :: MARGINAL_NONE = 4

   type, public :: marginal_prediction_result_t
      real(dp), allocatable :: fitted(:)
      real(dp) :: sigma_b = 0.0_dp
      integer :: parameter = 1
      integer :: method = MARGINAL_INTEGRATE
      integer :: integration_points = 0
      integer :: status = 0
   end type marginal_prediction_result_t

   public :: marginal_predict_eta, get_marginal_random_intercept, parameter_link_inverse

   real(dp), save :: ctx_eta = 0.0_dp, ctx_sigma = 1.0_dp
   integer, save :: ctx_family = 0, ctx_parameter = 1

contains

   subroutine marginal_predict_eta(eta_without_random,sigma_b,family,parameter,result,method,n_random)
      real(dp), intent(in) :: eta_without_random(:),sigma_b
      integer, intent(in) :: family,parameter
      type(marginal_prediction_result_t), intent(out) :: result
      integer, intent(in), optional :: method,n_random
      integer :: meth,nr,i,k,np
      real(dp), allocatable :: zq(:)
      real(dp) :: s,bound

      result%status=0
      meth=MARGINAL_INTEGRATE;if(present(method))meth=method
      nr=10000;if(present(n_random))nr=n_random
      np=family_npar(family)
      if(size(eta_without_random)==0.or.parameter<1.or.parameter>np.or.sigma_b<0.0_dp)then
         result%status=1;return
      end if
      if(meth<MARGINAL_INTEGRATE.or.meth>MARGINAL_NONE.or.nr<1)then
         result%status=2;return
      end if
      allocate(result%fitted(size(eta_without_random)))
      result%sigma_b=sigma_b;result%parameter=parameter;result%method=meth

      select case(meth)
      case(MARGINAL_NONE)
         do i=1,size(eta_without_random)
            result%fitted(i)=parameter_link_inverse(family,parameter,eta_without_random(i))
         end do
      case(MARGINAL_QFUNCTION)
         allocate(zq(999))
         do k=1,999
            zq(k)=normal_quantile(real(k,dp)/1000.0_dp)*sigma_b
         end do
         do i=1,size(eta_without_random)
            s=0.0_dp
            do k=1,999
               s=s+parameter_link_inverse(family,parameter,eta_without_random(i)+zq(k))
            end do
            result%fitted(i)=s/999.0_dp
         end do
         result%integration_points=999
      case(MARGINAL_RANDOM)
         do i=1,size(eta_without_random)
            s=0.0_dp
            do k=1,nr
               s=s+parameter_link_inverse(family,parameter,eta_without_random(i)+sigma_b*random_normal())
            end do
            result%fitted(i)=s/real(nr,dp)
         end do
         result%integration_points=nr
      case(MARGINAL_INTEGRATE)
         ! Upstream integrates from -Inf to Inf.  After standardization the
         ! normal tail outside this adaptive finite range is negligible; the
         ! extra sigma_b term also covers exponentially tilted inverse links.
         bound=max(9.0_dp,9.0_dp+sigma_b)
         ctx_sigma=sigma_b;ctx_family=family;ctx_parameter=parameter
         do i=1,size(eta_without_random)
            ctx_eta=eta_without_random(i)
            result%fitted(i)=adaptive_integral(marginal_integrand,-bound,bound,2.0e-10_dp,24)
         end do
      end select
   end subroutine marginal_predict_eta

   subroutine get_marginal_random_intercept(fit,group,result,method,n_random,sigma_b)
      type(random_intercept_result_t), intent(in) :: fit
      integer, intent(in) :: group(:)
      type(marginal_prediction_result_t), intent(out) :: result
      integer, intent(in), optional :: method,n_random
      real(dp), intent(in), optional :: sigma_b
      real(dp), allocatable :: eta(:),rt(:)
      real(dp) :: sb
      integer :: i,gidx,n

      result%status=0;n=size(group)
      if(fit%status/=0.or.n==0.or..not.allocated(fit%levels).or..not.allocated(fit%effects))then
         result%status=10;return
      end if
      allocate(eta(n),rt(n));rt=0.0_dp
      select case(fit%parameter)
      case(1)
         if(.not.allocated(fit%model%mu%eta).or.size(fit%model%mu%eta)/=n)then;result%status=11;return;end if
         eta=fit%model%mu%eta
      case(2)
         if(.not.allocated(fit%model%sigma%eta).or.size(fit%model%sigma%eta)/=n)then;result%status=11;return;end if
         eta=fit%model%sigma%eta
      case(3)
         if(.not.allocated(fit%model%nu%eta).or.size(fit%model%nu%eta)/=n)then;result%status=11;return;end if
         eta=fit%model%nu%eta
      case(4)
         if(.not.allocated(fit%model%tau%eta).or.size(fit%model%tau%eta)/=n)then;result%status=11;return;end if
         eta=fit%model%tau%eta
      case default
         result%status=12;return
      end select
      do i=1,n
         gidx=level_index(group(i),fit%levels)
         if(gidx==0)then;result%status=13;return;end if
         rt(i)=fit%effects(gidx)
      end do
      sb=fit%sigma_b;if(present(sigma_b))sb=sigma_b
      if(sb<0.0_dp)then;result%status=14;return;end if
      call marginal_predict_eta(eta-rt,sb,fit%model%family,fit%parameter,result,method,n_random)
   end subroutine get_marginal_random_intercept

   real(dp) function parameter_link_inverse(family,parameter,eta) result(value)
      integer, intent(in) :: family,parameter
      real(dp), intent(in) :: eta
      real(dp) :: e(4),a,b,c,d
      e=0.0_dp
      if(parameter<1.or.parameter>4)then;value=eta;return;end if
      e(parameter)=eta
      call map_parameters(family,e(1),e(2),e(3),e(4),a,b,c,d)
      select case(parameter)
      case(1);value=a
      case(2);value=b
      case(3);value=c
      case(4);value=d
      end select
   end function parameter_link_inverse

   real(dp) function marginal_integrand(z) result(v)
      real(dp), intent(in) :: z
      v=parameter_link_inverse(ctx_family,ctx_parameter,ctx_eta+ctx_sigma*z)*normal_pdf(z)
   end function marginal_integrand

   integer function level_index(level,levels) result(idx)
      integer, intent(in) :: level,levels(:)
      integer :: j
      idx=0
      do j=1,size(levels)
         if(levels(j)==level)then;idx=j;return;end if
      end do
   end function level_index

end module gamlss_marginal_v09
