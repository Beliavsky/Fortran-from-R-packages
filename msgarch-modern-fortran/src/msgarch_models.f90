! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_models
   use msgarch_kinds, only : dp
   use msgarch_types, only : regime_spec, msgarch_spec, regime_state
   use msgarch_distributions, only : distribution_valid, distribution_moments, innovation_logpdf
   implicit none
   private
   public :: create_spec, model_valid, regime_valid, spec_valid
   public :: initialize_state, update_state, conditional_log_density
   public :: parameter_count, pack_parameters, unpack_parameters
   public :: default_parameters, parameter_bounds, stationary_distribution
   public :: unconditional_variances, variance_target_intercept, normalize_transition, extract_regime
   public :: homogeneous_regimes, sort_parameters_by_variance
contains
   pure function lower_string(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i,code
      lower=text
      do i=1,len(text)
         code=iachar(text(i:i));if(code>=65.and.code<=90)lower(i:i)=achar(code+32)
      end do
   end function lower_string

   function create_spec(models,distributions,is_mixture) result(spec)
      character(len=*), intent(in) :: models(:),distributions(:)
      logical, intent(in), optional :: is_mixture
      type(msgarch_spec) :: spec
      integer :: i
      if(size(models)/=size(distributions).or.size(models)<1)error stop 'create_spec: incompatible model/distribution arrays'
      spec%k=size(models);spec%is_mixture=.false.;if(present(is_mixture))spec%is_mixture=is_mixture
      allocate(spec%regime(spec%k),spec%transition(spec%k,spec%k))
      do i=1,spec%k
         spec%regime(i)%model=trim(models(i));spec%regime(i)%distribution=trim(distributions(i))
         call set_default_regime(spec%regime(i))
      end do
      if(spec%k==1)then
         spec%transition=1.0_dp
      else if(spec%is_mixture)then
         spec%transition=1.0_dp/real(spec%k,dp)
      else
         spec%transition=0.1_dp/real(spec%k-1,dp)
         do i=1,spec%k;spec%transition(i,i)=0.9_dp;end do
      end if
   end function create_spec

   subroutine set_default_regime(regime)
      type(regime_spec), intent(inout) :: regime
      character(len=:),allocatable::model,dist
      model=trim(adjustl(lower_string(regime%model)));dist=trim(adjustl(lower_string(regime%distribution)))
      select case(model)
      case('sarch');regime%omega=0.1_dp;regime%alpha=0.1_dp;regime%gamma=0.0_dp;regime%beta=0.0_dp
      case('sgarch');regime%omega=0.1_dp;regime%alpha=0.1_dp;regime%gamma=0.0_dp;regime%beta=0.8_dp
      case('gjrgarch');regime%omega=0.1_dp;regime%alpha=0.05_dp;regime%gamma=0.1_dp;regime%beta=0.8_dp
      case('egarch');regime%omega=0.0_dp;regime%alpha=0.2_dp;regime%gamma=-0.1_dp;regime%beta=0.8_dp
      case('tgarch');regime%omega=0.125_dp;regime%alpha=0.05_dp;regime%gamma=0.01_dp;regime%beta=0.8_dp
      end select
      if(index(dist,'std')>0)regime%shape=8.0_dp
      if(index(dist,'ged')>0)regime%shape=2.0_dp
      regime%skew=1.0_dp
   end subroutine set_default_regime

   pure function model_valid(model) result(ok)
      character(len=*),intent(in)::model
      logical::ok
      select case(trim(adjustl(lower_string(model))))
      case('sarch','sgarch','egarch','gjrgarch','tgarch');ok=.true.
      case default;ok=.false.
      end select
   end function model_valid

   function regime_valid(regime) result(ok)
      type(regime_spec),intent(in)::regime
      logical::ok
      real(dp)::eabs,ezineg,ez2ineg,ineq
      character(len=:),allocatable::model
      model=trim(adjustl(lower_string(regime%model)))
      if(.not.model_valid(model).or..not.distribution_valid(regime%distribution,regime%shape,regime%skew))then;ok=.false.;return;end if
      call distribution_moments(regime%distribution,regime%shape,regime%skew,eabs,ezineg,ez2ineg)
      select case(model)
      case('sarch')
         ok=regime%omega>0.0_dp.and.regime%alpha>0.0_dp.and.regime%alpha<0.999999_dp
      case('sgarch')
         ok=regime%omega>0.0_dp.and.regime%alpha>0.0_dp.and.regime%beta>=0.0_dp.and. &
            regime%alpha+regime%beta<0.999999_dp
      case('gjrgarch')
         ineq=regime%alpha+ez2ineg*regime%gamma+regime%beta
         ok=regime%omega>0.0_dp.and.regime%alpha>0.0_dp.and.regime%gamma>=0.0_dp.and.regime%beta>=0.0_dp.and.ineq<0.999999_dp
      case('egarch')
         ok=abs(regime%beta)<0.999999_dp.and.abs(regime%omega)<50.0_dp.and.abs(regime%alpha)<5.0_dp.and.abs(regime%gamma)<5.0_dp
      case('tgarch')
         ineq=regime%alpha**2+regime%beta**2-2.0_dp*(regime%alpha+regime%gamma)*regime%beta*ezineg- &
            (regime%alpha**2-regime%gamma**2)*ez2ineg
         ok=regime%omega>0.0_dp.and.regime%alpha>0.0_dp.and.regime%gamma>0.0_dp.and.regime%beta>=0.0_dp.and.ineq<0.999999_dp
      case default;ok=.false.
      end select
   end function regime_valid

   function spec_valid(spec) result(ok)
      type(msgarch_spec),intent(in)::spec
      logical::ok
      integer::i
      ok=spec%k>=1.and.allocated(spec%regime).and.allocated(spec%transition)
      if(.not.ok)return
      if(size(spec%regime)/=spec%k.or.any(shape(spec%transition)/=[spec%k,spec%k]))then;ok=.false.;return;end if
      do i=1,spec%k;if(.not.regime_valid(spec%regime(i)))then;ok=.false.;return;end if;end do
      if(any(spec%transition<=0.0_dp).or.any(abs(sum(spec%transition,dim=2)-1.0_dp)>1.0e-8_dp))ok=.false.
      if(spec%is_mixture)then
         do i=2,spec%k;if(maxval(abs(spec%transition(i,:)-spec%transition(1,:)))>1.0e-10_dp)ok=.false.;end do
      end if
   end function spec_valid

   subroutine initialize_state(regime,state)
      type(regime_spec),intent(in)::regime
      type(regime_state),intent(out)::state
      real(dp)::eabs,ezineg,ez2ineg,denom
      character(len=:),allocatable::model
      call distribution_moments(regime%distribution,regime%shape,regime%skew,eabs,ezineg,ez2ineg)
      model=trim(adjustl(lower_string(regime%model)))
      select case(model)
      case('sarch');state%h=regime%omega/(1.0_dp-regime%alpha)
      case('sgarch');state%h=regime%omega/(1.0_dp-regime%alpha-regime%beta)
      case('gjrgarch');state%h=regime%omega/(1.0_dp-regime%alpha-ez2ineg*regime%gamma-regime%beta)
      case('egarch');state%lnh=regime%omega/(1.0_dp-regime%beta);state%h=exp(state%lnh)
      case('tgarch')
         denom=1.0_dp+(regime%alpha+regime%gamma)*ezineg-regime%beta
         state%fh=regime%omega/denom;state%h=state%fh**2
      end select
      state%h=max(state%h,1.0e-12_dp);state%lnh=log(state%h);if(model/='tgarch')state%fh=sqrt(state%h)
   end subroutine initialize_state

   subroutine update_state(regime,state,previous_y)
      type(regime_spec),intent(in)::regime
      type(regime_state),intent(inout)::state
      real(dp),intent(in)::previous_y
      real(dp)::eabs,ezineg,ez2ineg,z
      character(len=:),allocatable::model
      model=trim(adjustl(lower_string(regime%model)))
      select case(model)
      case('sarch')
         state%h=regime%omega+regime%alpha*previous_y**2
      case('sgarch')
         state%h=regime%omega+regime%alpha*previous_y**2+regime%beta*state%h
      case('gjrgarch')
         state%h=regime%omega+regime%alpha*previous_y**2+regime%beta*state%h
         if(previous_y<0.0_dp)state%h=state%h+regime%gamma*previous_y**2
      case('egarch')
         call distribution_moments(regime%distribution,regime%shape,regime%skew,eabs,ezineg,ez2ineg)
         z=previous_y/sqrt(max(state%h,1.0e-12_dp))
         state%lnh=regime%omega+regime%alpha*(abs(z)-eabs)+regime%gamma*z+regime%beta*state%lnh
         state%h=exp(max(-50.0_dp,min(50.0_dp,state%lnh)))
      case('tgarch')
         state%fh=regime%omega+regime%beta*state%fh
         if(previous_y>=0.0_dp)then;state%fh=state%fh+regime%alpha*previous_y
         else;state%fh=state%fh-regime%gamma*previous_y;end if
         state%fh=max(state%fh,1.0e-8_dp);state%h=state%fh**2
      end select
      state%h=max(state%h,1.0e-12_dp);state%lnh=log(state%h)
   end subroutine update_state

   pure function conditional_log_density(regime,state,y) result(value)
      type(regime_spec),intent(in)::regime
      type(regime_state),intent(in)::state
      real(dp),intent(in)::y
      real(dp)::value,z
      z=y/sqrt(state%h)
      value=innovation_logpdf(z,regime%distribution,regime%shape,regime%skew)-0.5_dp*log(state%h)
   end function conditional_log_density

   pure function model_parameter_count(model) result(n)
      character(len=*),intent(in)::model
      integer::n
      select case(trim(adjustl(lower_string(model))))
      case('sarch');n=2
      case('sgarch');n=3
      case('egarch','gjrgarch','tgarch');n=4
      case default;n=0
      end select
   end function model_parameter_count

   pure function distribution_parameter_count(distribution) result(n)
      character(len=*),intent(in)::distribution
      integer::n
      select case(trim(adjustl(lower_string(distribution))))
      case('norm');n=0
      case('std','ged','snorm');n=1
      case('sstd','sged');n=2
      case default;n=0
      end select
   end function distribution_parameter_count

   function parameter_count(spec) result(n)
      type(msgarch_spec),intent(in)::spec
      integer::n,i
      n=0
      do i=1,spec%k;n=n+model_parameter_count(spec%regime(i)%model)+distribution_parameter_count(spec%regime(i)%distribution);end do
      if(spec%k>1)then
         if(spec%is_mixture)then;n=n+spec%k-1;else;n=n+spec%k*(spec%k-1);end if
      end if
   end function parameter_count

   function pack_parameters(spec) result(theta)
      type(msgarch_spec),intent(in)::spec
      real(dp),allocatable::theta(:)
      integer::i,j,pos
      allocate(theta(parameter_count(spec)));pos=1
      do i=1,spec%k
         theta(pos)=spec%regime(i)%omega;theta(pos+1)=spec%regime(i)%alpha;pos=pos+2
         select case(trim(adjustl(lower_string(spec%regime(i)%model))))
         case('sgarch');theta(pos)=spec%regime(i)%beta;pos=pos+1
         case('egarch','gjrgarch','tgarch');theta(pos)=spec%regime(i)%gamma;theta(pos+1)=spec%regime(i)%beta;pos=pos+2
         end select
         select case(trim(adjustl(lower_string(spec%regime(i)%distribution))))
         case('std','ged');theta(pos)=spec%regime(i)%shape;pos=pos+1
         case('snorm');theta(pos)=spec%regime(i)%skew;pos=pos+1
         case('sstd','sged');theta(pos)=spec%regime(i)%shape;theta(pos+1)=spec%regime(i)%skew;pos=pos+2
         end select
      end do
      if(spec%k>1)then
         if(spec%is_mixture)then
            theta(pos:pos+spec%k-2)=spec%transition(1,1:spec%k-1)
         else
            do i=1,spec%k
               do j=1,spec%k-1;theta(pos)=spec%transition(i,j);pos=pos+1;end do
            end do
         end if
      end if
   end function pack_parameters

   subroutine unpack_parameters(template,theta,spec,valid)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::theta(:)
      type(msgarch_spec),intent(out)::spec
      logical,intent(out),optional::valid
      integer::i,j,pos
      real(dp)::last
      spec=template;pos=1
      if (size(theta) /= parameter_count(template)) then
         if (present(valid)) valid = .false.
         return
      else
         if (present(valid)) valid = .true.
      end if
      do i=1,spec%k
         spec%regime(i)%omega=theta(pos);spec%regime(i)%alpha=theta(pos+1);pos=pos+2
         select case(trim(adjustl(lower_string(spec%regime(i)%model))))
         case('sgarch');spec%regime(i)%beta=theta(pos);pos=pos+1
         case('egarch','gjrgarch','tgarch');spec%regime(i)%gamma=theta(pos);spec%regime(i)%beta=theta(pos+1);pos=pos+2
         end select
         select case(trim(adjustl(lower_string(spec%regime(i)%distribution))))
         case('std','ged');spec%regime(i)%shape=theta(pos);pos=pos+1
         case('snorm');spec%regime(i)%skew=theta(pos);pos=pos+1
         case('sstd','sged');spec%regime(i)%shape=theta(pos);spec%regime(i)%skew=theta(pos+1);pos=pos+2
         end select
      end do
      if(spec%k>1)then
         if(spec%is_mixture)then
            spec%transition(:,1:spec%k-1)=spread(theta(pos:pos+spec%k-2),1,spec%k)
            last=1.0_dp-sum(theta(pos:pos+spec%k-2));spec%transition(:,spec%k)=last
         else
            do i=1,spec%k
               do j=1,spec%k-1;spec%transition(i,j)=theta(pos);pos=pos+1;end do
               spec%transition(i,spec%k)=1.0_dp-sum(spec%transition(i,1:spec%k-1))
            end do
         end if
      end if
      if(present(valid))valid=spec_valid(spec)
   end subroutine unpack_parameters

   function default_parameters(spec) result(theta)
      type(msgarch_spec),intent(in)::spec
      real(dp),allocatable::theta(:)
      theta=pack_parameters(spec)
   end function default_parameters

   subroutine parameter_bounds(spec,lower,upper)
      type(msgarch_spec),intent(in)::spec
      real(dp),allocatable,intent(out)::lower(:),upper(:)
      integer::i,pos,n
      character(len=:),allocatable::model,dist
      n=parameter_count(spec);allocate(lower(n),upper(n));pos=1
      do i=1,spec%k
         model=trim(adjustl(lower_string(spec%regime(i)%model)));dist=trim(adjustl(lower_string(spec%regime(i)%distribution)))
         select case(model)
         case('egarch');lower(pos)=-50.0_dp;upper(pos)=50.0_dp;lower(pos+1)=-5.0_dp;upper(pos+1)=5.0_dp
         case default;lower(pos)=1.0e-7_dp;upper(pos)=100.0_dp;lower(pos+1)=1.0e-6_dp;upper(pos+1)=0.9999_dp
         end select
         pos=pos+2
         select case(model)
         case('sgarch');lower(pos)=0.0_dp;upper(pos)=0.9999_dp;pos=pos+1
         case('egarch');lower(pos)=-5.0_dp;upper(pos)=5.0_dp;lower(pos+1)=-0.9999_dp;upper(pos+1)=0.9999_dp;pos=pos+2
         case('gjrgarch');lower(pos)=1.0e-6_dp;upper(pos)=10.0_dp;lower(pos+1)=0.0_dp;upper(pos+1)=0.9999_dp;pos=pos+2
         case('tgarch');lower(pos)=1.0e-6_dp;upper(pos)=10.0_dp;lower(pos+1)=0.0_dp;upper(pos+1)=10.0_dp;pos=pos+2
         end select
         select case(dist)
         case('std');lower(pos)=2.1_dp;upper(pos)=100.0_dp;pos=pos+1
         case('ged');lower(pos)=0.7_dp;upper(pos)=20.0_dp;pos=pos+1
         case('snorm');lower(pos)=0.01_dp;upper(pos)=100.0_dp;pos=pos+1
         case('sstd');lower(pos)=2.1_dp;upper(pos)=100.0_dp;lower(pos+1)=0.01_dp;upper(pos+1)=100.0_dp;pos=pos+2
         case('sged');lower(pos)=0.7_dp;upper(pos)=20.0_dp;lower(pos+1)=0.01_dp;upper(pos+1)=100.0_dp;pos=pos+2
         end select
      end do
      if(pos<=n)then;lower(pos:n)=1.0e-5_dp;upper(pos:n)=0.99999_dp;end if
   end subroutine parameter_bounds

   subroutine normalize_transition(transition)
      real(dp),intent(inout)::transition(:,:)
      integer::i
      do i=1,size(transition,1)
         transition(i,:)=max(transition(i,:),1.0e-12_dp)
         transition(i,:)=transition(i,:)/sum(transition(i,:))
      end do
   end subroutine normalize_transition

   function stationary_distribution(transition) result(probability)
      real(dp),intent(in)::transition(:,:)
      real(dp),allocatable::probability(:),next(:)
      integer::i,n
      n=size(transition,1);allocate(probability(n),next(n));probability=1.0_dp/real(n,dp)
      do i=1,10000
         next=matmul(probability,transition)
         if(maxval(abs(next-probability))<1.0e-14_dp)exit
         probability=next
      end do
      probability=max(next,0.0_dp);probability=probability/sum(probability)
   end function stationary_distribution


   function extract_regime(spec,index) result(single)
      type(msgarch_spec), intent(in) :: spec
      integer, intent(in) :: index
      type(msgarch_spec) :: single
      if (index < 1 .or. index > spec%k) error stop 'extract_regime: index out of range'
      single%k = 1
      single%is_mixture = .false.
      allocate(single%regime(1),single%transition(1,1))
      single%regime(1) = spec%regime(index)
      single%transition(1,1) = 1.0_dp
   end function extract_regime

   pure function homogeneous_regimes(spec) result(homogeneous)
      type(msgarch_spec),intent(in)::spec
      logical::homogeneous
      integer::i
      homogeneous=spec%k>=1
      do i=2,spec%k
         if(trim(adjustl(lower_string(spec%regime(i)%model)))/=trim(adjustl(lower_string(spec%regime(1)%model))))homogeneous=.false.
         if(trim(adjustl(lower_string(spec%regime(i)%distribution)))/= &
            trim(adjustl(lower_string(spec%regime(1)%distribution))))homogeneous=.false.
      end do
   end function homogeneous_regimes

   function sort_parameters_by_variance(template,theta) result(sorted_theta)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::theta(:)
      real(dp),allocatable::sorted_theta(:)
      type(msgarch_spec)::spec,sorted_spec
      real(dp),allocatable::variance(:)
      integer,allocatable::order(:)
      integer::i,j,key
      logical::valid
      sorted_theta=theta
      if(template%k<=1.or..not.homogeneous_regimes(template))return
      call unpack_parameters(template,theta,spec,valid);if(.not.valid)return
      variance=unconditional_variances(spec);allocate(order(spec%k));order=[(i,i=1,spec%k)]
      do i=2,spec%k
         key=order(i);j=i-1
         do while(j>=1)
            if(variance(order(j))<=variance(key))exit
            order(j+1)=order(j);j=j-1
         end do
         order(j+1)=key
      end do
      sorted_spec=spec
      do i=1,spec%k
         sorted_spec%regime(i)=spec%regime(order(i))
         do j=1,spec%k
            sorted_spec%transition(i,j)=spec%transition(order(i),order(j))
         end do
      end do
      sorted_theta=pack_parameters(sorted_spec)
   end function sort_parameters_by_variance

   function variance_target_intercept(target_variance,regime) result(intercept)
      real(dp),intent(in)::target_variance
      type(regime_spec),intent(in)::regime
      real(dp)::intercept
      select case(trim(adjustl(lower_string(regime%model))))
      case('sarch')
         intercept=target_variance*(1.0_dp-regime%alpha)
      case('sgarch')
         intercept=target_variance*(1.0_dp-regime%alpha-regime%beta)
      case('gjrgarch')
         intercept=target_variance*(1.0_dp-regime%alpha-0.5_dp*regime%gamma-regime%beta)
      case('egarch')
         intercept=log(max(target_variance,tiny(1.0_dp)))*(1.0_dp-regime%beta)
      case('tgarch')
         intercept=target_variance*(1.0_dp+0.5_dp*(regime%alpha+regime%gamma)-regime%beta)
      case default
         intercept=0.0_dp
      end select
   end function variance_target_intercept

   function unconditional_variances(spec) result(variance)
      type(msgarch_spec),intent(in)::spec
      real(dp),allocatable::variance(:)
      type(regime_state)::state
      integer::i
      allocate(variance(spec%k))
      do i=1,spec%k;call initialize_state(spec%regime(i),state);variance(i)=state%h;end do
   end function unconditional_variances
end module msgarch_models
