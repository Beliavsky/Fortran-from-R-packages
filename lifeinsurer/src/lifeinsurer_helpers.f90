! SPDX-License-Identifier: GPL-2.0-or-later
module lifeinsurer_helpers
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
  use lifeinsurer_kinds, only : dp, lir_success, lir_invalid_argument
  use lifeinsurer_types, only : civil_date, frequency_correction, expense_loadings
  implicit none
  private
  public :: is_single_premium_contract, is_regular_premium_contract
  public :: death_benefit_linear_decreasing, death_benefit_annuity_decreasing
  public :: premium_refund_period_default, correction_payment_frequency
  public :: frequency_charge, round_value, round_vector, fill_na_gaps
  public :: rolling_mean, pad_zero, pad_last, head_zero, get_savings_premium
  public :: age_year_difference, age_exact_rounded, initialize_costs
  public :: costs_base_alpha, costs_scale_alpha

  interface round_value
    module procedure round_value_scalar
    module procedure round_value_vector
  end interface round_value
contains
  pure logical function is_single_premium_contract(premium_period)
    integer, intent(in) :: premium_period
    is_single_premium_contract = premium_period <= 1
  end function

  pure logical function is_regular_premium_contract(premium_period)
    integer, intent(in) :: premium_period
    is_regular_premium_contract = premium_period > 1
  end function

  pure integer function premium_refund_period_default(policy_period, deferral_period)
    integer, intent(in) :: policy_period, deferral_period
    if (deferral_period > 0) then
      premium_refund_period_default = deferral_period
    else
      premium_refund_period_default = policy_period
    end if
  end function

  function death_benefit_linear_decreasing(policy_period, deferral_period, n) result(b)
    integer, intent(in) :: policy_period, deferral_period, n
    real(dp), allocatable :: b(:)
    integer :: k, p
    allocate(b(max(0,n))); b = 0.0_dp
    p = policy_period - deferral_period
    if (p <= 0) return
    do k = 1, min(n,p+1)
      b(k) = real(p-k+1,dp)/real(p,dp)
    end do
  end function

  function death_benefit_annuity_decreasing(interest, policy_period, deferral_period, n) result(b)
    real(dp), intent(in) :: interest
    integer, intent(in) :: policy_period, deferral_period, n
    real(dp), allocatable :: b(:)
    integer :: k, p
    real(dp) :: v
    allocate(b(max(0,n))); b = 0.0_dp
    p = policy_period - deferral_period
    if (p <= 0) return
    if (abs(interest) <= epsilon(1.0_dp)) then
      do k=1,min(n,p+1)
        b(k)=real(p-k+1,dp)/real(p,dp)
      end do
    else
      v=1.0_dp/(1.0_dp+interest)
      do k=1,min(n,p+1)
        b(k)=(v**real(p-k+1,dp)-1.0_dp)/(v**real(p,dp)-1.0_dp)
      end do
    end if
  end function

  pure function correction_payment_frequency(i, m, order) result(c)
    real(dp), intent(in) :: i, order
    integer, intent(in) :: m
    type(frequency_correction) :: c
    real(dp) :: mm, d, im, dm
    c%alpha=1.0_dp; c%beta=0.0_dp
    if (m <= 0) return
    mm=real(m,dp)
    if (order >= 0.0_dp) c%beta=c%beta+(mm-1.0_dp)/(2.0_dp*mm)
    if (order >= 1.0_dp) c%beta=c%beta+(mm*mm-1.0_dp)/(6.0_dp*mm*mm)*i
    if (abs(order-1.5_dp) < 10.0_dp*epsilon(1.0_dp)) &
      c%beta=c%beta+(1.0_dp-mm*mm)/(12.0_dp*mm*mm)*i*i
    if (order >= 2.0_dp) then
      c%beta=c%beta+(1.0_dp-mm*mm)/(24.0_dp*mm*mm)*i*i
      c%alpha=c%alpha+(mm*mm-1.0_dp)/(12.0_dp*mm*mm)*i*i
    end if
    if (order > huge(1.0_dp)/2.0_dp .and. abs(i)>epsilon(1.0_dp)) then
      d=i/(1.0_dp+i); im=mm*((1.0_dp+i)**(1.0_dp/mm)-1.0_dp)
      dm=im/(1.0_dp+im/mm)
      if (abs(dm*im)>tiny(1.0_dp)) then
        c%alpha=d*i/(dm*im); c%beta=(i-im)/(dm*im)
      end if
    end if
  end function

  pure real(dp) function frequency_charge(frequency, monthly, quarterly, semiannually, yearly)
    integer, intent(in) :: frequency
    real(dp), intent(in) :: monthly, quarterly, semiannually, yearly
    real(dp) :: a,b,c,d
    a=monthly; b=quarterly; c=semiannually; d=yearly
    if (a>0.1_dp) a=a/100.0_dp
    if (b>0.1_dp) b=b/100.0_dp
    if (c>0.1_dp) c=c/100.0_dp
    if (d>0.1_dp) d=d/100.0_dp
    select case(frequency)
    case(12); frequency_charge=a
    case(4); frequency_charge=b
    case(2); frequency_charge=c
    case default; frequency_charge=d
    end select
  end function

  elemental real(dp) function round_value_scalar(x, digits)
    real(dp), intent(in) :: x
    integer, intent(in) :: digits
    real(dp) :: s
    s=10.0_dp**real(digits,dp)
    round_value_scalar=anint(x*s)/s
  end function

  function round_value_vector(x, digits) result(y)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: digits
    real(dp), allocatable :: y(:)
    allocate(y(size(x))); y=round_value_scalar(x,digits)
  end function

  subroutine round_vector(x, digits)
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: digits
    x=round_value_scalar(x,digits)
  end subroutine

  function pad_zero(v, n, value, start, value_start) result(out)
    real(dp), intent(in) :: v(:)
    integer, intent(in) :: n
    real(dp), intent(in), optional :: value, value_start
    integer, intent(in), optional :: start
    real(dp), allocatable :: out(:)
    real(dp) :: fill, lead
    integer :: s, k
    fill=0.0_dp; if(present(value)) fill=value
    lead=0.0_dp; if(present(value_start)) lead=value_start
    s=0; if(present(start)) s=max(0,start)
    allocate(out(max(0,n))); out=fill
    if(s>0 .and. n>0) out(1:min(s,n))=lead
    k=min(size(v),max(0,n-s))
    if(k>0) out(s+1:s+k)=v(1:k)
  end function

  function head_zero(v, start, value_start) result(out)
    real(dp), intent(in) :: v(:)
    integer, intent(in), optional :: start
    real(dp), intent(in), optional :: value_start
    real(dp), allocatable :: out(:)
    integer :: s
    real(dp) :: z
    s=0; if(present(start)) s=max(0,start)
    z=0.0_dp; if(present(value_start)) z=value_start
    allocate(out(size(v)+1)); out(1)=z; out(2:)=v
    if(s>0) out(1:min(s,size(out)))=z
  end function

  function pad_last(v, n) result(out)
    real(dp), intent(in) :: v(:)
    integer, intent(in) :: n
    real(dp), allocatable :: out(:)
    integer :: k
    allocate(out(max(0,n)))
    if(n<=0) return
    if(size(v)==0) then; out=0.0_dp; return; end if
    k=min(n,size(v)); out(1:k)=v(1:k)
    if(k<n) out(k+1:n)=v(size(v))
  end function

  function rolling_mean(x) result(y)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: y(:)
    integer :: n
    n=max(0,size(x)-1); allocate(y(n))
    if(n>0) y=0.5_dp*(x(1:n)+x(2:n+1))
  end function

  subroutine fill_na_gaps(x, first_back)
    real(dp), intent(inout) :: x(:)
    logical, intent(in), optional :: first_back
    logical :: back
    integer :: i, first
    back=.false.; if(present(first_back)) back=first_back
    first=0
    do i=1,size(x)
      if(ieee_is_finite(x(i)) .and. .not.ieee_is_nan(x(i))) then
        first=i; exit
      end if
    end do
    if(first==0) return
    if(back .and. first>1) x(1:first-1)=x(first)
    do i=first+1,size(x)
      if(.not.ieee_is_finite(x(i)) .or. ieee_is_nan(x(i))) x(i)=x(i-1)
    end do
  end subroutine fill_na_gaps

  function get_savings_premium(reserves, v, survival_advance, survival_arrears) result(sp)
    real(dp), intent(in) :: reserves(:), v
    real(dp), intent(in), optional :: survival_advance(:), survival_arrears(:)
    real(dp), allocatable :: sp(:)
    integer :: n
    n=size(reserves); allocate(sp(n)); sp=-reserves
    if(n>1) sp(1:n-1)=sp(1:n-1)+v*reserves(2:n)
    if(present(survival_advance)) sp=sp+survival_advance
    if(present(survival_arrears)) sp=sp+v*survival_arrears
  end function

  pure integer function age_year_difference(birth, closing)
    type(civil_date), intent(in) :: birth, closing
    age_year_difference=closing%year-birth%year
  end function

  pure integer function age_exact_rounded(birth, closing)
    type(civil_date), intent(in) :: birth, closing
    real(dp) :: years
    integer :: mdiff
    mdiff=(closing%year-birth%year)*12+closing%month-birth%month
    years=(real(mdiff,dp)+(real(closing%day-birth%day,dp)/30.4375_dp))/12.0_dp
    age_exact_rounded=nint(years)
  end function

  pure function initialize_costs(alpha,zillmer,beta,gamma,gamma_contract, &
      gamma_after_death,unit_cost,unit_cost_policy,security) result(c)
    real(dp), intent(in), optional :: alpha,zillmer,beta,gamma,gamma_contract
    real(dp), intent(in), optional :: gamma_after_death,unit_cost
    real(dp), intent(in), optional :: unit_cost_policy,security
    type(expense_loadings) :: c
    if(present(alpha)) c%alpha=alpha
    if(present(zillmer)) c%zillmer=zillmer
    if(present(beta)) c%beta=beta
    if(present(gamma)) c%gamma=gamma
    if(present(gamma_contract)) c%gamma_contract=gamma_contract
    if(present(gamma_after_death)) c%gamma_after_death=gamma_after_death
    if(present(unit_cost)) c%unit_cost=unit_cost
    if(present(unit_cost_policy)) c%unit_cost_policy=unit_cost_policy
    if(present(security)) c%security=security
  end function

  pure real(dp) function costs_base_alpha(alpha)
    real(dp), intent(in) :: alpha
    costs_base_alpha=alpha
  end function

  pure real(dp) function costs_scale_alpha(alpha, scale)
    real(dp), intent(in) :: alpha,scale
    costs_scale_alpha=alpha*scale
  end function
end module lifeinsurer_helpers
