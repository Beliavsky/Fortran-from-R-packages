! SPDX-License-Identifier: GPL-2.0-or-later
module lifeinsurer_profit
  use lifeinsurer_kinds, only : dp, lir_success, lir_dimension_error
  use lifeinsurer_types, only : profit_rate_table, profit_values, reserve_values, premium_values
  implicit none
  private
  public :: shift_by, pp_base_null, pp_base_previous_zillmer_reserve
  public :: pp_base_zillmer_reserve_t2, pp_base_contractual_reserve
  public :: pp_base_previous_contractual_reserve, pp_base_mean_contractual_reserve
  public :: pp_base_sum_insured, pp_calculate_rate_on_base
  public :: pp_calculate_rate_on_base_min0, pp_calculate_rate_plus_guarantee_on_base
  public :: pp_calculate_rate_on_base_sgf_factor, calculate_profit_participation
  public :: sum_profits
contains
  function shift_by(x,n) result(y)
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: n
    real(dp), allocatable :: y(:)
    integer :: k
    k=1; if(present(n)) k=max(0,n)
    allocate(y(size(x))); y=0.0_dp
    if(k<size(x)) y(k+1:)=x(1:size(x)-k)
  end function

  function pp_base_null(n) result(x)
    integer,intent(in)::n; real(dp),allocatable::x(:); allocate(x(max(0,n))); x=0.0_dp
  end function
  function pp_base_previous_zillmer_reserve(r) result(x)
    type(reserve_values),intent(in)::r; real(dp),allocatable::x(:); x=shift_by(r%zillmer,1)
  end function
  function pp_base_zillmer_reserve_t2(r) result(x)
    type(reserve_values),intent(in)::r; real(dp),allocatable::x(:); x=shift_by(r%zillmer,2)
  end function
  function pp_base_contractual_reserve(r) result(x)
    type(reserve_values),intent(in)::r; real(dp),allocatable::x(:); allocate(x(size(r%contractual))); x=max(0.0_dp,r%contractual)
  end function
  function pp_base_previous_contractual_reserve(r) result(x)
    type(reserve_values),intent(in)::r; real(dp),allocatable::x(:); x=shift_by(r%contractual,1)
  end function
  function pp_base_mean_contractual_reserve(r) result(x)
    type(reserve_values),intent(in)::r; real(dp),allocatable::x(:); real(dp),allocatable::z(:)
    integer::n
    n=size(r%contractual); allocate(z(n+1)); z(1)=0.0_dp; z(2:)=r%contractual
    allocate(x(n)); x=max(0.0_dp,0.5_dp*(z(1:n)+z(2:n+1)))
  end function
  function pp_base_sum_insured(sum_insured,n) result(x)
    real(dp),intent(in)::sum_insured; integer,intent(in)::n; real(dp),allocatable::x(:)
    allocate(x(max(0,n))); x=sum_insured
  end function

  pure function pp_calculate_rate_on_base(base,rate,waiting) result(x)
    real(dp),intent(in)::base(:),rate(:),waiting(:); real(dp)::x(size(base)); x=base*rate*waiting
  end function
  pure function pp_calculate_rate_on_base_min0(base,rate,waiting) result(x)
    real(dp),intent(in)::base(:),rate(:),waiting(:); real(dp)::x(size(base)); x=max(0.0_dp,base*rate*waiting)
  end function
  pure function pp_calculate_rate_plus_guarantee_on_base(base,rate,guaranteed,waiting) result(x)
    real(dp),intent(in)::base(:),rate(:),guaranteed(:),waiting(:); real(dp)::x(size(base)); x=base*(rate+guaranteed)*waiting
  end function
  pure function pp_calculate_rate_on_base_sgf_factor(base,rate,waiting,sgf) result(x)
    real(dp),intent(in)::base(:),rate(:),waiting(:),sgf(:); real(dp)::x(size(base)); x=base*rate*waiting*(1.0_dp-sgf)
  end function

  subroutine calculate_profit_participation(rates,res,prem,sum_insured,waiting,out,status)
    type(profit_rate_table),intent(in)::rates
    type(reserve_values),intent(in)::res
    type(premium_values),intent(in)::prem
    real(dp),intent(in)::sum_insured
    real(dp),intent(in),optional::waiting(:)
    type(profit_values),intent(out)::out
    integer,intent(out)::status
    real(dp),allocatable::w(:),bint(:),brisk(:),bexp(:),bsum(:)
    integer::n,i
    if(.not.allocated(res%contractual) .or. .not.allocated(rates%interest_profit)) then
      status=lir_dimension_error; return
    end if
    n=size(res%contractual)
    if(size(rates%interest_profit)/=n) then; status=lir_dimension_error; return; end if
    allocate(w(n)); w=1.0_dp; if(present(waiting)) then
      if(size(waiting)/=n) then; status=lir_dimension_error; return; end if
      w=waiting
    end if
    bint=pp_base_previous_contractual_reserve(res)
    brisk=shift_by(max(0.0_dp,prem%zillmer-res%zillmer),1)
    bexp=pp_base_previous_contractual_reserve(res)
    bsum=pp_base_sum_insured(sum_insured,n)
    allocate(out%interest(n),out%risk(n),out%expense(n),out%sum_component(n))
    allocate(out%interest_on_profit(n),out%terminal_bonus(n),out%terminal_bonus_fund(n))
    allocate(out%total_assignment(n),out%accumulated(n),out%benefit(n))
    out%interest=max(0.0_dp,bint*rates%interest_profit*w)
    out%risk=max(0.0_dp,brisk*rates%mortality_profit*w)
    out%expense=max(0.0_dp,bexp*rates%expense_profit*w)
    out%sum_component=max(0.0_dp,bsum*rates%sum_profit*w)
    out%interest_on_profit=0.0_dp; out%terminal_bonus=0.0_dp; out%terminal_bonus_fund=0.0_dp
    out%accumulated=0.0_dp
    do i=1,n
      if(i>1) out%interest_on_profit(i)=out%accumulated(i-1)*(rates%interest_profit(i)+rates%guaranteed_interest(i))*w(i)
      out%total_assignment(i)=out%interest(i)+out%risk(i)+out%expense(i)+out%sum_component(i)+out%interest_on_profit(i)
      if(allocated(rates%terminal_bonus_fund)) then
        out%terminal_bonus_fund(i)=out%total_assignment(i)*rates%terminal_bonus_fund(i)
      end if
      out%accumulated(i)=out%total_assignment(i)-out%terminal_bonus_fund(i)
      if(i>1) out%accumulated(i)=out%accumulated(i)+out%accumulated(i-1)
      if(allocated(rates%terminal_bonus)) out%terminal_bonus(i)=out%accumulated(i)*rates%terminal_bonus(i)
      out%benefit(i)=out%accumulated(i)+out%terminal_bonus(i)
    end do
    status=lir_success
  end subroutine

  pure real(dp) function sum_profits(values)
    real(dp), intent(in) :: values(:)
    sum_profits=sum(values)
  end function
end module lifeinsurer_profit
