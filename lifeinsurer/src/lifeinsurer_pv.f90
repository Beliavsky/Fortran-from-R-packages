! SPDX-License-Identifier: GPL-2.0-or-later
module lifeinsurer_pv
  use lifeinsurer_kinds, only : dp, lir_success, lir_dimension_error, lir_invalid_argument
  use lifeinsurer_types, only : transition_probabilities, frequency_correction, cash_flow_set, present_value_set
  implicit none
  private
  public :: pv_guaranteed, pv_survival, pv_death, pv_disease, pv_after_death
  public :: calculate_present_values
contains
  subroutine check_two(a,b,status)
    real(dp), intent(in) :: a(:),b(:)
    integer, intent(out) :: status
    status=lir_success
    if(size(a)/=size(b)) status=lir_dimension_error
  end subroutine

  function pv_guaranteed(advance,arrears,v,m,corr,status) result(res)
    real(dp), intent(in) :: advance(:),arrears(:),v
    integer, intent(in) :: m
    type(frequency_correction), intent(in) :: corr
    integer, intent(out), optional :: status
    real(dp), allocatable :: res(:)
    real(dp) :: nxt,advc,arrc
    integer :: i,st
    call check_two(advance,arrears,st); allocate(res(size(advance))); res=0.0_dp
    if(st/=lir_success .or. m<=0) then
      if(present(status)) status=merge(st,lir_invalid_argument,st/=lir_success)
      return
    end if
    advc=corr%alpha-corr%beta*(1.0_dp-v)
    arrc=corr%alpha-(corr%beta+1.0_dp/real(m,dp))*(1.0_dp-v)
    nxt=0.0_dp
    do i=size(res),1,-1
      res(i)=advance(i)*advc+arrears(i)*arrc+v*nxt; nxt=res(i)
    end do
    if(present(status)) status=lir_success
  end function

  function pv_survival(advance,arrears,px,v,m,corr,status) result(res)
    real(dp), intent(in) :: advance(:),arrears(:),px(:),v
    integer, intent(in) :: m
    type(frequency_correction), intent(in) :: corr
    integer, intent(out), optional :: status
    real(dp), allocatable :: res(:)
    real(dp) :: nxt,advc,arrc,p
    integer :: i,n
    n=size(advance); allocate(res(n)); res=0.0_dp
    if(size(arrears)/=n .or. size(px)<max(0,n-1) .or. m<=0) then
      if(present(status)) status=lir_dimension_error; return
    end if
    nxt=0.0_dp
    do i=n,1,-1
      p=0.0_dp; if(i<=size(px)) p=px(i)
      advc=corr%alpha-corr%beta*(1.0_dp-p*v)
      arrc=corr%alpha-(corr%beta+1.0_dp/real(m,dp))*(1.0_dp-p*v)
      res(i)=advance(i)*advc+arrears(i)*arrc+v*p*nxt; nxt=res(i)
    end do
    if(present(status)) status=lir_success
  end function

  function pv_death(benefits,qx,px,v,status) result(res)
    real(dp), intent(in) :: benefits(:),qx(:),px(:),v
    integer, intent(out), optional :: status
    real(dp), allocatable :: res(:)
    real(dp) :: nxt,q,p
    integer :: i,n
    n=size(benefits); allocate(res(n)); res=0.0_dp
    if(size(qx)<max(0,n-1) .or. size(px)<max(0,n-1)) then
      if(present(status)) status=lir_dimension_error; return
    end if
    nxt=0.0_dp
    do i=n,1,-1
      q=0.0_dp; p=0.0_dp
      if(i<=size(qx)) q=qx(i)
      if(i<=size(px)) p=px(i)
      res(i)=v*q*benefits(i)+v*p*nxt; nxt=res(i)
    end do
    if(present(status)) status=lir_success
  end function

  function pv_disease(benefits,ix,px,v,status) result(res)
    real(dp), intent(in) :: benefits(:),ix(:),px(:),v
    integer, intent(out), optional :: status
    real(dp), allocatable :: res(:)
    real(dp) :: nxt,q,p
    integer :: i,n
    n=size(benefits); allocate(res(n)); res=0.0_dp
    if(size(ix)<max(0,n-1) .or. size(px)<max(0,n-1)) then
      if(present(status)) status=lir_dimension_error; return
    end if
    nxt=0.0_dp
    do i=n,1,-1
      q=0.0_dp; p=0.0_dp
      if(i<=size(ix)) q=ix(i)
      if(i<=size(px)) p=px(i)
      res(i)=v*q*benefits(i)+v*p*nxt; nxt=res(i)
    end do
    if(present(status)) status=lir_success
  end function

  function pv_after_death(advance,arrears,qx,px,v,m,corr,status) result(res)
    real(dp), intent(in) :: advance(:),arrears(:),qx(:),px(:),v
    integer, intent(in) :: m
    type(frequency_correction), intent(in) :: corr
    integer, intent(out), optional :: status
    real(dp), allocatable :: res(:)
    real(dp) :: live_next,dead_next,q,p,advc,arrc,new_live,new_dead
    integer :: i,n
    n=size(advance); allocate(res(n)); res=0.0_dp
    if(size(arrears)/=n .or. size(qx)<max(0,n-1) .or. size(px)<max(0,n-1) .or. m<=0) then
      if(present(status)) status=lir_dimension_error; return
    end if
    advc=corr%alpha-corr%beta*(1.0_dp-v)
    arrc=corr%alpha-(corr%beta+1.0_dp/real(m,dp))*(1.0_dp-v)
    live_next=0.0_dp; dead_next=0.0_dp
    do i=n,1,-1
      q=0.0_dp; p=0.0_dp
      if(i<=size(qx)) q=qx(i)
      if(i<=size(px)) p=px(i)
      new_live=p*v*live_next+q*v*dead_next
      new_dead=advance(i)*advc+arrears(i)*arrc+v*dead_next
      res(i)=new_live; live_next=new_live; dead_next=new_dead
    end do
    if(present(status)) status=lir_success
  end function

  subroutine calculate_present_values(cf,tr,v,premium_m,premium_corr,benefit_m,benefit_corr,pv,status)
    type(cash_flow_set), intent(in) :: cf
    type(transition_probabilities), intent(in) :: tr
    real(dp), intent(in) :: v
    integer, intent(in) :: premium_m,benefit_m
    type(frequency_correction), intent(in) :: premium_corr,benefit_corr
    type(present_value_set), intent(out) :: pv
    integer, intent(out) :: status
    integer :: st
    pv%premiums=pv_survival(cf%premiums_advance,cf%premiums_arrears,tr%px,v,premium_m,premium_corr,st)
    if(st/=lir_success) then; status=st; return; end if
    pv%additional_capital=pv_survival(cf%additional_capital,0.0_dp*cf%additional_capital,tr%px,v,1,benefit_corr,st)
    pv%guaranteed=pv_guaranteed(cf%guaranteed_advance,cf%guaranteed_arrears,v,benefit_m,benefit_corr,st)
    pv%survival=pv_survival(cf%survival_advance,cf%survival_arrears,tr%px,v,benefit_m,benefit_corr,st)
    pv%death_sum_insured=pv_death(cf%death_sum_insured,tr%qx,tr%px,v,st)
    pv%disease_sum_insured=pv_disease(cf%disease_sum_insured,tr%ix,tr%px,v,st)
    pv%death_gross_premium=pv_death(cf%death_gross_premium,tr%qx,tr%px,v,st)
    pv%death_refund_past=pv_death(cf%death_refund_past,tr%qx,tr%px,v,st)
    pv%benefits=pv%guaranteed+pv%survival+pv%death_sum_insured+pv%disease_sum_insured
    status=lir_success
  end subroutine
end module lifeinsurer_pv
