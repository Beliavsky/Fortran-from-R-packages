! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
module rq_dates
  use rq_kinds, only: dp
  implicit none
  private
  public :: date_t, calendar_t, schedule_t
  public :: make_date, date_to_serial, serial_to_date, add_days, add_months
  public :: day_count, year_fraction, weekday, is_weekend, is_business_day
  public :: adjust_date, advance_date, business_days_between, make_schedule
  public :: following, modified_following, preceding, modified_preceding, unadjusted

  integer, parameter :: following=1, modified_following=2, preceding=3
  integer, parameter :: modified_preceding=4, unadjusted=5

  type :: date_t
    integer :: year=1970, month=1, day=1
  contains
    procedure :: serial => date_serial_method
  end type date_t

  type :: calendar_t
    type(date_t), allocatable :: holidays(:)
  contains
    procedure :: is_business_day => calendar_is_business
  end type calendar_t

  type :: schedule_t
    type(date_t), allocatable :: dates(:)
  end type schedule_t
contains
  pure function make_date(year,month,day) result(d)
    integer,intent(in)::year,month,day
    type(date_t)::d
    d%year=year; d%month=month; d%day=day
  end function make_date

  pure integer function date_serial_method(self) result(z)
    class(date_t),intent(in)::self
    z=date_to_serial(self)
  end function date_serial_method

  pure integer function date_to_serial(d) result(z)
    type(date_t),intent(in)::d
    integer :: y,m,era,yoe,doy,doe
    y=d%year; m=d%month
    if(m<=2) y=y-1
    era=floor_div(y,400)
    yoe=y-era*400
    if(m>2) then
      doy=(153*(m-3)+2)/5+d%day-1
    else
      doy=(153*(m+9)+2)/5+d%day-1
    end if
    doe=yoe*365+yoe/4-yoe/100+doy
    z=era*146097+doe-719468
  end function date_to_serial

  pure function serial_to_date(z) result(d)
    integer,intent(in)::z
    type(date_t)::d
    integer :: zz,era,doe,yoe,y,doy,mp
    zz=z+719468
    era=floor_div(zz,146097)
    doe=zz-era*146097
    yoe=(doe-doe/1460+doe/36524-doe/146096)/365
    y=yoe+era*400
    doy=doe-(365*yoe+yoe/4-yoe/100)
    mp=(5*doy+2)/153
    d%day=doy-(153*mp+2)/5+1
    if(mp<10) then
      d%month=mp+3
    else
      d%month=mp-9
    end if
    if(d%month<=2) y=y+1
    d%year=y
  end function serial_to_date

  pure integer function floor_div(a,b) result(q)
    integer,intent(in)::a,b
    q=a/b
    if(mod(a,b)/=0 .and. ((a<0).neqv.(b<0))) q=q-1
  end function floor_div

  pure function add_days(d,n) result(out)
    type(date_t),intent(in)::d
    integer,intent(in)::n
    type(date_t)::out
    out=serial_to_date(date_to_serial(d)+n)
  end function add_days

  pure logical function is_leap(year) result(ok)
    integer,intent(in)::year
    ok=mod(year,4)==0 .and. (mod(year,100)/=0 .or. mod(year,400)==0)
  end function is_leap

  pure integer function days_in_month(year,month) result(n)
    integer,intent(in)::year,month
    integer,parameter::dm(12)=[31,28,31,30,31,30,31,31,30,31,30,31]
    n=dm(month)
    if(month==2 .and. is_leap(year)) n=29
  end function days_in_month

  pure function add_months(d,n,end_of_month) result(out)
    type(date_t),intent(in)::d
    integer,intent(in)::n
    logical,intent(in),optional::end_of_month
    type(date_t)::out
    integer :: total,y,m,dd
    logical :: eom
    eom=.false.; if(present(end_of_month)) eom=end_of_month
    total=d%year*12+(d%month-1)+n
    y=floor_div(total,12); m=total-y*12+1
    dd=min(d%day,days_in_month(y,m))
    if(eom .and. d%day==days_in_month(d%year,d%month)) dd=days_in_month(y,m)
    out=make_date(y,m,dd)
  end function add_months

  pure integer function weekday(d) result(w)
    type(date_t),intent(in)::d
    ! Monday=1, ..., Sunday=7. 1970-01-01 was Thursday.
    w=modulo(date_to_serial(d)+3,7)+1
  end function weekday

  pure logical function is_weekend(d) result(ok)
    type(date_t),intent(in)::d
    ok=weekday(d)>=6
  end function is_weekend

  pure logical function same_date(a,b) result(ok)
    type(date_t),intent(in)::a,b
    ok=a%year==b%year .and. a%month==b%month .and. a%day==b%day
  end function same_date

  logical function is_business_day(d,calendar) result(ok)
    type(date_t),intent(in)::d
    type(calendar_t),intent(in),optional::calendar
    integer::i
    ok=.not.is_weekend(d)
    if(.not.ok .or. .not.present(calendar)) return
    if(allocated(calendar%holidays)) then
      do i=1,size(calendar%holidays)
        if(same_date(d,calendar%holidays(i))) then; ok=.false.; return; end if
      end do
    end if
  end function is_business_day

  logical function calendar_is_business(self,d) result(ok)
    class(calendar_t),intent(in)::self
    type(date_t),intent(in)::d
    ok=is_business_day(d,self)
  end function calendar_is_business

  function adjust_date(d,convention,calendar) result(out)
    type(date_t),intent(in)::d
    integer,intent(in)::convention
    type(calendar_t),intent(in),optional::calendar
    type(date_t)::out,tmp
    integer::m0
    out=d
    if(convention==unadjusted .or. is_business_day(out,calendar)) return
    select case(convention)
    case(following,modified_following)
      m0=d%month
      do while(.not.is_business_day(out,calendar)); out=add_days(out,1); end do
      if(convention==modified_following .and. out%month/=m0) then
        out=d
        do while(.not.is_business_day(out,calendar)); out=add_days(out,-1); end do
      end if
    case(preceding,modified_preceding)
      m0=d%month
      do while(.not.is_business_day(out,calendar)); out=add_days(out,-1); end do
      if(convention==modified_preceding .and. out%month/=m0) then
        out=d
        do while(.not.is_business_day(out,calendar)); out=add_days(out,1); end do
      end if
    case default
      tmp=d; out=tmp
    end select
  end function adjust_date

  function advance_date(d,amount,unit,convention,calendar,end_of_month) result(out)
    type(date_t),intent(in)::d
    integer,intent(in)::amount
    character(len=*),intent(in)::unit
    integer,intent(in),optional::convention
    type(calendar_t),intent(in),optional::calendar
    logical,intent(in),optional::end_of_month
    type(date_t)::out
    integer::c,step,count
    c=following; if(present(convention)) c=convention
    select case(trim(adjustl(unit)))
    case('day','days')
      out=add_days(d,amount)
    case('businessday','businessdays')
      out=d; step=merge(1,-1,amount>=0); count=0
      do while(count<abs(amount))
        out=add_days(out,step)
        if(is_business_day(out,calendar)) count=count+1
      end do
    case('week','weeks')
      out=add_days(d,7*amount)
    case('month','months')
      out=add_months(d,amount,end_of_month)
    case('year','years')
      out=add_months(d,12*amount,end_of_month)
    case default
      out=d
    end select
    out=adjust_date(out,c,calendar)
  end function advance_date

  integer function business_days_between(from,to,calendar,include_first,include_last) result(n)
    type(date_t),intent(in)::from,to
    type(calendar_t),intent(in),optional::calendar
    logical,intent(in),optional::include_first,include_last
    integer::s1,s2,k,step
    logical::ifirst,ilast
    ifirst=.true.; ilast=.false.
    if(present(include_first)) ifirst=include_first
    if(present(include_last)) ilast=include_last
    s1=date_to_serial(from); s2=date_to_serial(to); step=merge(1,-1,s2>=s1); n=0
    do k=s1,s2,step
      if(k==s1 .and. .not.ifirst) cycle
      if(k==s2 .and. .not.ilast) cycle
      if(is_business_day(serial_to_date(k),calendar)) n=n+step
    end do
  end function business_days_between

  pure integer function day_count(start_date,end_date,basis) result(n)
    type(date_t),intent(in)::start_date,end_date
    character(len=*),intent(in)::basis
    integer::d1,d2,m1,m2,y1,y2
    select case(trim(adjustl(basis)))
    case('30/360','Thirty360','30E/360')
      y1=start_date%year; m1=start_date%month; d1=min(start_date%day,30)
      y2=end_date%year; m2=end_date%month; d2=min(end_date%day,30)
      n=360*(y2-y1)+30*(m2-m1)+(d2-d1)
    case default
      n=date_to_serial(end_date)-date_to_serial(start_date)
    end select
  end function day_count

  pure real(dp) function year_fraction(start_date,end_date,basis) result(yf)
    type(date_t),intent(in)::start_date,end_date
    character(len=*),intent(in)::basis
    type(date_t)::jan1,nextjan
    integer::y,days
    yf=0.0_dp
    select case(trim(adjustl(basis)))
    case('Actual360')
      yf=real(day_count(start_date,end_date,'Actual'),dp)/360.0_dp
    case('Actual365','Actual365Fixed')
      yf=real(day_count(start_date,end_date,'Actual'),dp)/365.0_dp
    case('30/360','Thirty360','30E/360')
      yf=real(day_count(start_date,end_date,'30/360'),dp)/360.0_dp
    case('ActualActual','ActualActualISDA')
      if(start_date%year==end_date%year) then
        yf=real(day_count(start_date,end_date,'Actual'),dp)/merge(366.0_dp,365.0_dp,is_leap(start_date%year))
      else
        jan1=make_date(start_date%year+1,1,1)
        yf=real(day_count(start_date,jan1,'Actual'),dp)/merge(366.0_dp,365.0_dp,is_leap(start_date%year))
        do y=start_date%year+1,end_date%year-1; yf=yf+1.0_dp; end do
        nextjan=make_date(end_date%year,1,1)
        yf=yf+real(day_count(nextjan,end_date,'Actual'),dp)/merge(366.0_dp,365.0_dp,is_leap(end_date%year))
      end if
    case default
      days=day_count(start_date,end_date,'Actual'); yf=real(days,dp)/365.0_dp
    end select
  end function year_fraction

  subroutine make_schedule(effective,maturity,months,calendar,convention,termination_convention, &
                           forward_generation,end_of_month,schedule)
    type(date_t),intent(in)::effective,maturity
    integer,intent(in)::months
    type(calendar_t),intent(in),optional::calendar
    integer,intent(in),optional::convention,termination_convention
    logical,intent(in),optional::forward_generation,end_of_month
    type(schedule_t),intent(out)::schedule
    type(date_t),allocatable::tmp(:)
    type(date_t)::d
    integer::n,i,c,tc
    logical::fwd,eom
    c=following; tc=following; fwd=.true.; eom=.false.
    if(present(convention)) c=convention
    if(present(termination_convention)) tc=termination_convention
    if(present(forward_generation)) fwd=forward_generation
    if(present(end_of_month)) eom=end_of_month
    if(months<=0) error stop 'make_schedule: months must be positive'
    n=1
    if(fwd) then
      d=effective
      do while(date_to_serial(add_months(d,months,eom))<date_to_serial(maturity)); n=n+1; d=add_months(d,months,eom); end do
      n=n+1; allocate(tmp(n)); tmp(1)=adjust_date(effective,c,calendar); d=effective
      do i=2,n-1; d=add_months(d,months,eom); tmp(i)=adjust_date(d,c,calendar); end do
      tmp(n)=adjust_date(maturity,tc,calendar)
    else
      d=maturity
      do while(date_to_serial(add_months(d,-months,eom))>date_to_serial(effective)); n=n+1; d=add_months(d,-months,eom); end do
      n=n+1; allocate(tmp(n)); tmp(n)=adjust_date(maturity,tc,calendar); d=maturity
      do i=n-1,2,-1; d=add_months(d,-months,eom); tmp(i)=adjust_date(d,c,calendar); end do
      tmp(1)=adjust_date(effective,c,calendar)
    end if
    call move_alloc(tmp,schedule%dates)
  end subroutine make_schedule
end module rq_dates
