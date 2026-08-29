module forecast_calendar
   use forecast_kinds, only : dp
   implicit none
   private
   public :: is_leap_year, month_days_sequence, easter_gregorian, easter_effect, business_days_sequence
contains
   pure logical function is_leap_year(year) result(leap)
      integer,intent(in)::year
      leap=mod(year,4)==0 .and. (mod(year,100)/=0 .or. mod(year,400)==0)
   end function is_leap_year

   pure integer function days_in_month(year,month) result(nd)
      integer,intent(in)::year,month
      integer,parameter::base(12)=[31,28,31,30,31,30,31,31,30,31,30,31]
      nd=base(month)
      if(month==2 .and. is_leap_year(year))nd=29
   end function days_in_month

   pure integer function date_ordinal(year,month,day) result(ord)
      integer,intent(in)::year,month,day
      integer::y,m
      ord=365*(year-1)+(year-1)/4-(year-1)/100+(year-1)/400+day
      do m=1,month-1
         ord=ord+days_in_month(year,m)
      end do
   end function date_ordinal

   subroutine period_bounds(year,season,freq,lo,hi)
      integer,intent(in)::year,season,freq
      integer,intent(out)::lo,hi
      integer::m1,m2
      if(freq==12)then
         lo=date_ordinal(year,season,1)
         hi=date_ordinal(year,season,days_in_month(year,season))
      else if(freq==4)then
         m1=3*(season-1)+1
         m2=m1+2
         lo=date_ordinal(year,m1,1)
         hi=date_ordinal(year,m2,days_in_month(year,m2))
      else
         error stop 'period_bounds: frequency must be 12 or 4'
      end if
   end subroutine period_bounds

   function month_days_sequence(start_year,start_season,frequency,n) result(days)
      integer,intent(in)::start_year,start_season,frequency,n
      integer,allocatable::days(:)
      integer::i,y,s,lo,hi
      if(frequency/=12 .and. frequency/=4)error stop 'month_days_sequence: frequency must be 12 or 4'
      allocate(days(n))
      y=start_year
      s=start_season
      do i=1,n
         call period_bounds(y,s,frequency,lo,hi)
         days(i)=hi-lo+1
         s=s+1
         if(s>frequency)then
         s=1
         y=y+1
         end if
      end do
   end function month_days_sequence

   pure subroutine easter_gregorian(year,month,day)
      integer,intent(in)::year
      integer,intent(out)::month,day
      integer::a,b,c,d,e,f,g,h,i,k,l,m
      a=mod(year,19)
      b=year/100
      c=mod(year,100)
      d=b/4
      e=mod(b,4)
      f=(b+8)/25
      g=(b-f+1)/3
      h=mod(19*a+b-d-g+15,30)
      i=c/4
      k=mod(c,4)
      l=mod(32+2*e+2*i-h-k,7)
      m=(a+11*h+22*l)/451
      month=(h+l-7*m+114)/31
      day=mod(h+l-7*m+114,31)+1
   end subroutine easter_gregorian

   function easter_effect(start_year,start_season,frequency,n,easter_monday) result(effect)
      integer,intent(in)::start_year,start_season,frequency,n
      logical,intent(in),optional::easter_monday
      real(dp),allocatable::effect(:)
      integer::i,y,s,lo,hi,ey,em,ed,esun,elo,ehi,overlap,den
      logical::mon
      mon=.false.
      if(present(easter_monday))mon=easter_monday
      den=merge(4,3,mon)
      allocate(effect(n))
      effect=0.0_dp
      y=start_year
      s=start_season
      do i=1,n
         call period_bounds(y,s,frequency,lo,hi)
         do ey=y-1,y+1
            call easter_gregorian(ey,em,ed)
            esun=date_ordinal(ey,em,ed)
            elo=esun-2
            ehi=esun+merge(1,0,mon)
            overlap=max(0,min(hi,ehi)-max(lo,elo)+1)
            if(overlap>0)effect(i)=effect(i)+real(overlap,dp)/real(den,dp)
         end do
         s=s+1
         if(s>frequency)then
         s=1
         y=y+1
         end if
      end do
   end function easter_effect

   function business_days_sequence(start_year,start_season,frequency,n,holiday_dates) result(days)
      ! Numerical core of forecast::bizdays. Named financial-centre holiday
      ! calendars come from R's timeDate dependency; callers pass those dates
      ! explicitly here as rows [year,month,day]. Weekends are always excluded.
      integer,intent(in)::start_year,start_season,frequency,n
      integer,intent(in),optional::holiday_dates(:,:)
      integer,allocatable::days(:)
      integer::i,y,s,lo,hi,ord,h,j
      logical::holiday
      if(present(holiday_dates))then
         if(size(holiday_dates,2)/=3)error stop 'business_days_sequence: holidays must have three columns'
      end if
      allocate(days(n))
      y=start_year
      s=start_season
      do i=1,n
         call period_bounds(y,s,frequency,lo,hi)
         days(i)=0
         do ord=lo,hi
            ! 0001-01-01 is Monday in the proleptic Gregorian calendar.
            if(mod(ord-1,7)>=5)cycle
            holiday=.false.
            if(present(holiday_dates))then
               do j=1,size(holiday_dates,1)
                  h=date_ordinal(holiday_dates(j,1),holiday_dates(j,2),holiday_dates(j,3))
                  if(h==ord)then
                  holiday=.true.
                  exit
                  end if
               end do
            end if
            if(.not.holiday)days(i)=days(i)+1
         end do
         s=s+1
         if(s>frequency)then
         s=1
         y=y+1
         end if
      end do
   end function business_days_sequence
end module forecast_calendar
