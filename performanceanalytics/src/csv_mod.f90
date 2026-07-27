! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module csv_mod
  use kinds_mod, only: dp
  implicit none
  private
  public :: read_dated_csv
contains
  subroutine read_dated_csv(filename,data,status)
    character(len=*),intent(in)::filename
    real(dp),allocatable,intent(out)::data(:,:)
    integer,intent(out)::status
    character(len=8192)::line
    integer::unit,ios,nrow,nfield,ncol,i
    status=0;nrow=0;nfield=0
    open(newunit=unit,file=filename,status='old',action='read',iostat=ios)
    if(ios/=0)then;status=1;allocate(data(0,0));return;end if
    read(unit,'(a)',iostat=ios)line
    if(ios/=0)then;status=2;close(unit);allocate(data(0,0));return;end if
    nfield=1+count_commas(trim(line));ncol=max(0,nfield-1)
    do
      read(unit,'(a)',iostat=ios)line
      if(ios/=0)exit
      if(len_trim(line)>0)nrow=nrow+1
    end do
    rewind(unit);read(unit,'(a)')line
    allocate(data(nrow,ncol));data=0.0_dp
    do i=1,nrow
      read(unit,'(a)',iostat=ios)line
      if(ios/=0)then;status=3;exit;end if
      call parse_numeric_fields(line,data(i,:),ios)
      if(ios/=0)then;status=4;exit;end if
    end do
    close(unit)
  end subroutine read_dated_csv

  pure integer function count_commas(line) result(n)
    character(len=*),intent(in)::line
    integer::i
    n=0
    do i=1,len_trim(line)
      if(line(i:i)==',')n=n+1
    end do
  end function count_commas

  subroutine parse_numeric_fields(line,row,status)
    character(len=*),intent(in)::line
    real(dp),intent(out)::row(:)
    integer,intent(out)::status
    character(len=256)::token
    integer::start,stop,field,ios,n
    status=0;n=len_trim(line);start=1;field=0
    do
      stop=index(line(start:n),',')
      if(stop==0)then
        token=adjustl(line(start:n));field=field+1
        if(field>1 .and. field-1<=size(row))then
          read(token,*,iostat=ios)row(field-1)
          if(ios/=0)status=1
        end if
        exit
      else
        stop=start+stop-2
        token=adjustl(line(start:stop));field=field+1
        if(field>1 .and. field-1<=size(row))then
          read(token,*,iostat=ios)row(field-1)
          if(ios/=0)then;status=1;return;end if
        end if
        start=stop+2
        if(start>n+1)exit
      end if
    end do
  end subroutine parse_numeric_fields
end module csv_mod
