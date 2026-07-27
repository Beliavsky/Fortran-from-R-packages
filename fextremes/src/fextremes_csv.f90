! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software: you can redistribute it
! and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 2 of the License,
! or (at your option) any later version.
module fextremes_csv
  use fextremes_kinds, only: dp
  implicit none
  private
  public :: read_numeric_series
contains
  subroutine read_numeric_series(filename,x,ok)
    character(len=*),intent(in)::filename
    real(dp),allocatable,intent(out)::x(:)
    logical,intent(out)::ok
    character(len=2048)::line,clean,token
    real(dp),allocatable::tmp(:)
    real(dp)::v
    integer::unit,ios,n,i,p
    ok=.false.; n=0
    open(newunit=unit,file=filename,status='old',action='read',iostat=ios)
    if(ios/=0) then; allocate(x(0)); return; end if
    do
      read(unit,'(A)',iostat=ios) line
      if(ios/=0) exit
      clean=line
      do i=1,len_trim(clean); if(clean(i:i)==',' .or. clean(i:i)==';') clean(i:i)=' '; end do
      read(clean,*,iostat=ios) v
      if(ios==0) then; n=n+1; cycle; end if
      read(clean,*,iostat=ios) token,v
      if(ios==0) n=n+1
    end do
    rewind(unit); allocate(tmp(n)); p=0
    do
      read(unit,'(A)',iostat=ios) line
      if(ios/=0) exit
      clean=line
      do i=1,len_trim(clean); if(clean(i:i)==',' .or. clean(i:i)==';') clean(i:i)=' '; end do
      read(clean,*,iostat=ios) v
      if(ios/=0) read(clean,*,iostat=ios) token,v
      if(ios==0) then; p=p+1; tmp(p)=v; end if
    end do
    close(unit); allocate(x(p)); x=tmp(:p); ok=p>0
  end subroutine read_numeric_series
end module fextremes_csv
