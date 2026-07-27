! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_io
   use robustbase_kinds, only: dp
   implicit none
   private
   public :: read_two_column_csv
contains
   subroutine read_two_column_csv(filename,x,y)
      character(len=*),intent(in)::filename
      real(dp),allocatable,intent(out)::x(:),y(:)
      character(len=4096)::line
      integer::unit,ios,n,i
      real(dp)::a,b
      n=0
      open(newunit=unit,file=filename,status='old',action='read',iostat=ios)
      if(ios/=0) error stop "read_two_column_csv: cannot open file"
      do
         read(unit,'(a)',iostat=ios)line
         if(ios/=0)exit
         call commas_to_spaces(line)
         read(line,*,iostat=ios)a,b
         if(ios==0)n=n+1
      end do
      rewind(unit)
      allocate(x(n),y(n));i=0
      do
         read(unit,'(a)',iostat=ios)line
         if(ios/=0)exit
         call commas_to_spaces(line)
         read(line,*,iostat=ios)a,b
         if(ios==0)then;i=i+1;x(i)=a;y(i)=b;end if
      end do
      close(unit)
      if(n==0) error stop "read_two_column_csv: no numeric rows"
   contains
      subroutine commas_to_spaces(s)
         character(len=*),intent(inout)::s
         integer::j
         do j=1,len_trim(s)
            if(s(j:j)==',')s(j:j)=' '
         end do
      end subroutine
   end subroutine read_two_column_csv
end module robustbase_io
