! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
module fportfolio_io
  use fportfolio_kinds, only: dp
  implicit none
  private
  public :: read_returns_csv, write_weights_csv
contains
  subroutine read_returns_csv(filename,data,names,has_header)
    character(len=*),intent(in)::filename
    real(dp),allocatable,intent(out)::data(:,:)
    character(len=64),allocatable,intent(out)::names(:)
    logical,intent(in),optional::has_header
    character(len=4096)::line
    character(len=256), allocatable :: tokens(:)
    integer::u,ios,nrow,ncol,nt,start,j,row
    logical::header,first_numeric
    real(dp)::dummy
    allocate(tokens(512))
    header=.true.;if(present(has_header))header=has_header
    open(newunit=u,file=filename,status='old',action='read',iostat=ios)
    if(ios/=0)error stop "read_returns_csv: cannot open file"
    read(u,'(a)',iostat=ios)line
    if(ios/=0)error stop "read_returns_csv: empty file"
    call split_csv(line,tokens,nt)
    if(header)then
      ncol=nt
      read(tokens(1),*,iostat=ios)dummy
      first_numeric=(ios==0)
      start=merge(1,2,first_numeric)
      ncol=nt-start+1
      allocate(names(ncol))
      do j=1,ncol
        names(j)=adjustl(tokens(start+j-1)(:64))
      end do
    else
      read(tokens(1),*,iostat=ios)dummy;first_numeric=(ios==0);start=merge(1,2,first_numeric)
      ncol=nt-start+1;allocate(names(ncol))
      do j=1,ncol;write(names(j),'("Asset",i0)')j;end do
      rewind(u)
    end if
    nrow=0
    do
      read(u,'(a)',iostat=ios)line
      if(ios/=0)exit
      if(len_trim(line)>0)nrow=nrow+1
    end do
    allocate(data(nrow,ncol));rewind(u)
    if(header)read(u,'(a)')line
    row=0
    do
      read(u,'(a)',iostat=ios)line
      if(ios/=0)exit
      if(len_trim(line)==0)cycle
      call split_csv(line,tokens,nt);row=row+1
      do j=1,ncol
        read(tokens(start+j-1),*,iostat=ios)data(row,j)
        if(ios/=0)error stop "read_returns_csv: nonnumeric field"
      end do
    end do
    close(u)
  end subroutine read_returns_csv

  subroutine write_weights_csv(filename,names,weights)
    character(len=*),intent(in)::filename
    character(len=*),intent(in)::names(:)
    real(dp),intent(in)::weights(:)
    integer::u,i
    open(newunit=u,file=filename,status='replace',action='write')
    write(u,'(a)')'Asset,Weight'
    do i=1,size(weights);write(u,'(a,",",es24.16)')trim(names(i)),weights(i);end do
    close(u)
  end subroutine write_weights_csv

  subroutine split_csv(line,tokens,n)
    character(len=*),intent(in)::line
    character(len=*),intent(out)::tokens(:)
    integer,intent(out)::n
    integer::i,start,l
    l=len_trim(line);n=0;start=1
    do i=1,l+1
      if(i==l+1 .or. line(i:i)==',')then
        n=n+1
        if(n>size(tokens))error stop "split_csv: too many columns"
        tokens(n)=adjustl(line(start:i-1));start=i+1
      end if
    end do
  end subroutine split_csv
end module fportfolio_io
