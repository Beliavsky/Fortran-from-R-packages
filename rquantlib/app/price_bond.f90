! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
program price_bond
  use rq_kinds, only: dp
  use rq_bonds
  implicit none
  character(len=64) :: arg
  real(dp) :: face,coupon,maturity,yield
  integer :: frequency
  type(bond_result) :: result
  if(command_argument_count()<5) then
    write(*,'(a)') 'usage: price_bond face coupon_rate maturity frequency yield'
    error stop 1
  end if
  call get_command_argument(1,arg); read(arg,*) face
  call get_command_argument(2,arg); read(arg,*) coupon
  call get_command_argument(3,arg); read(arg,*) maturity
  call get_command_argument(4,arg); read(arg,*) frequency
  call get_command_argument(5,arg); read(arg,*) yield
  call fixed_rate_bond_from_yield(face,coupon,maturity,frequency,yield,result)
  write(*,'(a,f16.8)') 'npv:               ',result%npv
  write(*,'(a,f16.8)') 'clean price:       ',result%clean_price
  write(*,'(a,f16.8)') 'duration:          ',result%duration
  write(*,'(a,f16.8)') 'modified duration: ',result%modified_duration
  write(*,'(a,f16.8)') 'convexity:         ',result%convexity
end program price_bond
