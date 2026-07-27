! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
program price_option
  use rq_kinds, only: dp
  use rq_options
  implicit none
  character(len=32) :: style,option_type,arg
  real(dp) :: spot,strike,q,r,t,sigma
  type(option_result) :: result
  integer :: status
  if(command_argument_count()<8) then
    write(*,'(a)') 'usage: price_option european|american call|put spot strike q r maturity volatility'
    error stop 1
  end if
  call get_command_argument(1,style)
  call get_command_argument(2,option_type)
  call get_command_argument(3,arg); read(arg,*) spot
  call get_command_argument(4,arg); read(arg,*) strike
  call get_command_argument(5,arg); read(arg,*) q
  call get_command_argument(6,arg); read(arg,*) r
  call get_command_argument(7,arg); read(arg,*) t
  call get_command_argument(8,arg); read(arg,*) sigma
  select case(trim(style))
  case('european')
    result=european_option(option_type,spot,strike,q,r,t,sigma)
  case('american')
    result=american_option(option_type,spot,strike,q,r,t,sigma,600)
  case default
    write(*,'(a)') 'unknown style: '//trim(style)
    error stop 2
  end select
  status=0
  write(*,'(a,f16.8)') 'value: ',result%value
  write(*,'(a,f16.8)') 'delta: ',result%delta
  write(*,'(a,f16.8)') 'gamma: ',result%gamma
  write(*,'(a,f16.8)') 'vega:  ',result%vega
  if(status/=0) error stop 3
end program price_option
