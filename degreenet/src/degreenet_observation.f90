! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet_observation
  use degreenet_kinds, only : dp, huge_neg
  use degreenet_models, only : model_pmf, grouped_probability
  implicit none
  private
  public :: grouped_loglik, rounded_loglik, rounded_probability, rounded_bin

contains
  real(dp) function grouped_loglik(model,par,x,cutoff) result(ll)
    integer,intent(in)::model,x(:),cutoff
    real(dp),intent(in)::par(:)
    integer::i
    real(dp)::p
    ll=0.0_dp
    do i=1,size(x)
      p=grouped_probability(model,par,x(i),cutoff)
      if(p<=0.0_dp)then;ll=huge_neg;return;end if
      ll=ll+log(p)
    end do
  end function grouped_loglik

  subroutine rounded_bin(obs,lo,hi)
    integer,intent(in)::obs
    integer,intent(out)::lo,hi
    lo=obs;hi=obs
    select case(obs)
    case(5);lo=4;hi=6
    case(10);lo=7;hi=13
    case(12);lo=8;hi=16
    case(15);lo=12;hi=18
    case(20);lo=16;hi=24
    case(25);lo=20;hi=29
    case(30);lo=24;hi=36
    case(35);lo=25;hi=44
    case(40);lo=30;hi=49
    case(45);lo=35;hi=54
    case(50);lo=35;hi=64
    case(55);lo=40;hi=69
    case(60);lo=45;hi=74
    case(65);lo=50;hi=79
    case(70);lo=55;hi=84
    case(75);lo=60;hi=89
    case(80);lo=65;hi=94
    case(85);lo=70;hi=99
    case(90);lo=75;hi=104
    case(95);lo=80;hi=109
    case(100);lo=75;hi=124
    case(120);lo=95;hi=144
    case(130);lo=105;hi=154
    case(150);lo=115;hi=184
    case(200);lo=150;hi=249
    case(300);lo=250;hi=349
    case(560);lo=500;hi=639
    case(800);lo=500;hi=1100
    end select
  end subroutine rounded_bin

  real(dp) function rounded_probability(model,par,obs,cutoff) result(p)
    integer,intent(in)::model,obs,cutoff
    real(dp),intent(in)::par(:)
    integer::lo,hi,k
    call rounded_bin(obs,lo,hi)
    lo=max(lo,cutoff);p=0.0_dp
    do k=lo,hi;p=p+model_pmf(model,par,k,cutoff);end do
  end function rounded_probability

  real(dp) function rounded_loglik(model,par,x,cutoff,cutabove) result(ll)
    integer,intent(in)::model,x(:),cutoff,cutabove
    real(dp),intent(in)::par(:)
    integer::i
    real(dp)::p
    ll=0.0_dp
    do i=1,size(x)
      if(x(i)<cutoff.or.x(i)>cutabove)cycle
      p=rounded_probability(model,par,x(i),cutoff)
      if(p<=0.0_dp)then;ll=huge_neg;return;end if
      ll=ll+log(p)
    end do
  end function rounded_loglik
end module degreenet_observation
