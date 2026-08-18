! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
module rmkdiscrete_manaclash
  use rmkdiscrete_kinds, only : dp
  use rmkdiscrete_math, only : qnan, ninf, real_equal, log_choose, dbinom_prob, rbinom_simple
  implicit none
  private
  public :: dmanaclash_dmg, dmanaclash_xyn, dmanaclash_net, rmanaclash, rmanaclash_sample
contains

  pure subroutine normalize_probs(pa,pb,pc,pd,a,b,c,d,ok)
    real(dp), intent(in) :: pa,pb,pc,pd
    real(dp), intent(out) :: a,b,c,d
    logical, intent(out) :: ok
    real(dp) :: s
    ok = pa > 0.0_dp .and. pb > 0.0_dp .and. pc > 0.0_dp .and. pd > 0.0_dp
    if (.not. ok) then
      a=qnan()
      b=qnan()
      c=qnan()
      d=qnan()
      return
    end if
    s=pa+pb+pc+pd
    if (.not. (s > 0.0_dp)) then
      ok=.false.
      a=qnan()
      b=qnan()
      c=qnan()
      d=qnan()
      return
    end if
    a=pa/s
    b=pb/s
    c=pc/s
    d=pd/s
  end subroutine normalize_probs

  subroutine get_probs(pa,pb,pc,pd,a,b,c,d,ok)
    real(dp), intent(in), optional :: pa,pb,pc,pd
    real(dp), intent(out) :: a,b,c,d
    logical, intent(out) :: ok
    real(dp) :: aa,bb,cc,dd
    aa=0.25_dp
    bb=0.25_dp
    cc=0.25_dp
    dd=0.25_dp
    if(present(pa)) aa=pa
    if(present(pb)) bb=pb
    if(present(pc)) cc=pc
    if(present(pd)) dd=pd
    call normalize_probs(aa,bb,cc,dd,a,b,c,d,ok)
  end subroutine get_probs

  real(dp) function dmanaclash_xyn(x,y,n,pa,pb,pc,pd,give_log) result(v)
    integer, intent(in) :: x,y,n
    real(dp), intent(in), optional :: pa,pb,pc,pd
    logical, intent(in), optional :: give_log
    real(dp) :: a,b,c,d,lv
    logical :: ok,gl
    gl=.false.
    if(present(give_log))gl=give_log
    call get_probs(pa,pb,pc,pd,a,b,c,d,ok)
    if(.not.ok) then
    v=qnan()
    return
    end if
    if(x<0 .or. y<0 .or. n<0 .or. x>n .or. y>n .or. x+y-n<0 .or. x+y-n>x) then
      v=merge(ninf(),0.0_dp,gl)
      return
    end if
    lv=log(d)-real(x,dp)*log(b)+real(n,dp)*log(b)-real(n,dp)*log(a) &
       +real(x+y,dp)*log(a)-real(y,dp)*log(c)+real(n,dp)*log(c) &
       +log_choose(n,x)+log_choose(x,x+y-n)
    v=merge(lv,exp(lv),gl)
  end function dmanaclash_xyn

  real(dp) function dmanaclash_dmg(x,y,n,pa,pb,pc,pd,give_log) result(v)
    integer, intent(in) :: x,y
    integer, intent(in), optional :: n
    real(dp), intent(in), optional :: pa,pb,pc,pd
    logical, intent(in), optional :: give_log
    real(dp) :: a,b,c,d,r,s,lv
    logical :: ok,gl
    integer :: i
    gl=.false.
    if(present(give_log))gl=give_log
    call get_probs(pa,pb,pc,pd,a,b,c,d,ok)
    if(.not.ok) then
    v=qnan()
    return
    end if
    if(x<0 .or. y<0) then
    v=merge(ninf(),0.0_dp,gl)
    return
    end if
    r=a+b+c
    if(present(n)) then
      lv=dbinom_prob(x,n,(a+c)/r,.true.) + dbinom_prob(x+y-n,x,a/(a+c),.true.)
      if(gl) then
      v=lv
      else
      v=exp(lv)
      end if
      return
    end if
    s=0.0_dp
    do i=max(x,y),x+y
      s=s+dmanaclash_xyn(x,y,i,a,b,c,d)
    end do
    if(gl) then
      if(s>0.0_dp) then
      v=log(s)
      else
      v=ninf()
      end if
    else
      v=s
    end if
  end function dmanaclash_dmg

  real(dp) function dmanaclash_net(z,pa,pb,pc,pd,rel_eps,give_log) result(v)
    integer, intent(in) :: z
    real(dp), intent(in), optional :: pa,pb,pc,pd,rel_eps
    logical, intent(in), optional :: give_log
    real(dp) :: a,b,c,d,eps,s,newterm
    logical :: ok,gl
    integer :: x,y,iter
    gl=.false.
    if(present(give_log))gl=give_log
    eps=1.0e-8_dp
    if(present(rel_eps))eps=rel_eps
    call get_probs(pa,pb,pc,pd,a,b,c,d,ok)
    if(.not.ok .or. eps<=0.0_dp) then
    v=qnan()
    return
    end if
    if(z>0) then
    x=z
    y=0
    else if(z<0) then
    y=-z
    x=0
    else
    x=0
    y=0
    end if
    s=dmanaclash_dmg(x,y,pa=a,pb=b,pc=c,pd=d)
    do iter=1,1000000
      x=x+1
      y=y+1
      newterm=dmanaclash_dmg(x,y,pa=a,pb=b,pc=c,pd=d)
      if(s>0.0_dp) then
        if(newterm/s<=eps) then
          s=s+newterm
          exit
        end if
      else if(real_equal(newterm,0.0_dp)) then
        exit
      end if
      s=s+newterm
    end do
    if(gl) then
      if(s>0.0_dp) then
      v=log(s)
      else
      v=ninf()
      end if
    else
      v=s
    end if
  end function dmanaclash_net

  subroutine get_probs_rng(pa,pb,pc,pd,a,b,c,d,ok)
    real(dp), intent(in), optional :: pa,pb,pc,pd
    real(dp), intent(out) :: a,b,c,d
    logical, intent(out) :: ok
    real(dp) :: aa,bb,cc,dd,s
    aa=0.25_dp
    bb=0.25_dp
    cc=0.25_dp
    dd=0.25_dp
    if(present(pa)) aa=pa
    if(present(pb)) bb=pb
    if(present(pc)) cc=pc
    if(present(pd)) dd=pd
    ok=aa>=0.0_dp .and. bb>=0.0_dp .and. cc>=0.0_dp .and. dd>0.0_dp
    if(.not.ok) then
      a=qnan()
      b=qnan()
      c=qnan()
      d=qnan()
      return
    end if
    s=aa+bb+cc+dd
    a=aa/s
    b=bb/s
    c=cc/s
    d=dd/s
  end subroutine get_probs_rng

  subroutine rmanaclash(pa,pb,pc,pd,out,n)
    real(dp), intent(in), optional :: pa,pb,pc,pd
    integer, intent(out) :: out(3)
    integer, intent(in), optional :: n
    real(dp) :: a,b,c,d,r,u
    logical :: ok
    integer :: x,y,nn
    call get_probs_rng(pa,pb,pc,pd,a,b,c,d,ok)
    if(.not.ok) then
    out=-huge(1)
    return
    end if
    r=a+b+c
    if(present(n)) then
      if(n<0) then
      out=-huge(1)
      return
      end if
      x=rbinom_simple(n,(a+c)/r)
      y=rbinom_simple(x,a/(a+c))+n-x
      out=[x,y,n]
      return
    end if
    x=0
    y=0
    nn=0
    do
      call random_number(u)
      if(u<=a) then
        x=x+1
        y=y+1
        nn=nn+1
      else if(u<=a+b) then
        y=y+1
        nn=nn+1
      else if(u<=a+b+c) then
        x=x+1
        nn=nn+1
      else
        exit
      end if
    end do
    out=[x,y,nn]
  end subroutine rmanaclash

  function rmanaclash_sample(nsim,pa,pb,pc,pd,n) result(out)
    integer, intent(in) :: nsim
    real(dp), intent(in), optional :: pa,pb,pc,pd
    integer, intent(in), optional :: n
    integer, allocatable :: out(:,:)
    integer :: i
    allocate(out(max(0,nsim),3))
    do i=1,nsim
      if(present(n)) then
        call rmanaclash(pa,pb,pc,pd,out(i,:),n)
      else
        call rmanaclash(pa,pb,pc,pd,out(i,:))
      end if
    end do
  end function rmanaclash_sample
end module rmkdiscrete_manaclash
