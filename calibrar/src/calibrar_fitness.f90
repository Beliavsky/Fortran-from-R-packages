! SPDX-License-Identifier: GPL-2.0-only
module calibrar_fitness
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use calibrar_kinds, only : dp
  implicit none
  private
  public :: fitness_norm2, fitness_normp, fitness_penalty, fitness_pois
  public :: fitness_lnorm2, fitness_lnorm3, fitness_lnorm4, fitness_lnorm4b
  public :: fitness_rangeq, fitness_multinom, weighted_sum_fitness

contains

  pure function weighted_sum_fitness(x, w) result(v)
    real(dp), intent(in) :: x(:), w(:)
    real(dp) :: v
    v = sum(x*w)
  end function weighted_sum_fitness

  function fitness_norm2(obs, sim) result(v)
    real(dp), intent(in) :: obs(:), sim(:)
    real(dp) :: v
    integer :: i
    v = 0.0_dp
    do i = 1, min(size(obs),size(sim))
      if (ieee_is_finite(obs(i)) .and. ieee_is_finite(sim(i))) v = v + (obs(i)-sim(i))**2
    end do
  end function fitness_norm2

  function fitness_normp(sim) result(v)
    real(dp), intent(in) :: sim(:)
    real(dp) :: v
    integer :: i
    v = 0.0_dp
    do i = 1, size(sim)
      if (ieee_is_finite(sim(i))) v = v + sim(i)*sim(i)
    end do
  end function fitness_normp

  function fitness_penalty(sim, n) result(v)
    real(dp), intent(in) :: sim(:)
    integer, intent(in), optional :: n
    real(dp) :: v
    integer :: i, nn, count
    nn = 100
    if (present(n)) nn = n
    v = 0.0_dp; count = 0
    do i = 1, size(sim)
      if (ieee_is_finite(sim(i))) then
        v = v + sim(i)*sim(i)
        count = count + 1
      end if
    end do
    if (count > 0) v = real(nn,dp)*v/real(count,dp)
  end function fitness_penalty

  function fitness_pois(obs, sim) result(v)
    real(dp), intent(in) :: obs(:), sim(:)
    real(dp) :: v
    integer :: i
    v = 0.0_dp
    do i = 1, min(size(obs),size(sim))
      if (ieee_is_finite(obs(i)) .and. ieee_is_finite(sim(i)) .and. sim(i) > 0.0_dp) then
        v = v - (obs(i)*log(sim(i)) - sim(i))
      end if
    end do
  end function fitness_pois

  function fitness_lnorm2(obs, sim, tiny_value) result(v)
    real(dp), intent(in) :: obs(:), sim(:)
    real(dp), intent(in), optional :: tiny_value
    real(dp) :: v, t
    integer :: i
    t = 1.0e-2_dp
    if (present(tiny_value)) t = tiny_value
    v = 0.0_dp
    do i = 1, min(size(obs),size(sim))
      if (ieee_is_finite(obs(i)) .and. ieee_is_finite(sim(i)) .and. obs(i)+t > 0.0_dp .and. sim(i)+t > 0.0_dp) then
        v = v + (log(obs(i)+t)-log(sim(i)+t))**2
      end if
    end do
  end function fitness_lnorm2

  function fitness_lnorm3(obs, sim, tiny_value) result(v)
    real(dp), intent(in) :: obs(:), sim(:)
    real(dp), intent(in), optional :: tiny_value
    real(dp) :: v, t, q, sr
    integer :: i, n
    t = 1.0e-2_dp
    if (present(tiny_value)) t = tiny_value
    sr = 0.0_dp; n = 0
    do i = 1, min(size(obs),size(sim))
      if (ieee_is_finite(obs(i)) .and. ieee_is_finite(sim(i)) .and. abs(sim(i)) > tiny(1.0_dp)) then
        sr = sr + obs(i)/sim(i); n = n+1
      end if
    end do
    if (n == 0) then
      v = huge(1.0_dp); return
    end if
    q = sr/real(n,dp)
    if (q <= 0.0_dp) then
      v = huge(1.0_dp); return
    end if
    v = 0.0_dp
    do i = 1, min(size(obs),size(sim))
      if (ieee_is_finite(obs(i)) .and. ieee_is_finite(sim(i)) .and. obs(i)+t > 0.0_dp .and. sim(i)+t > 0.0_dp) then
        v = v + (log(obs(i)+t)-log(sim(i)+t)-log(q))**2
      end if
    end do
  end function fitness_lnorm3

  function fitness_rangeq(obs, sim, b, c, dump, qout) result(pen)
    real(dp), intent(in) :: obs(:), sim(:)
    real(dp), intent(in), optional :: b, c
    logical, intent(in), optional :: dump
    real(dp), intent(out), optional :: qout
    real(dp) :: pen, bb, cc, q, lr, sr
    logical :: dd
    integer :: i, n
    bb=1.0_dp; cc=2.0_dp; dd=.true.
    if (present(b)) bb=b
    if (present(c)) cc=c
    if (present(dump)) dd=dump
    sr=0.0_dp; n=0
    do i=1,min(size(obs),size(sim))
      if (ieee_is_finite(obs(i)) .and. ieee_is_finite(sim(i)) .and. abs(sim(i)) > tiny(1.0_dp)) then
        sr=sr+obs(i)/sim(i); n=n+1
      end if
    end do
    if (n==0) then
      pen=huge(1.0_dp); q=1.0_dp
      if (present(qout)) qout=q
      return
    end if
    q=sr/real(n,dp)
    if (present(qout)) qout=q
    if (q <= 0.0_dp) then
      pen=huge(1.0_dp); return
    end if
    if (dd) then
      lr=abs(log(q)/log(2.0_dp))
      pen=real(n,dp)*(max(lr,bb)**cc-bb**cc)
    else
      pen=0.0_dp
      do i=1,min(size(obs),size(sim))
        if (ieee_is_finite(obs(i)) .and. ieee_is_finite(sim(i)) .and. abs(obs(i))>tiny(1.0_dp) .and. abs(sim(i))>tiny(1.0_dp)) then
          if (obs(i)/sim(i) > 0.0_dp) then
            lr=abs(log(obs(i)/sim(i))/log(2.0_dp))
            pen=pen+(max(lr,bb)**cc-bb**cc)
          end if
        end if
      end do
    end if
  end function fitness_rangeq

  function fitness_lnorm4(obs, sim, tiny_value, b, c) result(v)
    real(dp), intent(in) :: obs(:), sim(:)
    real(dp), intent(in), optional :: tiny_value, b, c
    real(dp) :: v, t, bb, cc, q, pen
    integer :: i
    t=1.0e-2_dp; bb=1.0_dp; cc=2.0_dp
    if(present(tiny_value)) t=tiny_value
    if(present(b)) bb=b
    if(present(c)) cc=c
    pen=fitness_rangeq(obs,sim,bb,cc,.true.,q)
    if(q<=0.0_dp) then; v=huge(1.0_dp); return; end if
    v=pen
    do i=1,min(size(obs),size(sim))
      if (ieee_is_finite(obs(i)) .and. ieee_is_finite(sim(i)) .and. obs(i)+t>0.0_dp .and. sim(i)+t>0.0_dp) &
        v=v+(log(obs(i)+t)-log(sim(i)+t)-log(q))**2
    end do
  end function fitness_lnorm4

  function fitness_lnorm4b(obs, sim, tiny_value, b, c) result(v)
    real(dp), intent(in) :: obs(:), sim(:)
    real(dp), intent(in), optional :: tiny_value, b, c
    real(dp) :: v, t, bb, cc, q, pen
    integer :: i
    t=1.0e-2_dp; bb=1.0_dp; cc=2.0_dp
    if(present(tiny_value)) t=tiny_value
    if(present(b)) bb=b
    if(present(c)) cc=c
    pen=fitness_rangeq(obs,sim,bb,cc,.false.,q)
    if(q<=0.0_dp) then; v=huge(1.0_dp); return; end if
    v=pen
    do i=1,min(size(obs),size(sim))
      if (ieee_is_finite(obs(i)) .and. ieee_is_finite(sim(i)) .and. obs(i)+t>0.0_dp .and. sim(i)+t>0.0_dp) &
        v=v+(log(obs(i)+t)-log(sim(i)+t)-log(q))**2
    end do
  end function fitness_lnorm4b

  function fitness_multinom(obs, sim, sample_size, tiny_value) result(v)
    real(dp), intent(in) :: obs(:,:), sim(:,:)
    integer, intent(in), optional :: sample_size
    real(dp), intent(in), optional :: tiny_value
    real(dp) :: v, tinyv, os, ss, po, ps, sig2, err
    integer :: i,j,nc,sz
    if(any(shape(obs)/=shape(sim))) error stop "fitness_multinom: shape mismatch"
    sz=20; tinyv=1.0e-3_dp
    if(present(sample_size)) sz=sample_size
    if(present(tiny_value)) tinyv=tiny_value
    nc=size(obs,2); v=0.0_dp
    do i=1,size(obs,1)
      os=sum(obs(i,:),mask=[(ieee_is_finite(obs(i,j)),j=1,nc)])
      ss=sum(sim(i,:),mask=[(ieee_is_finite(sim(i,j)),j=1,nc)])
      if(os<=0.0_dp) cycle
      if(ss<=0.0_dp) ss=real(nc,dp)
      do j=1,nc
        if(.not.ieee_is_finite(obs(i,j)) .or. .not.ieee_is_finite(sim(i,j))) cycle
        po=obs(i,j)/os
        if(sum(abs(sim(i,:)))<=tiny(1.0_dp)) then
          ps=1.0_dp/real(nc,dp)
        else
          ps=sim(i,j)/ss
        end if
        sig2=((1.0_dp-po)*po+1.0_dp/real(nc,dp))/real(sz,dp)
        err=log(exp(-(po-ps)**2/(2.0_dp*sig2))+tinyv)
        v=v-real(sz,dp)*err
      end do
    end do
  end function fitness_multinom
end module calibrar_fitness
