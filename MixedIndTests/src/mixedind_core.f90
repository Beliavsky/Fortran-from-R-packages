! SPDX-License-Identifier: GPL-3.0-only
module mixedind_core
  use mixedind_kinds, only : dp
  use mixedind_types, only : prepared_data_result, pair_dependence_result, sn_result, &
    bootstrap_result, mixedind_success, mixedind_invalid_argument, mixedind_numerical_error
  use mixedind_special, only : normal_quantile
  implicit none
  private

  public :: choose_int, count_nonserial_subsets, count_serial_subsets
  public :: preparedata_core, stat_dep_core, stat_dep_serial_core, ddn_scalar
  public :: sn_nonserial_core, sn_serial_core, sn_serial_vector_core
  public :: bootstrap_core, moebius_nonserial_core, moebius_serial_core

contains

  pure integer function choose_int(n, k) result(v)
    integer, intent(in) :: n, k
    integer :: i, kk
    if (k < 0 .or. k > n) then
      v = 0
      return
    end if
    kk = min(k, n-k)
    v = 1
    do i = 1, kk
      v = v * (n-kk+i) / i
    end do
  end function choose_int

  pure integer function count_nonserial_subsets(d, trunc_level) result(v)
    integer, intent(in) :: d, trunc_level
    integer :: k
    v = 0
    do k = 2, min(d,trunc_level)
      v = v + choose_int(d,k)
    end do
  end function count_nonserial_subsets

  pure integer function count_serial_subsets(p, trunc_level) result(v)
    integer, intent(in) :: p, trunc_level
    integer :: k
    v = 0
    do k = 2, min(p,trunc_level)
      v = v + choose_int(p-1,k-1)
    end do
  end function count_serial_subsets

  subroutine make_subsets(d, trunc_level, serial, subsets, cardinality)
    integer, intent(in) :: d, trunc_level
    logical, intent(in) :: serial
    integer, allocatable, intent(out) :: subsets(:,:), cardinality(:)
    integer :: mask, nmask, j, card, count, nsets

    if (serial) then
      nsets = count_serial_subsets(d,trunc_level)
      nmask = 2**(d-1)-1
    else
      nsets = count_nonserial_subsets(d,trunc_level)
      nmask = 2**d-1
    end if
    allocate(subsets(nsets,d), cardinality(nsets))
    subsets = 0
    cardinality = 0
    count = 0
    do mask = 1, nmask
      if (serial) then
        card = 1 + popcnt(mask)
        if (card > trunc_level) cycle
        count = count + 1
        subsets(count,1) = 1
        do j = 2, d
          subsets(count,j) = merge(1,0,btest(mask,j-2))
        end do
      else
        card = popcnt(mask)
        if (card < 2 .or. card > trunc_level) cycle
        count = count + 1
        do j = 1, d
          subsets(count,j) = merge(1,0,btest(mask,j-1))
        end do
      end if
      cardinality(count) = card
    end do
  end subroutine make_subsets

  function preparedata_core(x) result(out)
    real(dp), intent(in) :: x(:)
    type(prepared_data_result) :: out
    real(dp), allocatable :: tmp(:)
    integer :: i, j, m, n, ncount

    n = size(x)
    if (n < 1) then
      out%status = mixedind_invalid_argument
      return
    end if
    allocate(tmp(n))
    tmp = x
    call insertion_sort(tmp)
    m = 1
    do i = 2, n
      if (tmp(i) > tmp(i-1)) m = m + 1
    end do
    allocate(out%values(m), out%cdf(m), out%pdf(m))
    out%values(1) = tmp(1)
    j = 1
    do i = 2, n
      if (tmp(i) > tmp(i-1)) then
        j = j + 1
        out%values(j) = tmp(i)
      end if
    end do
    do j = 1, m
      ncount = count(x <= out%values(j))
      out%cdf(j) = real(ncount,dp)/real(n,dp)
    end do
    out%pdf(1) = out%cdf(1)
    if (m > 1) out%pdf(2:m) = out%cdf(2:m)-out%cdf(1:m-1)
  end function preparedata_core

  subroutine insertion_sort(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key
    do i = 2, size(x)
      key = x(i)
      j = i-1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j+1) = x(j)
        j = j-1
      end do
      x(j+1) = key
    end do
  end subroutine insertion_sort

  function stat_dep_core(x, y) result(out)
    real(dp), intent(in) :: x(:), y(:)
    type(pair_dependence_result) :: out
    integer :: n, i, j
    real(dp) :: sum0, sum1, sum2, sum3, a1, a2, c1, c2, s1, s2

    n = size(x)
    if (n < 2 .or. size(y) /= n) then
      out%status = mixedind_invalid_argument
      return
    end if
    sum0 = 0.0_dp
    sum3 = 0.0_dp
    s1 = 0.0_dp
    s2 = 0.0_dp
    do i = 1, n
      sum1 = 0.0_dp
      sum2 = 0.0_dp
      do j = 1, n
        a1 = real(merge(1,0,x(j) <= x(i)) + merge(1,0,x(j) < x(i)),dp)
        a2 = real(merge(1,0,y(j) <= y(i)) + merge(1,0,y(j) < y(i)),dp)
        sum0 = sum0 + a1*a2
        sum1 = sum1 + a1
        sum2 = sum2 + a2
      end do
      c1 = sum1/real(n,dp)-1.0_dp
      c2 = sum2/real(n,dp)-1.0_dp
      sum3 = sum3+c1*c2
      s1 = s1+c1*c1
      s2 = s2+c2*c2
    end do
    s1 = s1/real(n,dp)
    s2 = s2/real(n,dp)
    out%scale = sqrt(max(0.0_dp,s1*s2))
    out%tau = -1.0_dp+sum0/real(n*n,dp)
    if (out%scale > tiny(1.0_dp)) then
      out%rho = sum3/real(n,dp)/out%scale
    else
      out%rho = 0.0_dp
      out%status = mixedind_numerical_error
    end if
  end function stat_dep_core

  function stat_dep_serial_core(x, lag) result(out)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lag
    type(pair_dependence_result) :: out
    real(dp), allocatable :: y(:)
    integer :: n, i

    n = size(x)
    if (lag < 0 .or. lag >= n) then
      out%status = mixedind_invalid_argument
      return
    end if
    allocate(y(n))
    do i = 1, n
      y(i) = x(mod(i-1+lag,n)+1)
    end do
    out = stat_dep_core(x,y)
    if (out%status == mixedind_success .or. out%status == mixedind_numerical_error) then
      out%scale = population_midrank_variance(y)
      if (out%scale > tiny(1.0_dp)) then
        out%rho = raw_midrank_covariance(x,y)/out%scale
      else
        out%rho = 0.0_dp
      end if
    end if
  end function stat_dep_serial_core

  function population_midrank_variance(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v
    integer :: i, j, n
    real(dp) :: a, c
    n = size(x)
    v = 0.0_dp
    do i = 1, n
      a = 0.0_dp
      do j = 1, n
        a = a + real(merge(1,0,x(j)<=x(i))+merge(1,0,x(j)<x(i)),dp)
      end do
      c = a/real(n,dp)-1.0_dp
      v = v+c*c
    end do
    v = v/real(n,dp)
  end function population_midrank_variance

  function raw_midrank_covariance(x,y) result(v)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: v
    integer :: i, j, n
    real(dp) :: ax, ay, cx, cy
    n = size(x)
    v = 0.0_dp
    do i = 1, n
      ax = 0.0_dp
      ay = 0.0_dp
      do j = 1, n
        ax = ax+real(merge(1,0,x(j)<=x(i))+merge(1,0,x(j)<x(i)),dp)
        ay = ay+real(merge(1,0,y(j)<=y(i))+merge(1,0,y(j)<y(i)),dp)
      end do
      cx = ax/real(n,dp)-1.0_dp
      cy = ay/real(n,dp)-1.0_dp
      v = v+cx*cy
    end do
    v = v/real(n,dp)
  end function raw_midrank_covariance

  subroutine ifun_scalar(x, i1, i1point, i4)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: i1(:,:), i1point(:), i4(:,:)
    type(prepared_data_result) :: pd
    integer :: n, m, i, j, k
    real(dp) :: a, am, b, bm, s

    n = size(x)
    pd = preparedata_core(x)
    m = size(pd%values)
    do j = 1, n
      do i = 1, n
        s = 0.0_dp
        do k = 1, m
          a = real(merge(1,0,x(i)<=pd%values(k)),dp)
          am = real(merge(1,0,x(i)<pd%values(k)),dp)
          b = real(merge(1,0,x(j)<=pd%values(k)),dp)
          bm = real(merge(1,0,x(j)<pd%values(k)),dp)
          s = s+pd%pdf(k)*((a+am)*(b+bm)+am*bm+a*b)/6.0_dp
        end do
        i1(i,j) = s
      end do
      i1point(j) = sum(i1(:,j))/real(n,dp)
    end do
    do j = 1, n
      do i = 1, n
        i4(i,j) = i1(i,j)-i1point(j)-i1point(i)+1.0_dp/3.0_dp
      end do
    end do
  end subroutine ifun_scalar

  subroutine ddn_scalar(x, values, pdf, cdf, iv)
    real(dp), intent(in) :: x(:), values(:), pdf(:), cdf(:)
    real(dp), intent(out) :: iv(:,:)
    integer :: n, m, i, j, k
    real(dp) :: a, am, b, bm, aa, bb, s

    n = size(x)
    m = size(values)
    do j = 1, n
      do i = 1, n
        s = 1.0_dp/3.0_dp
        do k = 1, m
          a = real(merge(1,0,x(i)<=values(k)),dp)
          am = real(merge(1,0,x(i)<values(k)),dp)
          aa = a+am
          b = real(merge(1,0,x(j)<=values(k)),dp)
          bm = real(merge(1,0,x(j)<values(k)),dp)
          bb = b+bm
          s = s+pdf(k)*(-0.5_dp*cdf(k)*(aa+bb) + &
            (aa*bb+am*bm+a*b+pdf(k)*(aa+bb+am+bm))/6.0_dp)
        end do
        iv(i,j) = s
      end do
    end do
  end subroutine ddn_scalar

  subroutine sn_base(i1, i1point, cte, sn, jmat, i4)
    real(dp), intent(in) :: i1(:,:,:), i1point(:,:), cte
    real(dp), intent(out) :: sn
    real(dp), intent(out), optional :: jmat(:,:)
    real(dp), intent(in), optional :: i4(:,:,:)
    integer :: n, d, i, j, k, k1, k2, r, i1idx, j1idx, first1, first2
    integer, allocatable :: allsets(:,:), cards(:), first(:)
    real(dp), allocatable :: prodpoint(:), sumpoint(:)
    real(dp) :: ss, c0, c1, c2, s1, s2, s12, value

    n = size(i1,1)
    d = size(i1,3)
    c0 = cte**d
    allocate(prodpoint(n),sumpoint(n))
    do j = 1, n
      prodpoint(j) = product(i1point(j,:))
      sumpoint(j) = sum(i1point(j,:))
    end do
    sn = 0.0_dp
    do j = 1, n
      do i = 1, n
        ss = product(i1(i,j,:))
        sn = sn+ss-prodpoint(i)-prodpoint(j)+c0
      end do
    end do
    sn = sn/real(n,dp)
    if (.not. present(jmat)) return

    if (.not. present(i4)) then
      c1 = c0/cte
      c2 = c1/cte
      do j = 1, n
        do i = 1, n
          ss = product(i1(i,j,:))
          s1 = 0.0_dp
          s2 = 0.0_dp
          s12 = 0.0_dp
          do k = 1, d
            s12 = s12+i1point(i,k)*i1point(j,k)
            if (abs(i1point(j,k)) > tiny(1.0_dp)) &
              s1 = s1+prodpoint(j)*i1(i,j,k)/i1point(j,k)
            if (abs(i1point(i,k)) > tiny(1.0_dp)) &
              s1 = s1+i1(i,j,k)*prodpoint(i)/i1point(i,k)
            s2 = s2+i1(i,j,k)
          end do
          jmat(i,j) = ss-s1+s2*c1+(sumpoint(i)*sumpoint(j)-s12)*c2
        end do
      end do
      return
    end if

    call make_subsets(d,d,.false.,allsets,cards)
    allocate(first(size(cards)))
    do k = 1, size(cards)
      first(k) = 1
      do while (first(k) <= d .and. allsets(k,first(k)) == 0)
        first(k) = first(k)+1
      end do
    end do
    do j = 1, n
      do i = 1, n
        value = 0.0_dp
        do k1 = 1, size(cards)
          first1 = first(k1)-1
          do k2 = 1, size(cards)
            first2 = first(k2)-1
            ss = 1.0_dp
            do r = 0, d-1
              i1idx = modulo(i-1+first1-r,n)+1
              j1idx = modulo(j-1+first2-r,n)+1
              if (allsets(k1,r+1)==1 .and. allsets(k2,r+1)==1) then
                ss = ss*i4(i1idx,j1idx,1)
              else if (allsets(k1,r+1)==1) then
                ss = ss*(i1point(i1idx,1)-cte)
              else if (allsets(k2,r+1)==1) then
                ss = ss*(i1point(j1idx,1)-cte)
              else
                ss = ss*cte
              end if
            end do
            value = value+ss
          end do
        end do
        jmat(i,j) = value
      end do
    end do
  end subroutine sn_base

  subroutine subset_statistics(i4, subsets, stats, multiplier)
    real(dp), intent(in) :: i4(:,:,:)
    integer, intent(in) :: subsets(:,:)
    real(dp), intent(out) :: stats(:)
    real(dp), intent(out) :: multiplier(:,:,:)
    integer :: n, ns, d, i, j, k, l
    real(dp) :: ss

    n = size(i4,1)
    ns = size(subsets,1)
    d = size(subsets,2)
    do k = 1, ns
      stats(k) = 0.0_dp
      do j = 1, n
        do i = 1, n
          ss = 1.0_dp
          do l = 1, d
            if (subsets(k,l)==1) ss = ss*i4(i,j,l)
          end do
          multiplier(k,i,j) = ss
          stats(k) = stats(k)+ss
        end do
      end do
      stats(k) = stats(k)/real(n,dp)
    end do
  end subroutine subset_statistics

  function sn_nonserial_core(x, trunc_level) result(out)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: trunc_level
    type(sn_result) :: out
    real(dp), allocatable :: i1(:,:,:), i1point(:,:), i4(:,:,:)
    integer :: n, d, k

    n = size(x,1)
    d = size(x,2)
    if (n < 2 .or. d < 2 .or. trunc_level < 2) then
      out%status = mixedind_invalid_argument
      return
    end if
    call make_subsets(d,min(d,trunc_level),.false.,out%subsets,out%cardinality)
    allocate(out%stats(size(out%cardinality)),out%multiplier(size(out%cardinality),n,n), &
             out%sn_multiplier(n,n),i1(n,n,d),i1point(n,d),i4(n,n,d))
    do k = 1, d
      call ifun_scalar(x(:,k),i1(:,:,k),i1point(:,k),i4(:,:,k))
    end do
    call sn_base(i1,i1point,1.0_dp/3.0_dp,out%sn,out%sn_multiplier)
    call subset_statistics(i4,out%subsets,out%stats,out%multiplier)
  end function sn_nonserial_core

  function sn_serial_core(x, p, trunc_level, need_multipliers) result(out)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: p, trunc_level
    logical, intent(in), optional :: need_multipliers
    type(sn_result) :: out
    real(dp), allocatable :: lagx(:), i1(:,:,:), i1point(:,:), i4(:,:,:)
    integer :: n, k, i
    logical :: need

    n = size(x)
    need = .true.
    if (present(need_multipliers)) need = need_multipliers
    if (n < 2 .or. p < 2 .or. p > n .or. trunc_level < 2) then
      out%status = mixedind_invalid_argument
      return
    end if
    call make_subsets(p,min(p,trunc_level),.true.,out%subsets,out%cardinality)
    allocate(out%stats(size(out%cardinality)),out%multiplier(size(out%cardinality),n,n), &
             i1(n,n,p),i1point(n,p),i4(n,n,p),lagx(n))
    if (need) allocate(out%sn_multiplier(n,n))
    do k = 0, p-1
      do i = 1, n
        lagx(i) = x(modulo(i-1-k,n)+1)
      end do
      call ifun_scalar(lagx,i1(:,:,k+1),i1point(:,k+1),i4(:,:,k+1))
    end do
    if (need) then
      call sn_base(i1,i1point,1.0_dp/3.0_dp,out%sn,out%sn_multiplier,i4)
    else
      call sn_base(i1,i1point,1.0_dp/3.0_dp,out%sn)
    end if
    call subset_statistics(i4,out%subsets,out%stats,out%multiplier)
  end function sn_serial_core

  subroutine ifun_vector(x, i1, i1point, i4, d00)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: i1(:,:), i1point(:), i4(:,:), d00
    integer :: n, d, i, j, k
    real(dp), allocatable :: t1(:,:), tp(:), t4(:,:)

    n = size(x,1)
    d = size(x,2)
    i1 = 1.0_dp
    allocate(t1(n,n),tp(n),t4(n,n))
    do k = 1, d
      call ifun_scalar(x(:,k),t1,tp,t4)
      i1 = i1*t1
    end do
    do j = 1, n
      i1point(j) = sum(i1(:,j))/real(n,dp)
    end do
    d00 = sum(i1point)/real(n,dp)
    do j = 1, n
      do i = 1, n
        i4(i,j) = i1(i,j)-i1point(j)-i1point(i)+d00
      end do
    end do
  end subroutine ifun_vector

  function sn_serial_vector_core(x, p, trunc_level) result(out)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: p, trunc_level
    type(sn_result) :: out
    real(dp), allocatable :: lagx(:,:), i1(:,:,:), i1point(:,:), i4(:,:,:)
    real(dp) :: d00, d00_first
    integer :: n, d, k, i

    n = size(x,1)
    d = size(x,2)
    if (n < 2 .or. d < 1 .or. p < 2 .or. p > n .or. trunc_level < 2) then
      out%status = mixedind_invalid_argument
      return
    end if
    call make_subsets(p,min(p,trunc_level),.true.,out%subsets,out%cardinality)
    allocate(out%stats(size(out%cardinality)),out%multiplier(size(out%cardinality),n,n), &
             out%sn_multiplier(n,n),i1(n,n,p),i1point(n,p),i4(n,n,p),lagx(n,d))
    d00_first = 0.0_dp
    do k = 0, p-1
      do i = 1, n
        lagx(i,:) = x(modulo(i-1-k,n)+1,:)
      end do
      call ifun_vector(lagx,i1(:,:,k+1),i1point(:,k+1),i4(:,:,k+1),d00)
      if (k == 0) d00_first = d00
    end do
    call sn_base(i1,i1point,d00_first,out%sn,out%sn_multiplier,i4)
    call subset_statistics(i4,out%subsets,out%stats,out%multiplier)
  end function sn_serial_vector_core

  function bootstrap_core(multiplier, sn_multiplier, xi) result(out)
    real(dp), intent(in) :: multiplier(:,:,:), sn_multiplier(:,:), xi(:)
    type(bootstrap_result) :: out
    real(dp), allocatable :: xc(:)
    integer :: n, k

    n = size(xi)
    if (size(multiplier,2)/=n .or. size(multiplier,3)/=n .or. &
        size(sn_multiplier,1)/=n .or. size(sn_multiplier,2)/=n) then
      out%status = mixedind_invalid_argument
      return
    end if
    allocate(out%cvm(size(multiplier,1)),xc(n))
    xc = xi-sum(xi)/real(n,dp)
    do k = 1, size(out%cvm)
      out%cvm(k) = dot_product(xc,matmul(multiplier(k,:,:),xc))/real(n,dp)
    end do
    out%sn = dot_product(xc,matmul(sn_multiplier,xc))/real(n,dp)
  end function bootstrap_core

  subroutine cdf_at_observations(x, fn_le, fn_lt, mass)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: fn_le(:), fn_lt(:), mass(:)
    integer :: n, i
    n = size(x)
    do i = 1, n
      fn_le(i) = real(count(x<=x(i)),dp)/real(n,dp)
      fn_lt(i) = real(count(x<x(i)),dp)/real(n,dp)
      mass(i) = fn_le(i)-fn_lt(i)
    end do
  end subroutine cdf_at_observations

  elemental function lg_integral(u) result(v)
    real(dp), intent(in) :: u
    real(dp) :: v, z
    if (u > 0.0_dp .and. u < 1.0_dp) then
      z = normal_quantile(u)
      v = exp(-0.5_dp*z*z)/sqrt(2.0_dp*acos(-1.0_dp))
    else
      v = 0.0_dp
    end if
  end function lg_integral

  elemental function le_integral(u) result(v)
    real(dp), intent(in) :: u
    real(dp) :: v
    if (u > 0.0_dp) then
      v = u-u*log(u)
    else
      v = 0.0_dp
    end if
  end function le_integral

  subroutine moebius_scores(x, score_s, score_g, score_e, sd_s, sd_g, sd_e)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: score_s(:),score_g(:),score_e(:)
    real(dp), intent(out) :: sd_s,sd_g,sd_e
    real(dp), allocatable :: f1(:),f0(:),mass(:)
    integer :: n
    n = size(x)
    allocate(f1(n),f0(n),mass(n))
    call cdf_at_observations(x,f1,f0,mass)
    score_s = 0.5_dp*(f1+f0-1.0_dp)
    score_g = (lg_integral(f1)-lg_integral(f0))/mass
    score_e = (le_integral(f1)-le_integral(f0))/mass-1.0_dp
    sd_s = sqrt(sum((score_s-sum(score_s)/real(n,dp))**2)/real(n,dp))
    sd_g = sqrt(sum((score_g-sum(score_g)/real(n,dp))**2)/real(n,dp))
    sd_e = sqrt(sum((score_e-sum(score_e)/real(n,dp))**2)/real(n,dp))
  end subroutine moebius_scores

  subroutine moebius_nonserial_core(x,trunc_level,stat_s,stat_g,stat_e,cards,subsets,status)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: trunc_level
    real(dp), allocatable, intent(out) :: stat_s(:),stat_g(:),stat_e(:)
    integer, allocatable, intent(out) :: cards(:),subsets(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: ss(:,:),sg(:,:),se(:,:),sds(:),sdg(:),sde(:)
    real(dp) :: ps,pg,pe,ns,ng,ne
    integer :: n,d,i,j,k

    n=size(x,1); d=size(x,2)
    if(n<2 .or. d<2 .or. trunc_level<2) then
      status=mixedind_invalid_argument
      return
    end if
    call make_subsets(d,min(d,trunc_level),.false.,subsets,cards)
    allocate(stat_s(size(cards)),stat_g(size(cards)),stat_e(size(cards)), &
      ss(n,d),sg(n,d),se(n,d),sds(d),sdg(d),sde(d))
    do j=1,d
      call moebius_scores(x(:,j),ss(:,j),sg(:,j),se(:,j),sds(j),sdg(j),sde(j))
    end do
    do k=1,size(cards)
      ns=1.0_dp; ng=1.0_dp; ne=1.0_dp
      do j=1,d
        if(subsets(k,j)==1) then
          ns=ns*sds(j); ng=ng*sdg(j); ne=ne*sde(j)
        end if
      end do
      stat_s(k)=0.0_dp; stat_g(k)=0.0_dp; stat_e(k)=0.0_dp
      do i=1,n
        ps=1.0_dp; pg=1.0_dp; pe=1.0_dp
        do j=1,d
          if(subsets(k,j)==1) then
            ps=ps*ss(i,j); pg=pg*sg(i,j); pe=pe*se(i,j)
          end if
        end do
        stat_s(k)=stat_s(k)+ps
        stat_g(k)=stat_g(k)+pg
        stat_e(k)=stat_e(k)+pe
      end do
      stat_s(k)=stat_s(k)/real(n,dp)/ns
      stat_g(k)=stat_g(k)/real(n,dp)/ng
      stat_e(k)=stat_e(k)/real(n,dp)/ne
    end do
    status=mixedind_success
  end subroutine moebius_nonserial_core

  subroutine moebius_serial_core(y,p,trunc_level,stat_s,stat_g,stat_e,cards,subsets,status)
    real(dp), intent(in) :: y(:)
    integer, intent(in) :: p,trunc_level
    real(dp), allocatable, intent(out) :: stat_s(:),stat_g(:),stat_e(:)
    integer, allocatable, intent(out) :: cards(:),subsets(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: lagy(:),ss(:,:),sg(:,:),se(:,:)
    real(dp), allocatable :: sds(:),sdg(:),sde(:)
    real(dp) :: ps,pg,pe,ns,ng,ne
    integer :: n,i,j,k

    n=size(y)
    if(n<2 .or. p<2 .or. p>n .or. trunc_level<2) then
      status=mixedind_invalid_argument
      return
    end if
    call make_subsets(p,min(p,trunc_level),.true.,subsets,cards)
    allocate(stat_s(size(cards)),stat_g(size(cards)),stat_e(size(cards)), &
      lagy(n),ss(n,p),sg(n,p),se(n,p),sds(p),sdg(p),sde(p))
    do j=0,p-1
      do i=1,n
        lagy(i)=y(modulo(i-1-j,n)+1)
      end do
      call moebius_scores(lagy,ss(:,j+1),sg(:,j+1),se(:,j+1),sds(j+1),sdg(j+1),sde(j+1))
    end do
    do k=1,size(cards)
      ns=1.0_dp; ng=1.0_dp; ne=1.0_dp
      do j=1,p
        if(subsets(k,j)==1) then
          ns=ns*sds(j); ng=ng*sdg(j); ne=ne*sde(j)
        end if
      end do
      stat_s(k)=0.0_dp; stat_g(k)=0.0_dp; stat_e(k)=0.0_dp
      do i=1,n
        ps=1.0_dp; pg=1.0_dp; pe=1.0_dp
        do j=1,p
          if(subsets(k,j)==1) then
            ps=ps*ss(i,j); pg=pg*sg(i,j); pe=pe*se(i,j)
          end if
        end do
        stat_s(k)=stat_s(k)+ps
        stat_g(k)=stat_g(k)+pg
        stat_e(k)=stat_e(k)+pe
      end do
      stat_s(k)=stat_s(k)/real(n,dp)/ns
      stat_g(k)=stat_g(k)/real(n,dp)/ng
      stat_e(k)=stat_e(k)/real(n,dp)/ne
    end do
    status=mixedind_success
  end subroutine moebius_serial_core

end module mixedind_core
