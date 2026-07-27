! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from PeerPerformance 2.4.0, copyright 2012-2023 David Ardia and Kris Boudt.
module peerperformance_bootstrap
  use peerperformance_kinds, only: dp
  use peerperformance_math, only: finite_value, missing_value, set_random_seed, &
                                  random_integer, random_geometric, clamp_probability
  use peerperformance_stats, only: sharpe, modified_sharpe, sharpe_difference, &
                                   modified_sharpe_difference, sharpe_standard_error, &
                                   modified_sharpe_standard_error
  use peerperformance_linalg, only: ols_fit
  use peerperformance_types, only: test_result
  implicit none
  private
  public :: bootstrap_indices, sharpe_testing_bootstrap
  public :: modified_sharpe_testing_bootstrap
  public :: sharpe_block_size, modified_sharpe_block_size

contains

  subroutine bootstrap_indices(n, n_boot, block_length, indices, seed)
    integer, intent(in) :: n, n_boot, block_length
    integer, allocatable, intent(out) :: indices(:,:)
    integer, intent(in), optional :: seed
    integer :: b, j, block, start, pos, length_use, sd
    if (n < 1 .or. n_boot < 1 .or. block_length < 1 .or. block_length > n) then
      allocate(indices(0,0))
      return
    end if
    sd=12345; if (present(seed)) sd=seed
    call set_random_seed(sd)
    allocate(indices(n,n_boot))
    if (block_length == 1) then
      do b=1,n_boot
        do j=1,n
          indices(j,b)=random_integer(n)
        end do
      end do
    else
      do b=1,n_boot
        pos=1
        do while (pos <= n)
          start=random_integer(n)
          length_use=min(block_length,n-pos+1)
          do block=0,length_use-1
            indices(pos+block,b)=1+modulo(start-1+block,n)
          end do
          pos=pos+length_use
        end do
      end do
    end if
  end subroutine bootstrap_indices

  subroutine sharpe_testing_bootstrap(x, y, result, n_boot, block_length, ttype, &
                                      p_boot, seed, min_obs, null_difference)
    real(dp), intent(in) :: x(:), y(:)
    type(test_result), intent(out) :: result
    integer, intent(in), optional :: n_boot, block_length, ttype, p_boot, seed, min_obs
    real(dp), intent(in), optional :: null_difference
    logical, allocatable :: mask(:)
    real(dp), allocatable :: xx(:), yy(:), bx(:), by(:), bstats(:), m(:,:)
    integer, allocatable :: indices(:,:)
    real(dp) :: d, d0, se, bd, bse, tvalue, vals(2)
    integer :: nb, bl, kind_test, ptype, sd, minimum, n, b, valid, upper, lower, counts(2)
    result%status=0; result%message=''
    nb=499; if (present(n_boot)) nb=n_boot
    bl=1; if (present(block_length)) bl=block_length
    kind_test=2; if (present(ttype)) kind_test=ttype
    ptype=1; if (present(p_boot)) ptype=p_boot
    sd=12345; if (present(seed)) sd=seed
    minimum=10; if (present(min_obs)) minimum=min_obs
    d0=0.0_dp; if (present(null_difference)) d0=null_difference
    if (size(x) /= size(y)) then
      result%status=1; result%message='x and y must have equal length'; return
    end if
    allocate(mask(size(x))); mask=finite_value(x) .and. finite_value(y)
    n=count(mask); result%n=n
    if (n < minimum .or. bl < 1 .or. bl > n) then
      result%status=2; result%message='invalid sample or block length'; return
    end if
    xx=pack(x,mask); yy=pack(y,mask)
    d=sharpe_difference(xx,yy,kind_test)-d0
    se=sharpe_standard_error(xx,yy,.false.,kind_test,bl)
    if (.not. finite_value(d) .or. .not. finite_value(se) .or. se <= tiny(1.0_dp)) then
      result%status=3; result%message='degenerate Sharpe comparison'; return
    end if
    call bootstrap_indices(n,nb,bl,indices,sd)
    allocate(bx(n),by(n),bstats(nb)); bstats=missing_value()
    valid=0
    do b=1,nb
      bx=xx(indices(:,b)); by=yy(indices(:,b))
      bd=sharpe_difference(bx,by,kind_test)
      bse=sharpe_standard_error(bx,by,.false.,kind_test,bl)
      if (.not. finite_value(bd) .or. .not. finite_value(bse) .or. bse <= tiny(1.0_dp)) cycle
      valid=valid+1
      if (ptype == 1) then
        bstats(valid)=abs(bd-(d+d0))/bse
      else
        bstats(valid)=(bd-(d+d0))/bse
      end if
    end do
    if (valid < 1) then
      result%status=4; result%message='all bootstrap replicates were degenerate'; return
    end if
    tvalue=d/se
    if (ptype == 1) then
      upper=count(bstats(1:valid) >= abs(tvalue))
      d0=real(upper+1,dp)/real(valid+1,dp)
    else
      upper=count(bstats(1:valid) > tvalue)
      lower=count(bstats(1:valid) < tvalue)
      d0=2.0_dp*real(min(upper+1,lower+1),dp)/real(valid+1,dp)
    end if
    allocate(m(n,2)); m(:,1)=xx; m(:,2)=yy
    call sharpe(m,vals,counts)
    allocate(result%estimate(1,2),result%difference(1),result%standard_error(1), &
             result%tstat(1),result%pvalue(1))
    result%estimate(1,:)=vals
    result%difference(1)=vals(1)-vals(2)
    result%standard_error(1)=se
    result%tstat(1)=tvalue
    result%pvalue(1)=clamp_probability(d0)
  end subroutine sharpe_testing_bootstrap

  subroutine modified_sharpe_testing_bootstrap(x, y, level, result, na_negative, &
                                               n_boot, block_length, ttype, p_boot, &
                                               seed, min_obs, null_difference)
    real(dp), intent(in) :: x(:), y(:), level
    type(test_result), intent(out) :: result
    logical, intent(in), optional :: na_negative
    integer, intent(in), optional :: n_boot, block_length, ttype, p_boot, seed, min_obs
    real(dp), intent(in), optional :: null_difference
    logical, allocatable :: mask(:)
    real(dp), allocatable :: xx(:), yy(:), bx(:), by(:), bstats(:), m(:,:)
    integer, allocatable :: indices(:,:)
    real(dp) :: d, d0, se, bd, bse, tvalue, vals(2), pval
    integer :: nb, bl, kind_test, ptype, sd, minimum, n, b, valid, upper, lower, counts(2)
    logical :: reject_negative
    result%status=0; result%message=''
    reject_negative=.true.; if (present(na_negative)) reject_negative=na_negative
    nb=499; if (present(n_boot)) nb=n_boot
    bl=1; if (present(block_length)) bl=block_length
    kind_test=2; if (present(ttype)) kind_test=ttype
    ptype=1; if (present(p_boot)) ptype=p_boot
    sd=12345; if (present(seed)) sd=seed
    minimum=10; if (present(min_obs)) minimum=min_obs
    d0=0.0_dp; if (present(null_difference)) d0=null_difference
    if (size(x) /= size(y)) then
      result%status=1; result%message='x and y must have equal length'; return
    end if
    allocate(mask(size(x))); mask=finite_value(x) .and. finite_value(y)
    n=count(mask); result%n=n
    if (n < minimum .or. bl < 1 .or. bl > n) then
      result%status=2; result%message='invalid sample or block length'; return
    end if
    xx=pack(x,mask); yy=pack(y,mask)
    d=modified_sharpe_difference(xx,yy,level,reject_negative,kind_test)-d0
    se=modified_sharpe_standard_error(xx,yy,level,.false.,kind_test,bl)
    if (.not. finite_value(d) .or. .not. finite_value(se) .or. se <= tiny(1.0_dp)) then
      result%status=3; result%message='degenerate modified-Sharpe comparison'; return
    end if
    call bootstrap_indices(n,nb,bl,indices,sd)
    allocate(bx(n),by(n),bstats(nb)); bstats=missing_value()
    valid=0
    do b=1,nb
      bx=xx(indices(:,b)); by=yy(indices(:,b))
      bd=modified_sharpe_difference(bx,by,level,.false.,kind_test)
      bse=modified_sharpe_standard_error(bx,by,level,.false.,kind_test,bl)
      if (.not. finite_value(bd) .or. .not. finite_value(bse) .or. bse <= tiny(1.0_dp)) cycle
      valid=valid+1
      if (ptype == 1) then
        bstats(valid)=abs(bd-(d+d0))/bse
      else
        bstats(valid)=(bd-(d+d0))/bse
      end if
    end do
    if (valid < 1) then
      result%status=4; result%message='all bootstrap replicates were degenerate'; return
    end if
    tvalue=d/se
    if (ptype == 1) then
      upper=count(bstats(1:valid) >= abs(tvalue))
      pval=real(upper+1,dp)/real(valid+1,dp)
    else
      upper=count(bstats(1:valid) > tvalue)
      lower=count(bstats(1:valid) < tvalue)
      pval=2.0_dp*real(min(upper+1,lower+1),dp)/real(valid+1,dp)
    end if
    allocate(m(n,2)); m(:,1)=xx; m(:,2)=yy
    call modified_sharpe(m,level,vals,counts,reject_negative)
    allocate(result%estimate(1,2),result%difference(1),result%standard_error(1), &
             result%tstat(1),result%pvalue(1))
    result%estimate(1,:)=vals
    result%difference(1)=vals(1)-vals(2)
    result%standard_error(1)=se
    result%tstat(1)=tvalue
    result%pvalue(1)=clamp_probability(pval)
  end subroutine modified_sharpe_testing_bootstrap

  subroutine stationary_indices(n_source, n_output, average_block, indices)
    integer, intent(in) :: n_source, n_output, average_block
    integer, intent(out) :: indices(n_output)
    integer :: current, start, length_block, j
    current=1
    do while (current <= n_output)
      start=random_integer(n_source)
      length_block=1+random_geometric(1.0_dp/real(max(1,average_block),dp))
      do j=0,min(length_block,n_output-current+1)-1
        indices(current+j)=1+modulo(start-1+j,n_source)
      end do
      current=current+min(length_block,n_output-current+1)
    end do
  end subroutine stationary_indices

  subroutine fit_var1(x, y, coefficients, residuals, ok)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(out) :: coefficients(3,2)
    real(dp), allocatable, intent(out) :: residuals(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: design(:,:), beta(:), resid(:), cov(:,:), se(:), ts(:), pv(:)
    logical :: ok1, ok2
    integer :: n
    n=size(x)
    ok=size(y)==n .and. n>4
    if (.not. ok) then
      allocate(residuals(0,0)); return
    end if
    allocate(design(n-1,3),beta(3),resid(n-1),cov(3,3),se(3),ts(3),pv(3),residuals(n-1,2))
    design(:,1)=1.0_dp; design(:,2)=x(1:n-1); design(:,3)=y(1:n-1)
    call ols_fit(design,x(2:n),beta,resid,cov,se,ts,pv,.false.,ok1)
    if (ok1) then
      coefficients(:,1)=beta; residuals(:,1)=resid
    end if
    call ols_fit(design,y(2:n),beta,resid,cov,se,ts,pv,.false.,ok2)
    if (ok2) then
      coefficients(:,2)=beta; residuals(:,2)=resid
    end if
    ok=ok1 .and. ok2
  end subroutine fit_var1

  integer function sharpe_block_size(x, y, candidates, alpha, m_boot, k_sim, &
                                     average_block, start_length, ttype, seed) result(best_block)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in), optional :: candidates(:), m_boot, k_sim, average_block, start_length, ttype, seed
    real(dp), intent(in), optional :: alpha
    integer, parameter :: default_candidates(4)=[1,3,6,10]
    integer, allocatable :: cand(:), ids(:)
    real(dp), allocatable :: residuals(:,:), simulated(:,:)
    real(dp), allocatable :: rejection(:)
    real(dp) :: coefficients(3,2), target, sig_level
    type(test_result) :: test = test_result()
    logical :: ok
    integer :: mb, ks, bav, warm, kind_test, sd, n, k, j, t, best
    n=size(x); best_block=1
    if (size(y)/=n .or. n<5) return
    if (present(candidates)) then
      allocate(cand(size(candidates)))
      cand = candidates
    else
      allocate(cand(size(default_candidates)))
      cand = default_candidates
    end if
    mb=199; if (present(m_boot)) mb=m_boot
    ks=500; if (present(k_sim)) ks=k_sim
    bav=5; if (present(average_block)) bav=average_block
    warm=50; if (present(start_length)) warm=start_length
    kind_test=2; if (present(ttype)) kind_test=ttype
    sd=12345; if (present(seed)) sd=seed
    sig_level=0.05_dp; if (present(alpha)) sig_level=alpha
    target=sharpe_difference(x,y,kind_test)
    call fit_var1(x,y,coefficients,residuals,ok)
    if (.not. ok) return
    allocate(rejection(size(cand)),simulated(warm+n,2),ids(warm+n-1)); rejection=0.0_dp
    call set_random_seed(sd)
    do k=1,ks
      call stationary_indices(n-1,warm+n-1,bav,ids)
      simulated=0.0_dp; simulated(1,:)=[x(1),y(1)]
      do t=2,warm+n
        simulated(t,1)=coefficients(1,1)+coefficients(2,1)*simulated(t-1,1)+ &
                       coefficients(3,1)*simulated(t-1,2)+residuals(ids(t-1),1)
        simulated(t,2)=coefficients(1,2)+coefficients(2,2)*simulated(t-1,1)+ &
                       coefficients(3,2)*simulated(t-1,2)+residuals(ids(t-1),2)
      end do
      do j=1,size(cand)
        if (cand(j)>n) cycle
        call sharpe_testing_bootstrap(simulated(warm+1:,1),simulated(warm+1:,2),test, &
             mb,cand(j),kind_test,1,sd+k*100+j,5,target)
        if (test%status == 0) then
          if (allocated(test%pvalue)) then
            if (test%pvalue(1) <= sig_level) rejection(j) = rejection(j)+1.0_dp
          end if
        end if
      end do
    end do
    rejection=rejection/real(max(1,ks),dp)
    best=1
    do j=2,size(cand)
      if (abs(rejection(j)-sig_level)<abs(rejection(best)-sig_level)) best=j
    end do
    best_block=cand(best)
  end function sharpe_block_size

  integer function modified_sharpe_block_size(x, y, level, candidates, alpha, m_boot, &
                                              k_sim, average_block, start_length, &
                                              ttype, seed, na_negative) result(best_block)
    real(dp), intent(in) :: x(:), y(:), level
    integer, intent(in), optional :: candidates(:), m_boot, k_sim, average_block, start_length, ttype, seed
    real(dp), intent(in), optional :: alpha
    logical, intent(in), optional :: na_negative
    integer, parameter :: default_candidates(4)=[1,3,6,10]
    integer, allocatable :: cand(:), ids(:)
    real(dp), allocatable :: residuals(:,:), simulated(:,:), rejection(:)
    real(dp) :: coefficients(3,2), target, sig_level
    type(test_result) :: test = test_result()
    logical :: ok, reject_negative
    integer :: mb, ks, bav, warm, kind_test, sd, n, k, j, t, best
    n=size(x); best_block=1
    if (size(y)/=n .or. n<5) return
    if (present(candidates)) then
      allocate(cand(size(candidates)))
      cand = candidates
    else
      allocate(cand(size(default_candidates)))
      cand = default_candidates
    end if
    mb=199; if (present(m_boot)) mb=m_boot
    ks=500; if (present(k_sim)) ks=k_sim
    bav=5; if (present(average_block)) bav=average_block
    warm=50; if (present(start_length)) warm=start_length
    kind_test=2; if (present(ttype)) kind_test=ttype
    sd=12345; if (present(seed)) sd=seed
    sig_level=0.05_dp; if (present(alpha)) sig_level=alpha
    reject_negative=.true.; if (present(na_negative)) reject_negative=na_negative
    target=modified_sharpe_difference(x,y,level,reject_negative,kind_test)
    call fit_var1(x,y,coefficients,residuals,ok)
    if (.not. ok .or. .not. finite_value(target)) return
    allocate(rejection(size(cand)),simulated(warm+n,2),ids(warm+n-1)); rejection=0.0_dp
    call set_random_seed(sd)
    do k=1,ks
      call stationary_indices(n-1,warm+n-1,bav,ids)
      simulated=0.0_dp; simulated(1,:)=[x(1),y(1)]
      do t=2,warm+n
        simulated(t,1)=coefficients(1,1)+coefficients(2,1)*simulated(t-1,1)+ &
                       coefficients(3,1)*simulated(t-1,2)+residuals(ids(t-1),1)
        simulated(t,2)=coefficients(1,2)+coefficients(2,2)*simulated(t-1,1)+ &
                       coefficients(3,2)*simulated(t-1,2)+residuals(ids(t-1),2)
      end do
      do j=1,size(cand)
        if (cand(j)>n) cycle
        call modified_sharpe_testing_bootstrap(simulated(warm+1:,1),simulated(warm+1:,2), &
             level,test,reject_negative,mb,cand(j),kind_test,1,sd+k*100+j,5,target)
        if (test%status == 0) then
          if (allocated(test%pvalue)) then
            if (test%pvalue(1) <= sig_level) rejection(j) = rejection(j)+1.0_dp
          end if
        end if
      end do
    end do
    rejection=rejection/real(max(1,ks),dp)
    best=1
    do j=2,size(cand)
      if (abs(rejection(j)-sig_level)<abs(rejection(best)-sig_level)) best=j
    end do
    best_block=cand(best)
  end function modified_sharpe_block_size

end module peerperformance_bootstrap
