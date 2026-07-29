! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_bounds
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use qrmtools_kinds, only : dp
  use qrmtools_types, only : rearrangement_result, ra_bounds_result
  use qrmtools_stats, only : stable_order, variance_value, random_uniform
  use qrmtools_risk, only : es_np
  implicit none
  private
  integer, parameter, public :: bound_worst_var=1, bound_best_var=2, bound_best_es=3
  public :: indices_opp_ordered_to, num_of_opp_ordered_cols
  public :: rearrange_matrix, block_rearrange_matrix
  public :: ra_bounds, adaptive_ra_bounds
  public :: crude_var_bounds_hom, pareto_var_bounds_hom, wang_pareto_bounds
  public :: dual_bound_value, dual_worst_var
  abstract interface
    function scalar_callback(x) result(value)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: value
    end function scalar_callback
    function quantile_callback(p, margin) result(value)
      import dp
      real(dp), intent(in) :: p
      integer, intent(in) :: margin
      real(dp) :: value
    end function quantile_callback
  end interface
contains
  function indices_opp_ordered_to(x) result(indices)
    real(dp), intent(in) :: x(:)
    integer, allocatable :: indices(:),order_desc(:)
    integer :: i
    order_desc=stable_order(x,.true.); allocate(indices(size(x)))
    do i=1,size(x); indices(order_desc(i))=i; end do
  end function indices_opp_ordered_to

  integer function num_of_opp_ordered_cols(x,tolerance) result(number)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in), optional :: tolerance
    real(dp) :: tol
    real(dp), allocatable :: rows(:),sorted(:)
    integer, allocatable :: idx(:)
    integer :: j
    tol=0.0_dp; if(present(tolerance))tol=tolerance; rows=sum(x,dim=2); number=0
    do j=1,size(x,2)
      sorted=x(:,j); call sort_local(sorted); idx=indices_opp_ordered_to(rows-x(:,j))
      if(maxval(abs(sorted(idx)-x(:,j)))<=tol)number=number+1
    end do
  end function num_of_opp_ordered_cols

  function rearrange_matrix(input,method,level,tolerance,relative_tolerance,n_lookback,&
      max_iterations,sample_columns,already_sorted) result(output)
    real(dp), intent(in) :: input(:,:)
    integer, intent(in), optional :: method
    real(dp), intent(in), optional :: level,tolerance
    logical, intent(in), optional :: relative_tolerance,sample_columns,already_sorted
    integer, intent(in), optional :: n_lookback,max_iterations
    type(rearrangement_result) :: output
    real(dp), allocatable :: x(:,:),sorted(:,:),row_sums(:),new_column(:),old_opt(:),history(:)
    integer, allocatable :: idx(:),perm(:)
    integer :: n,d,j,iter,meth,lookback,maxit,i,nopt
    real(dp) :: tol,current,previous,delta,lev
    logical :: relative,sample,sorted_input,changed
    n=size(input,1); d=size(input,2); meth=bound_worst_var; if(present(method))meth=method
    tol=0.0_dp; if(present(tolerance))tol=tolerance; relative=.false.
    if(present(relative_tolerance))relative=relative_tolerance
    lookback=d; if(present(n_lookback))lookback=n_lookback; maxit=max(64,20*d)
    if(present(max_iterations))maxit=max_iterations; sample=.false.; if(present(sample_columns))sample=sample_columns
    sorted_input=.false.; if(present(already_sorted))sorted_input=already_sorted
    lev=0.95_dp; if(present(level))lev=level
    if(n<2 .or. d<2 .or. lookback<1 .or. maxit<1) then; output%message='Invalid rearrangement dimensions.'; return; end if
    x=input; allocate(sorted(n,d),row_sums(n),new_column(n),old_opt(lookback),history(maxit))
    do j=1,d
      sorted(:,j)=input(:,j); if(.not.sorted_input)call sort_local(sorted(:,j))
      if(sample) then
        call random_permutation(n,perm); x(:,j)=x(perm,j)
      end if
    end do
    row_sums=sum(x,dim=2); old_opt=huge(1.0_dp); j=0; changed=.true.
    do iter=1,maxit
      j=mod(j,d)+1; idx=indices_opp_ordered_to(row_sums-x(:,j)); new_column=sorted(idx,j)
      changed=maxval(abs(new_column-x(:,j)))>0.0_dp
      row_sums=row_sums-x(:,j)+new_column; x(:,j)=new_column
      current=objective(row_sums,meth,lev); history(iter)=current
      i=mod(iter-1,lookback)+1; previous=old_opt(i); old_opt(i)=current
      if(iter>lookback) then
        if(relative) then; delta=abs(current-previous)/max(abs(previous),tiny(1.0_dp))
        else; delta=abs(current-previous); end if
        if(delta<=tol) then; output%converged=.true.; exit; end if
      else if(.not.changed .and. iter>=d) then
        output%converged=.true.; delta=0.0_dp; exit
      end if
      delta=huge(1.0_dp)
    end do
    output%iterations=min(iter,maxit); output%bound=objective(row_sums,meth,lev); output%tolerance=delta
    output%rearranged=x; output%optimal_values=history(:output%iterations)
    current=output%bound; nopt=count(abs(row_sums-current)<=1.0e-12_dp*max(1.0_dp,abs(current)))
    allocate(output%optimal_row(d),source=0.0_dp)
    if(nopt>0) then
      do i=1,n
        if(abs(row_sums(i)-current)<=1.0e-12_dp*max(1.0_dp,abs(current))) then
          output%optimal_row=output%optimal_row+x(i,:)
        end if
      end do
      output%optimal_row=output%optimal_row/real(nopt,dp)
    end if
    output%ok=.true.
  end function rearrange_matrix

  function block_rearrange_matrix(input,method,level,tolerance,relative_tolerance,n_lookback,&
      max_iterations,sample_columns) result(output)
    real(dp), intent(in) :: input(:,:)
    integer, intent(in), optional :: method
    real(dp), intent(in), optional :: level,tolerance
    logical, intent(in), optional :: relative_tolerance,sample_columns
    integer, intent(in), optional :: n_lookback,max_iterations
    type(rearrangement_result) :: output
    real(dp), allocatable :: x(:,:),row_sums(:),block_sum(:),complement(:),ordered_sum(:),history(:),old_var(:)
    real(dp), allocatable :: block_matrix(:,:)
    integer, allocatable :: block(:),ord(:),opp(:),perm(:)
    integer :: n,d,meth,lookback,maxit,iter,bsize,i,j,nopt
    real(dp) :: tol,lev,current,var_current,var_old,delta
    logical :: relative,sample
    n=size(input,1); d=size(input,2); meth=bound_worst_var; if(present(method))meth=method
    tol=0.0_dp; if(present(tolerance))tol=tolerance; relative=.false.
    if(present(relative_tolerance))relative=relative_tolerance
    lookback=d; if(present(n_lookback))lookback=n_lookback; maxit=max(128,40*d)
    if(present(max_iterations))maxit=max_iterations; sample=.false.; if(present(sample_columns))sample=sample_columns
    lev=0.95_dp; if(present(level))lev=level
    x=input
    if(sample) then; do j=1,d; call random_permutation(n,perm); x(:,j)=x(perm,j); end do; end if
    allocate(row_sums(n),block_sum(n),complement(n),ordered_sum(n),history(maxit),old_var(lookback))
    row_sums=sum(x,dim=2); old_var=huge(1.0_dp)
    do iter=1,maxit
      bsize=1+int(random_uniform()*real(d-1,dp)); if(bsize>d/2)bsize=d-bsize
      call sample_indices(d,bsize,block); block_sum=0.0_dp
      do i=1,bsize; block_sum=block_sum+x(:,block(i)); end do
      complement=row_sums-block_sum; ord=stable_order(block_sum); opp=indices_opp_ordered_to(complement)
      allocate(block_matrix(n,bsize))
      block_matrix=x(ord,block)
      block_matrix=block_matrix(opp,:)
      x(:,block)=block_matrix
      deallocate(block_matrix)
      ordered_sum=block_sum(ord); ordered_sum=ordered_sum(opp); row_sums=complement+ordered_sum
      current=objective(row_sums,meth,lev); history(iter)=current; var_current=variance_value(row_sums)
      i=mod(iter-1,lookback)+1; var_old=old_var(i); old_var(i)=var_current
      if(iter>lookback) then
        if(relative) then; delta=abs(var_current-var_old)/max(abs(var_old),tiny(1.0_dp))
        else; delta=abs(var_current-var_old); end if
        if(delta<=tol) then; output%converged=.true.; exit; end if
      end if
      delta=huge(1.0_dp)
    end do
    output%iterations=min(iter,maxit); output%bound=objective(row_sums,meth,lev); output%tolerance=delta
    output%rearranged=x; output%optimal_values=history(:output%iterations); current=output%bound
    nopt=count(abs(row_sums-current)<=1.0e-12_dp*max(1.0_dp,abs(current)))
    allocate(output%optimal_row(d),source=0.0_dp)
    if(nopt>0) then
      do i=1,n
        if(abs(row_sums(i)-current)<=1.0e-12_dp*max(1.0_dp,abs(current))) then
          output%optimal_row=output%optimal_row+x(i,:)
        end if
      end do
      output%optimal_row=output%optimal_row/real(nopt,dp)
    end if
    output%ok=.true.
  end function block_rearrange_matrix


  function ra_bounds(level, dimension, n_grid, quantile, method, tolerance, &
      relative_tolerance, n_lookback, max_iterations, sample_columns, &
      use_blocks) result(output)
    real(dp), intent(in) :: level
    integer, intent(in) :: dimension
    integer, intent(in) :: n_grid
    procedure(quantile_callback) :: quantile
    integer, intent(in), optional :: method
    real(dp), intent(in), optional :: tolerance
    logical, intent(in), optional :: relative_tolerance
    integer, intent(in), optional :: n_lookback
    integer, intent(in), optional :: max_iterations
    logical, intent(in), optional :: sample_columns
    logical, intent(in), optional :: use_blocks
    type(ra_bounds_result) :: output
    real(dp), allocatable :: lower_matrix(:,:)
    real(dp), allocatable :: upper_matrix(:,:)
    real(dp) :: probability
    real(dp) :: replacement_probability
    real(dp) :: tol
    integer :: meth
    integer :: lookback
    integer :: maxit
    integer :: i
    integer :: j
    logical :: relative
    logical :: sample
    logical :: blocks

    if(level <= 0.0_dp .or. level >= 1.0_dp .or. dimension < 2 .or. n_grid < 2) then
      output%message = 'Invalid RA dimensions or confidence level.'
      return
    end if

    meth = bound_worst_var
    if(present(method)) meth = method
    tol = 0.0_dp
    if(present(tolerance)) tol = tolerance
    relative = .false.
    if(present(relative_tolerance)) relative = relative_tolerance
    lookback = dimension
    if(present(n_lookback)) lookback = n_lookback
    maxit = max(64, 20*dimension)
    if(present(max_iterations)) maxit = max_iterations
    sample = .true.
    if(present(sample_columns)) sample = sample_columns
    blocks = .false.
    if(present(use_blocks)) blocks = use_blocks

    allocate(lower_matrix(n_grid,dimension))
    allocate(upper_matrix(n_grid,dimension))
    do i = 1, n_grid
      select case(meth)
      case(bound_worst_var)
        probability = level + (1.0_dp-level)*real(i-1,dp)/real(n_grid,dp)
      case(bound_best_var)
        probability = level*real(i-1,dp)/real(n_grid,dp)
      case(bound_best_es)
        probability = real(i-1,dp)/real(n_grid,dp)
      case default
        output%message = 'Unknown rearrangement method.'
        return
      end select
      do j = 1, dimension
        lower_matrix(i,j) = quantile(probability,j)
      end do

      select case(meth)
      case(bound_worst_var)
        probability = level + (1.0_dp-level)*real(i,dp)/real(n_grid,dp)
      case(bound_best_var)
        probability = level*real(i,dp)/real(n_grid,dp)
      case(bound_best_es)
        probability = real(i,dp)/real(n_grid,dp)
      end select
      do j = 1, dimension
        upper_matrix(i,j) = quantile(probability,j)
      end do
    end do

    if(meth == bound_best_var .or. meth == bound_best_es) then
      replacement_probability = level/(2.0_dp*real(n_grid,dp))
      do j = 1, dimension
        if(.not. ieee_is_finite(lower_matrix(1,j))) then
          lower_matrix(1,j) = quantile(replacement_probability,j)
        end if
      end do
    end if

    if(meth == bound_worst_var .or. meth == bound_best_es) then
      replacement_probability = level + (1.0_dp-level) * &
        (1.0_dp-1.0_dp/(2.0_dp*real(n_grid,dp)))
      if(meth == bound_best_es) then
        replacement_probability = 1.0_dp-1.0_dp/(2.0_dp*real(n_grid,dp))
      end if
      do j = 1, dimension
        if(.not. ieee_is_finite(upper_matrix(n_grid,j))) then
          upper_matrix(n_grid,j) = quantile(replacement_probability,j)
        end if
      end do
    end if

    if(blocks) then
      output%lower = block_rearrange_matrix(lower_matrix,meth,level,tol, &
        relative,lookback,maxit,sample)
      output%upper = block_rearrange_matrix(upper_matrix,meth,level,tol, &
        relative,lookback,maxit,sample)
    else
      output%lower = rearrange_matrix(lower_matrix,meth,level,tol,relative, &
        lookback,maxit,sample,.true.)
      output%upper = rearrange_matrix(upper_matrix,meth,level,tol,relative, &
        lookback,maxit,sample,.true.)
    end if

    if(.not. output%lower%ok .or. .not. output%upper%ok) then
      output%message = 'A rearrangement calculation failed.'
      return
    end if

    output%bounds = [output%lower%bound, output%upper%bound]
    output%relative_gap = abs(output%bounds(2)-output%bounds(1)) / &
      max(abs(output%bounds(2)),tiny(1.0_dp))
    output%tolerances(1:2) = [output%lower%tolerance, output%upper%tolerance]
    output%tolerances(3) = output%relative_gap
    output%converged(1:2) = [output%lower%converged, output%upper%converged]
    output%converged(3) = .true.
    output%n_used = n_grid
    output%ok = .true.
  end function ra_bounds

  function adaptive_ra_bounds(level, dimension, exponents, quantile, method, &
      joint_tolerance, individual_tolerance, n_lookback, max_iterations, &
      sample_columns, use_blocks) result(output)
    real(dp), intent(in) :: level
    integer, intent(in) :: dimension
    integer, intent(in) :: exponents(:)
    procedure(quantile_callback) :: quantile
    integer, intent(in), optional :: method
    real(dp), intent(in), optional :: joint_tolerance
    real(dp), intent(in), optional :: individual_tolerance
    integer, intent(in), optional :: n_lookback
    integer, intent(in), optional :: max_iterations
    logical, intent(in), optional :: sample_columns
    logical, intent(in), optional :: use_blocks
    type(ra_bounds_result) :: output
    real(dp) :: joint_tol
    real(dp) :: individual_tol
    integer :: meth
    integer :: lookback
    integer :: maxit
    integer :: k
    integer :: n_grid
    logical :: sample
    logical :: blocks

    if(size(exponents) < 1 .or. any(exponents < 1)) then
      output%message = 'At least one positive discretization exponent is required.'
      return
    end if

    meth = bound_worst_var
    if(present(method)) meth = method
    joint_tol = 0.01_dp
    if(present(joint_tolerance)) joint_tol = joint_tolerance
    individual_tol = 0.0_dp
    if(present(individual_tolerance)) individual_tol = individual_tolerance
    lookback = dimension
    if(present(n_lookback)) lookback = n_lookback
    maxit = max(64,20*dimension)
    if(present(max_iterations)) maxit = max_iterations
    sample = .true.
    if(present(sample_columns)) sample = sample_columns
    blocks = .false.
    if(present(use_blocks)) blocks = use_blocks

    do k = 1, size(exponents)
      n_grid = 2**exponents(k)
      output = ra_bounds(level,dimension,n_grid,quantile,meth, &
        individual_tol,.not.blocks,lookback,maxit,sample,blocks)
      if(.not.output%ok) return
      output%converged(3) = output%relative_gap <= joint_tol
      if(output%lower%converged .and. output%upper%converged .and. &
          output%converged(3)) exit
    end do
  end function adaptive_ra_bounds

  function crude_var_bounds_hom(level,dimension,quantile) result(bounds)
    real(dp), intent(in) :: level
    integer, intent(in) :: dimension
    procedure(scalar_callback) :: quantile
    real(dp) :: bounds(2)
    bounds=real(dimension,dp)*[quantile(level/real(dimension,dp)),&
      quantile((real(dimension-1,dp)+level)/real(dimension,dp))]
  end function crude_var_bounds_hom

  function pareto_var_bounds_hom(level,dimension,shape) result(bounds)
    real(dp), intent(in) :: level,shape
    integer, intent(in) :: dimension
    real(dp) :: bounds(2)
    if(dimension==2) then
      bounds=[(1.0_dp-level)**(-1.0_dp/shape)-1.0_dp,&
        2.0_dp*((0.5_dp*(1.0_dp-level))**(-1.0_dp/shape)-1.0_dp)]
    else
      bounds=wang_pareto_bounds(level,dimension,shape)
    end if
  end function pareto_var_bounds_hom

  function wang_pareto_bounds(level,dimension,shape) result(bounds)
    real(dp), intent(in) :: level,shape
    integer, intent(in) :: dimension
    real(dp) :: bounds(2),ibar,lo,hi,mid,flo,fmid,c,t1
    integer :: i
    if(abs(shape-1.0_dp)<1.0e-12_dp) then
      ibar=-log(1.0_dp-level)-level
    else
      ibar=((1.0_dp-level)**(1.0_dp-1.0_dp/shape)-1.0_dp)/(1.0_dp-1.0_dp/shape)-level
    end if
    bounds(1)=max((1.0_dp-level)**(-1.0_dp/shape)-1.0_dp,real(dimension,dp)*ibar)
    lo=max(tiny(1.0_dp),(1.0_dp-level)*1.0e-10_dp); hi=(1.0_dp-level)/real(dimension,dp)*(1.0_dp-1.0e-10_dp)
    flo=wang_h_pareto(lo,level,dimension,shape)
    do i=1,200
      mid=0.5_dp*(lo+hi); fmid=wang_h_pareto(mid,level,dimension,shape)
      if(flo*fmid<=0.0_dp) then; hi=mid; else; lo=mid; flo=fmid; end if
    end do
    c=0.5_dp*(lo+hi); t1=(1.0_dp-level)/c-real(dimension-1,dp)
    bounds(2)=c**(-1.0_dp/shape)*((real(dimension-1,dp))*t1**(-1.0_dp/shape)+1.0_dp)-real(dimension,dp)
  end function wang_pareto_bounds

  real(dp) function dual_bound_value(s,dimension,cdf) result(value)
    real(dp), intent(in) :: s
    integer, intent(in) :: dimension
    procedure(scalar_callback) :: cdf
    real(dp) :: lo,hi,mid,flo,fmid
    integer :: i
    if(s<=0.0_dp) then; value=real(dimension,dp); return; end if
    lo=0.0_dp; hi=s/real(dimension,dp)*(1.0_dp-1.0e-10_dp); flo=dual_foc(s,lo,dimension,cdf)
    do i=1,150
      mid=0.5_dp*(lo+hi); fmid=dual_foc(s,mid,dimension,cdf)
      if(flo*fmid<=0.0_dp) then; hi=mid; else; lo=mid; flo=fmid; end if
    end do
    mid=0.5_dp*(lo+hi)
    value=1.0_dp-cdf(mid)+real(dimension-1,dp)*(1.0_dp-cdf(s-real(dimension-1,dp)*mid))
  end function dual_bound_value

  real(dp) function dual_worst_var(level,dimension,cdf,lower,upper) result(value)
    real(dp), intent(in) :: level,lower,upper
    integer, intent(in) :: dimension
    procedure(scalar_callback) :: cdf
    real(dp) :: lo,hi,mid,flo,fmid
    integer :: i
    lo=lower; hi=upper; flo=dual_bound_value(lo,dimension,cdf)-(1.0_dp-level)
    do i=1,160
      mid=0.5_dp*(lo+hi); fmid=dual_bound_value(mid,dimension,cdf)-(1.0_dp-level)
      if(flo*fmid<=0.0_dp) then; hi=mid; else; lo=mid; flo=fmid; end if
    end do
    value=0.5_dp*(lo+hi)
  end function dual_worst_var

  real(dp) function dual_foc(s,t,dimension,cdf) result(value)
    real(dp), intent(in) :: s,t
    integer, intent(in) :: dimension
    procedure(scalar_callback) :: cdf
    real(dp) :: dval
    dval=dual_d(s,t,dimension,cdf)
    value=dval-(1.0_dp-cdf(t)+real(dimension-1,dp)*(1.0_dp-cdf(s-real(dimension-1,dp)*t)))
  end function dual_foc

  real(dp) function dual_d(s,t,dimension,cdf) result(value)
    real(dp), intent(in) :: s,t
    integer, intent(in) :: dimension
    procedure(scalar_callback) :: cdf
    real(dp) :: a,b,h,sumv,x
    integer, parameter :: nint=2048
    integer :: i
    a=t; b=s-real(dimension-1,dp)*t
    if(abs(s-real(dimension,dp)*t)<=1.0e-12_dp) then
      value=real(dimension,dp)*(1.0_dp-cdf(s/real(dimension,dp))); return
    end if
    h=(b-a)/real(nint,dp); sumv=0.5_dp*(cdf(a)+cdf(b))
    do i=1,nint-1; x=a+h*real(i,dp); sumv=sumv+cdf(x); end do
    value=real(dimension,dp)*(1.0_dp-h*sumv/(s-real(dimension,dp)*t))
  end function dual_d

  pure real(dp) function wang_h_pareto(c,level,dimension,shape) result(value)
    real(dp), intent(in) :: c,level,shape
    integer, intent(in) :: dimension
    real(dp) :: t1,t2,ibar,aux
    t1=(1.0_dp-level)/c-real(dimension-1,dp); t2=1.0_dp-level-real(dimension,dp)*c
    if(abs(shape-1.0_dp)<1.0e-12_dp) then; ibar=log(t1)/t2-1.0_dp
    else; ibar=shape/(1.0_dp-shape)*c**(1.0_dp-1.0_dp/shape)*(1.0_dp-t1**(1.0_dp-1.0_dp/shape))/t2-1.0_dp; end if
    aux=c**(-1.0_dp/shape)/real(dimension,dp)*(real(dimension-1,dp)*t1**(-1.0_dp/shape)+1.0_dp)-1.0_dp
    value=ibar-aux
  end function wang_h_pareto

  real(dp) function objective(x,method,level) result(value)
    real(dp), intent(in) :: x(:),level
    integer, intent(in) :: method
    select case(method)
    case(bound_worst_var); value=minval(x)
    case(bound_best_var); value=maxval(x)
    case(bound_best_es); value=es_np(x,level)
    case default; value=minval(x)
    end select
  end function objective

  subroutine sort_local(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: key
    do i=2,size(x); key=x(i); j=i-1; do while(j>=1); if(x(j)<=key)exit; x(j+1)=x(j); j=j-1; end do; x(j+1)=key; end do
  end subroutine sort_local

  subroutine random_permutation(n,perm)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: perm(:)
    integer :: i,j,tmp
    allocate(perm(n)); perm=[(i,i=1,n)]
    do i=n,2,-1; j=1+int(random_uniform()*real(i,dp)); if(j>i)j=i; tmp=perm(i); perm(i)=perm(j); perm(j)=tmp; end do
  end subroutine random_permutation

  subroutine sample_indices(n,k,indices)
    integer, intent(in) :: n,k
    integer, allocatable, intent(out) :: indices(:)
    integer, allocatable :: perm(:)
    call random_permutation(n,perm); indices=perm(1:k)
  end subroutine sample_indices
end module qrmtools_bounds
