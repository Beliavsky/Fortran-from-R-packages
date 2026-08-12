module flsss_api
  use flsss_kinds, only : dp, i8
  use flsss_types
  use flsss_search, only : search_1d, search_md, search_md_i8
  use flsss_fast_search, only : search_md_i8_packed, search_md_i8_pat
  use flsss_mpat, only : search_md_i8_mpat
  use flsss_packed, only : zero_minimum_i8, sort_rows_for_comonotonicity, is_comonotonic_i8
  use flsss_parallel, only : openmp_enabled
  use flsss_util, only : append_solution, timer_type
  use flsss_integerize, only : integerize_problem
  use flsss_knapsack, only : mm_knapsack, mm_knapsack_integerized, aux_knapsack01bb, aux_knapsack01dp
  use flsss_gap, only : gap_solve, aux_gap_bb, aux_gap_bbdp, aux_gap_ga
  use flsss_arbitrary, only : arb_flsss, decompose_arb_flsss, arb_flsss_obj_run, &
    arb_flsss_decomp_run, build_ksum_hash
  use flsss_bigint, only : add_num_strings
  implicit none
  private
  public :: dp, i8
  public :: subset_solution, subset_solutions, real_bucket, multiset_solution, multiset_solutions
  public :: integerized_problem, integerized_search_result, knapsack_result, knapsack_multi_result, gap_result
  public :: mflsss_object, mflsss_decomposition, arb_flsss_object, arb_flsss_decomposition, ksum_table
  public :: flsss, flsss_multiset, mflsss_par, mflsss_par_impose_bounds
  public :: mflsss_par_integerized, mflsss_par_integerized_parallel, mflsss_par_impose_bounds_integerized
  public :: decompose_mflsss, mflsss_obj_run, mflsss_decomp_run
  public :: mm_knapsack, mm_knapsack_integerized, aux_knapsack01bb, aux_knapsack01dp
  public :: gap_solve, aux_gap_bb, aux_gap_bbdp, aux_gap_ga
  public :: arb_flsss, decompose_arb_flsss, arb_flsss_obj_run, arb_flsss_decomp_run, build_ksum_hash, add_num_strings

contains

  function flsss(len, v, target, me, solution_need, lb, ub, via_conjugate, tlimit, n_fraction_digits) result(r)
    integer, intent(in) :: len
    real(dp), intent(in) :: v(:), target, me
    integer, intent(in), optional :: solution_need, lb(:), ub(:), n_fraction_digits
    logical, intent(in), optional :: via_conjugate
    real(dp), intent(in), optional :: tlimit
    type(subset_solutions) :: r
    type(subset_solutions) :: tmp
    real(dp), allocatable :: vv(:)
    real(dp) :: tt, ee, lim, scale
    integer :: need, k, remain, j, q
    logical :: use_conj
    integer, allocatable :: mark(:), comp(:)

    need=1; if(present(solution_need)) need=max(1,solution_need)
    lim=huge(1.0_dp); if(present(tlimit)) lim=tlimit
    allocate(vv(size(v))); vv=v; tt=target; ee=me
    use_conj=.false.; if(present(via_conjugate)) use_conj=via_conjugate
    if(present(n_fraction_digits)) then
      scale=10.0_dp**n_fraction_digits
      vv=real(nint(vv*scale,kind=i8),dp)
      tt=real(nint(tt*scale,kind=i8),dp); ee=real(nint(ee*scale,kind=i8),dp)
    end if
    allocate(r%sol(0),tmp%sol(0))
    if(len==0) then
      remain=need
      do k=1,size(v)
        call search_1d(vv,k,tt,ee,tmp,remain,lim)
        call merge_results(r,tmp,need)
        if(r%timed_out .or. size(r%sol)>=need) exit
        remain=need-size(r%sol)
      end do
    else
      if (use_conj .and. .not.present(lb) .and. .not.present(ub)) then
        call search_1d(vv,size(v)-len,sum(vv)-tt,ee,r,need,lim)
        allocate(mark(size(v)))
        do j=1,r%size()
          mark=0; mark(r%sol(j)%idx)=1
          allocate(comp(len)); q=0
          do k=1,size(v)
            if(mark(k)==0) then; q=q+1; comp(q)=k; end if
          end do
          deallocate(r%sol(j)%idx); allocate(r%sol(j)%idx(len)); r%sol(j)%idx=comp
          deallocate(comp)
        end do
      else if(present(lb) .and. present(ub)) then
        call search_1d(vv,len,tt,ee,r,need,lim,lb,ub)
      else if(present(lb)) then
        call search_1d(vv,len,tt,ee,r,need,lim,lb=lb)
      else if(present(ub)) then
        call search_1d(vv,len,tt,ee,r,need,lim,ub=ub)
      else
        call search_1d(vv,len,tt,ee,r,need,lim)
      end if
    end if
  end function flsss

  function flsss_multiset(len, buckets, target, me, solution_need, tlimit) result(r)
    integer, intent(in) :: len(:)
    type(real_bucket), intent(in) :: buckets(:)
    real(dp), intent(in) :: target, me
    integer, intent(in), optional :: solution_need
    real(dp), intent(in), optional :: tlimit
    type(multiset_solutions) :: r
    integer :: b,totaln,totalk,off,p,q,nsol
    real(dp), allocatable :: v(:)
    integer, allocatable :: lb(:),ub(:),whichb(:),base(:)
    type(subset_solutions) :: raw
    real(dp)::lim
    integer::need
    if(size(len)/=size(buckets)) error stop "flsss_multiset: len/buckets mismatch"
    totaln=0;totalk=sum(len)
    do b=1,size(buckets); totaln=totaln+size(buckets(b)%value); end do
    allocate(v(totaln),lb(totalk),ub(totalk),whichb(totalk),base(size(buckets)))
    off=0;p=0
    do b=1,size(buckets)
      base(b)=off
      v(off+1:off+size(buckets(b)%value))=buckets(b)%value
      do q=1,len(b)
        p=p+1; lb(p)=off+q; ub(p)=off+size(buckets(b)%value)-len(b)+q; whichb(p)=b
      end do
      off=off+size(buckets(b)%value)
    end do
    need=1;if(present(solution_need))need=max(1,solution_need)
    lim=huge(1.0_dp);if(present(tlimit))lim=tlimit
    call search_1d(v,totalk,target,me,raw,need,lim,lb,ub)
    allocate(r%sol(raw%size())); r%timed_out=raw%timed_out;r%nodes=raw%nodes
    do nsol=1,raw%size()
      allocate(r%sol(nsol)%bucket(size(buckets)))
      do b=1,size(buckets)
        allocate(r%sol(nsol)%bucket(b)%idx(len(b)));p=0
        do q=1,totalk
          if(whichb(q)==b) then;p=p+1;r%sol(nsol)%bucket(b)%idx(p)=raw%sol(nsol)%idx(q)-base(b);end if
        end do
      end do
    end do
  end function flsss_multiset

  function mflsss_par(len,v,target,me,solution_need,tlimit,dl,du) result(r)
    integer,intent(in)::len
    real(dp),intent(in)::v(:,:),target(:),me(:)
    integer,intent(in),optional::solution_need,dl,du
    real(dp),intent(in),optional::tlimit
    type(subset_solutions)::r,tmp
    integer::need,k,remain,nlo,nup
    real(dp)::lim
    need=1;if(present(solution_need))need=max(1,solution_need)
    lim=huge(1.0_dp);if(present(tlimit))lim=tlimit
    nlo=size(v,2);if(present(dl))nlo=dl
    nup=size(v,2);if(present(du))nup=du
    allocate(r%sol(0),tmp%sol(0))
    if(len==0)then
      remain=need
      do k=1,size(v,1)
        call search_md(v,k,target,me,tmp,remain,lim,dl=nlo,du=nup)
        call merge_results(r,tmp,need)
        if(r%timed_out.or.size(r%sol)>=need)exit
        remain=need-size(r%sol)
      end do
    else
      call search_md(v,len,target,me,r,need,lim,dl=nlo,du=nup)
    end if
  end function mflsss_par

  function mflsss_par_impose_bounds(len,v,target,me,lb,ub,solution_need,tlimit,dl,du) result(r)
    integer,intent(in)::len,lb(:),ub(:)
    real(dp),intent(in)::v(:,:),target(:),me(:)
    integer,intent(in),optional::solution_need,dl,du
    real(dp),intent(in),optional::tlimit
    type(subset_solutions)::r
    integer::need,nlo,nup
    real(dp)::lim
    need=1;if(present(solution_need))need=max(1,solution_need)
    lim=huge(1.0_dp);if(present(tlimit))lim=tlimit
    nlo=size(v,2);if(present(dl))nlo=dl;nup=size(v,2);if(present(du))nup=du
    call search_md(v,len,target,me,r,need,lim,lb,ub,nlo,nup)
  end function mflsss_par_impose_bounds

  function mflsss_par_integerized(len,v,target,me,precision_level,solution_need,tlimit,dl,du, &
                                   return_before_mining,engine) result(out)
    integer,intent(in)::len
    real(dp),intent(in)::v(:,:),target(:),me(:)
    integer,intent(in),optional::precision_level(:),solution_need,dl,du
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::return_before_mining
    character(len=*),intent(in),optional::engine
    type(integerized_search_result)::out
    integer(i8),allocatable::wv(:,:),wt(:),sorted(:,:)
    integer,allocatable::ord(:)
    integer::k
    logical::rbm,shiftok,mono
    character(len=16)::mode

    k=len;if(k==0)k=size(v,1)
    out%integerized=integerize_problem(k,v,target,me,precision_level)
    rbm=.false.;if(present(return_before_mining))rbm=return_before_mining
    if(rbm)then;allocate(out%solution%sol(0));return;end if
    if (len <= 0) then
      out%solution=mflsss_par(len,v,target,me,solution_need,tlimit,dl,du)
      return
    end if

    mode='auto'
    if(present(engine)) mode=lower_ascii(trim(engine))
    if(trim(mode)=='dfs') then
      call search_md_i8(out%integerized%v,len,out%integerized%target,out%integerized%me, &
        out%solution,solution_need,tlimit,dl=dl,du=du)
      out%solution%engine='dfs'
      return
    end if

    call zero_minimum_i8(out%integerized%v,len,out%integerized%target,wv,wt,shiftok)
    if(.not.shiftok) then
      call search_md_i8(out%integerized%v,len,out%integerized%target,out%integerized%me, &
        out%solution,solution_need,tlimit,dl=dl,du=du)
      out%solution%engine='dfs-shift-fail'
      return
    end if

    select case(trim(mode))
    case('packed')
      call search_md_i8_packed(wv,len,wt,out%integerized%me,out%solution,solution_need,tlimit,dl=dl,du=du)
    case('pat')
      call sort_rows_for_comonotonicity(wv,ord,sorted,mono)
      if(mono) then
        call search_md_i8_pat(sorted,len,wt,out%integerized%me,out%solution,solution_need,tlimit,dl=dl,du=du)
        call remap_solutions(out%solution,ord)
      else
        call search_md_i8_packed(wv,len,wt,out%integerized%me,out%solution,solution_need,tlimit,dl=dl,du=du)
      end if
    case('mpat')
      call sort_rows_for_comonotonicity(wv,ord,sorted,mono)
      if(mono) then
        call search_md_i8_mpat(sorted,len,wt,out%integerized%me,out%solution,solution_need,tlimit,dl=dl,du=du)
        call remap_solutions(out%solution,ord)
      else
        call search_md_i8_packed(wv,len,wt,out%integerized%me,out%solution,solution_need,tlimit,dl=dl,du=du)
        out%solution%engine='mpat-fallback'
      end if
    case('auto')
      call sort_rows_for_comonotonicity(wv,ord,sorted,mono)
      if(mono .and. mpat_recommended(sorted,len,wt,out%integerized%me,dl,du)) then
        call search_md_i8_mpat(sorted,len,wt,out%integerized%me,out%solution,solution_need,tlimit,dl=dl,du=du)
        call remap_solutions(out%solution,ord)
      else
        call search_md_i8_packed(wv,len,wt,out%integerized%me,out%solution,solution_need,tlimit,dl=dl,du=du)
      end if
    case default
      error stop 'mflsss_par_integerized: engine must be auto, mpat, pat, packed, or dfs'
    end select
  end function mflsss_par_integerized

  function mflsss_par_integerized_parallel(len,v,target,me,max_threads,precision_level,solution_need, &
                                               tlimit,dl,du,return_before_mining) result(out)
    integer,intent(in)::len,max_threads
    real(dp),intent(in)::v(:,:),target(:),me(:)
    integer,intent(in),optional::precision_level(:),solution_need,dl,du
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::return_before_mining
    type(integerized_search_result)::out
    type(subset_solutions),allocatable::part(:)
    integer(i8),allocatable::wv(:,:),wt(:)
    integer::need,nfirst,p,nthreads,k
    real(dp)::lim
    logical::rbm,shiftok

    k=len;if(k==0)k=size(v,1)
    out%integerized=integerize_problem(k,v,target,me,precision_level)
    rbm=.false.;if(present(return_before_mining))rbm=return_before_mining
    if(rbm)then;allocate(out%solution%sol(0));return;end if
    if(len<=0 .or. max_threads<=1) then
      out=mflsss_par_integerized(len,v,target,me,precision_level,solution_need,tlimit,dl,du,.false.,'auto')
      return
    end if
    call zero_minimum_i8(out%integerized%v,len,out%integerized%target,wv,wt,shiftok)
    if(.not.shiftok) then
      out=mflsss_par_integerized(len,v,target,me,precision_level,solution_need,tlimit,dl,du,.false.,'dfs')
      return
    end if
    need=1;if(present(solution_need))need=max(1,solution_need)
    lim=huge(1.0_dp);if(present(tlimit))lim=tlimit
    nthreads=max(1,max_threads)
    nfirst=max(0,size(v,1)-len+1)
    allocate(part(nfirst),out%solution%sol(0))
!$omp parallel do schedule(dynamic) num_threads(nthreads)
    do p=1,nfirst
      call search_md_i8_packed(wv,len,wt,out%integerized%me,part(p),need,lim,dl=dl,du=du,prefix=p)
    end do
!$omp end parallel do
    do p=1,nfirst
      call merge_results(out%solution,part(p),need)
      out%solution%partitions_run=out%solution%partitions_run+1_i8
      if(out%solution%size()>=need)exit
    end do
    if(openmp_enabled()) then
      out%solution%engine='packed-omp'
    else
      out%solution%engine='packed-decomp-serial'
    end if
  end function mflsss_par_integerized_parallel

  function mflsss_par_impose_bounds_integerized(len,v,target,me,lb,ub,precision_level,solution_need, &
                                                tlimit,dl,du,return_before_mining,engine) result(out)
    integer,intent(in)::len,lb(:),ub(:)
    real(dp),intent(in)::v(:,:),target(:),me(:)
    integer,intent(in),optional::precision_level(:),solution_need,dl,du
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::return_before_mining
    character(len=*),intent(in),optional::engine
    type(integerized_search_result)::out
    integer(i8),allocatable::wv(:,:),wt(:)
    logical::rbm,shiftok
    character(len=16)::mode

    out%integerized=integerize_problem(len,v,target,me,precision_level)
    rbm=.false.;if(present(return_before_mining))rbm=return_before_mining
    if(rbm)then;allocate(out%solution%sol(0));return;end if
    mode='auto';if(present(engine))mode=lower_ascii(trim(engine))
    if(trim(mode)=='dfs') then
      call search_md_i8(out%integerized%v,len,out%integerized%target,out%integerized%me, &
        out%solution,solution_need,tlimit,lb,ub,dl,du)
      out%solution%engine='dfs'
      return
    end if
    call zero_minimum_i8(out%integerized%v,len,out%integerized%target,wv,wt,shiftok)
    if(.not.shiftok) then
      call search_md_i8(out%integerized%v,len,out%integerized%target,out%integerized%me, &
        out%solution,solution_need,tlimit,lb,ub,dl,du)
      out%solution%engine='dfs-shift-fail'
      return
    end if
    if(trim(mode)=='pat') then
      call search_md_i8_pat(wv,len,wt,out%integerized%me,out%solution,solution_need,tlimit,lb,ub,dl,du)
    else if(trim(mode)=='mpat') then
      call search_md_i8_mpat(wv,len,wt,out%integerized%me,out%solution,solution_need,tlimit,lb,ub,dl,du)
    else if(trim(mode)=='auto' .and. is_comonotonic_i8(wv)) then
      if(mpat_recommended(wv,len,wt,out%integerized%me,dl,du)) then
        call search_md_i8_mpat(wv,len,wt,out%integerized%me,out%solution,solution_need,tlimit,lb,ub,dl,du)
      else
        call search_md_i8_packed(wv,len,wt,out%integerized%me,out%solution,solution_need,tlimit,lb,ub,dl,du)
      end if
    else if(trim(mode)=='packed' .or. trim(mode)=='auto') then
      call search_md_i8_packed(wv,len,wt,out%integerized%me,out%solution,solution_need,tlimit,lb,ub,dl,du)
    else
      error stop 'mflsss_par_impose_bounds_integerized: invalid engine'
    end if
  end function mflsss_par_impose_bounds_integerized

  function decompose_mflsss(len,v,target,me,dl,du,approx_ninstance) result(dec)
    integer,intent(in)::len
    real(dp),intent(in)::v(:,:),target(:),me(:)
    integer,intent(in),optional::dl,du,approx_ninstance
    type(mflsss_decomposition)::dec
    integer::m,i,d,q,ng,lo,hi
    m=max(0,size(v,1)-len+1);d=size(v,2)
    ng=m
    if (present(approx_ninstance)) then
      if (approx_ninstance < 1) error stop "decompose_mflsss: approx_ninstance must be positive"
      ng=min(m,approx_ninstance)
    end if
    allocate(dec%object(ng),dec%solutions_found%sol(0))
    do i=1,ng
      lo=1+(i-1)*m/max(1,ng)
      hi=i*m/max(1,ng)
      dec%object(i)%len=len;dec%object(i)%v=v;dec%object(i)%target=target;dec%object(i)%me=me
      allocate(dec%object(i)%lb(len),dec%object(i)%ub(len))
      dec%object(i)%lb=[(q,q=1,len)];dec%object(i)%ub=[(size(v,1)-len+q,q=1,len)]
      dec%object(i)%prefix_lo=lo;dec%object(i)%prefix_hi=hi
      if(lo==hi) dec%object(i)%prefix=lo
      dec%object(i)%dl=d;dec%object(i)%du=d
      if(present(dl))dec%object(i)%dl=dl
      if(present(du))dec%object(i)%du=du
    end do
  end function decompose_mflsss

  function mflsss_obj_run(obj,solution_need,tlimit) result(r)
    type(mflsss_object),intent(in)::obj
    integer,intent(in),optional::solution_need
    real(dp),intent(in),optional::tlimit
    type(subset_solutions)::r,tmp
    integer::need,p,plo,phi
    real(dp)::lim,remain
    type(timer_type)::timer
    need=1;if(present(solution_need))need=max(1,solution_need)
    lim=huge(1.0_dp);if(present(tlimit))lim=tlimit
    call timer%start(lim)
    allocate(r%sol(0),tmp%sol(0))
    plo=obj%prefix_lo;phi=obj%prefix_hi
    if(plo<=0 .or. phi<plo) then
      if(obj%prefix>0) then
        plo=obj%prefix;phi=obj%prefix
      else
        plo=1;phi=max(0,size(obj%v,1)-obj%len+1)
      end if
    end if
    do p=plo,phi
      remain=max(0.0_dp,lim-timer%elapsed())
      if(remain<=0.0_dp) then
        r%timed_out=.true.;exit
      end if
      call search_md(obj%v,obj%len,obj%target,obj%me,tmp,need-r%size(),remain, &
        obj%lb,obj%ub,obj%dl,obj%du,p)
      call merge_results(r,tmp,need)
      if(r%timed_out .or. r%size()>=need) exit
    end do
  end function mflsss_obj_run

  function mflsss_decomp_run(dec,solution_need,tlimit,parallel,max_threads) result(r)
    type(mflsss_decomposition),intent(in)::dec
    integer,intent(in),optional::solution_need,max_threads
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::parallel
    type(subset_solutions)::r
    type(subset_solutions),allocatable::part(:)
    integer::i,need,nthreads
    real(dp)::lim,remain
    logical::usepar
    type(timer_type)::timer
    need=1;if(present(solution_need))need=max(1,solution_need)
    lim=huge(1.0_dp);if(present(tlimit))lim=tlimit
    usepar=.false.;if(present(parallel))usepar=parallel
    nthreads=max(1,size(dec%object));if(present(max_threads))nthreads=max(1,max_threads)
    allocate(part(size(dec%object)),r%sol(0))
    if(usepar) then
!$omp parallel do schedule(dynamic) num_threads(nthreads)
      do i=1,size(dec%object)
        part(i)=mflsss_obj_run(dec%object(i),need,lim)
      end do
!$omp end parallel do
    else
      call timer%start(lim)
      do i=1,size(dec%object)
        remain=max(0.0_dp,lim-timer%elapsed())
        if(remain<=0.0_dp) then
          part(i)%timed_out=.true.
          cycle
        end if
        part(i)=mflsss_obj_run(dec%object(i),need,remain)
      end do
    end if
    do i=1,size(part)
      call merge_results(r,part(i),need)
      if(r%size()>=need)exit
    end do
    r%partitions_run=int(size(part),i8)
    if(usepar .and. openmp_enabled()) then
      r%engine='decomp-omp'
    else
      r%engine='decomp-serial'
    end if
  end function mflsss_decomp_run

  subroutine remap_solutions(r,order)
    type(subset_solutions),intent(inout)::r
    integer,intent(in)::order(:)
    integer::i,j,k,tmp
    do i=1,r%size()
      do j=1,size(r%sol(i)%idx)
        r%sol(i)%idx(j)=order(r%sol(i)%idx(j))
      end do
      do j=2,size(r%sol(i)%idx)
        tmp=r%sol(i)%idx(j);k=j-1
        do while(k>=1)
          if(r%sol(i)%idx(k)<=tmp)exit
          r%sol(i)%idx(k+1)=r%sol(i)%idx(k);k=k-1
        end do
        r%sol(i)%idx(k+1)=tmp
      end do
    end do
  end subroutine remap_solutions

  logical function mpat_recommended(v,len,target,me,dl,du) result(use_mpat)
    integer(i8),intent(in)::v(:,:),target(:),me(:)
    integer,intent(in)::len
    integer,intent(in),optional::dl,du
    integer::d,nlo,nup,j
    integer(i8)::lo,hi,width,dist
    real(dp)::edge_score
    use_mpat=.false.
    if(len<3 .or. size(v,1)-len<4) return
    if(.not.is_comonotonic_i8(v)) return
    d=size(v,2);nlo=d;nup=d
    if(present(dl))nlo=max(0,min(d,dl))
    if(present(du))nup=max(0,min(d,du))
    edge_score=0.0_dp
    do j=1,max(nlo,nup)
      lo=sum(v(1:len,j));hi=sum(v(size(v,1)-len+1:size(v,1),j))
      width=max(1_i8,hi-lo)
      if(j<=nlo .and. j<=nup)then
        dist=min(abs((target(j)-me(j))-lo),abs(hi-(target(j)+me(j))))
      else if(j<=nlo)then
        dist=abs((target(j)-me(j))-lo)
      else
        dist=abs(hi-(target(j)+me(j)))
      end if
      edge_score=edge_score+real(dist,dp)/real(width,dp)
    end do
    if(max(nlo,nup)>0)edge_score=edge_score/real(max(nlo,nup),dp)
    ! mPAT is cheap enough to use on moderate/large comonotonic searches even
    ! away from the boundary.  Tiny central searches remain on packed DFS.
    use_mpat=(len>=6 .or. size(v,1)>=32 .or. edge_score<0.32_dp)
  end function mpat_recommended

  function lower_ascii(s) result(r)
    character(len=*),intent(in)::s
    character(len=len(s))::r
    integer::i,c
    r=s
    do i=1,len(s)
      c=iachar(r(i:i))
      if(c>=iachar('A') .and. c<=iachar('Z'))r(i:i)=achar(c+32)
    end do
  end function lower_ascii

  subroutine merge_results(a,b,need)
    type(subset_solutions),intent(inout)::a
    type(subset_solutions),intent(in)::b
    integer,intent(in)::need
    integer::i
    do i=1,b%size()
      if(a%size()>=need)exit
      call append_solution(a,b%sol(i)%idx)
    end do
    a%nodes=a%nodes+b%nodes; a%pruned=a%pruned+b%pruned
    a%hash_lookups=a%hash_lookups+b%hash_lookups
    a%hash_candidates=a%hash_candidates+b%hash_candidates
    a%bound_states=a%bound_states+b%bound_states
    a%timed_out=a%timed_out.or.b%timed_out
  end subroutine merge_results

end module flsss_api
