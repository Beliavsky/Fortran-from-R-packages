module flsss_compat
  use flsss_kinds, only : dp, i8
  use flsss_types
  use flsss_parallel, only : openmp_enabled
  use flsss_api, only : flsss, flsss_multiset, mflsss_par, mflsss_par_impose_bounds, &
    mflsss_par_integerized, mflsss_par_integerized_parallel, mflsss_par_impose_bounds_integerized, &
    decompose_mflsss, mflsss_obj_run, mflsss_decomp_run, &
    mm_knapsack, mm_knapsack_integerized, aux_knapsack01bb, aux_knapsack01dp, gap_solve, &
    aux_gap_bb, aux_gap_bbdp, aux_gap_ga, arb_flsss, decompose_arb_flsss, arb_flsss_obj_run, &
    arb_flsss_decomp_run, &
    build_ksum_hash, add_num_strings
  implicit none
  private
  public :: flsssmultiset, mflssspar, mflsssparimposebounds, mflsssparintegerized
  public :: mflsssparimposeboundsintegerized, decomposemflsss, mflsssobjrun
  public :: mmknapsack, mmknapsackintegerized, gap, auxgapbb, auxgapbbdp, auxgapga
  public :: auxknapsack01bb, auxknapsack01dp, arbflsss, decomposearbflsss, arbflsssobjrun
  public :: ksumhash, addnumstrings

contains

  function flsssmultiset(len,buckets,target,me,solution_need,tlimit) result(r)
    integer,intent(in)::len(:)
    type(real_bucket),intent(in)::buckets(:)
    real(dp),intent(in)::target,me
    integer,intent(in),optional::solution_need
    real(dp),intent(in),optional::tlimit
    type(multiset_solutions)::r
    r=flsss_multiset(len,buckets,target,me,solution_need,tlimit)
  end function flsssmultiset

  function mflssspar(max_core,len,mv,mtarget,mme,solution_need,tlimit,dl,du,use_bisrch_in_fb,avg_thread_load) result(r)
    integer,intent(in)::max_core,len
    real(dp),intent(in)::mv(:,:),mtarget(:),mme(:)
    integer,intent(in),optional::solution_need,dl,du,avg_thread_load
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::use_bisrch_in_fb
    type(subset_solutions)::r
    type(mflsss_decomposition)::dec
    integer::load
    call consume_parallel_controls(max_core,use_bisrch_in_fb,avg_thread_load)
    load=8;if(present(avg_thread_load))load=avg_thread_load
    if(max_core>1 .and. len>0 .and. openmp_enabled()) then
      dec=decompose_mflsss(len,mv,mtarget,mme,dl,du,max_core*load)
      r=mflsss_decomp_run(dec,solution_need,tlimit,.true.,max_core)
    else
      r=mflsss_par(len,mv,mtarget,mme,solution_need,tlimit,dl,du)
    end if
  end function mflssspar

  function mflsssparimposebounds(max_core,len,mv,mtarget,mme,lb,ub,solution_need,tlimit,dl,du, &
                                 use_bisrch_in_fb,avg_thread_load) result(r)
    integer,intent(in)::max_core,len,lb(:),ub(:)
    real(dp),intent(in)::mv(:,:),mtarget(:),mme(:)
    integer,intent(in),optional::solution_need,dl,du,avg_thread_load
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::use_bisrch_in_fb
    type(subset_solutions)::r
    call consume_parallel_controls(max_core,use_bisrch_in_fb,avg_thread_load)
    r=mflsss_par_impose_bounds(len,mv,mtarget,mme,lb,ub,solution_need,tlimit,dl,du)
  end function mflsssparimposebounds

  function mflsssparintegerized(max_core,len,mv,mtarget,mme,precision_level,solution_need,tlimit,dl,du, &
                                return_before_mining,use_bisrch_in_fb,avg_thread_load) result(r)
    integer,intent(in)::max_core,len
    real(dp),intent(in)::mv(:,:),mtarget(:),mme(:)
    integer,intent(in),optional::precision_level(:),solution_need,dl,du,avg_thread_load
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::return_before_mining,use_bisrch_in_fb
    type(integerized_search_result)::r
    call consume_parallel_controls(max_core,use_bisrch_in_fb,avg_thread_load)
    if(max_core>1 .and. len>0 .and. openmp_enabled()) then
      r=mflsss_par_integerized_parallel(len,mv,mtarget,mme,max_core,precision_level,solution_need, &
        tlimit,dl,du,return_before_mining)
    else
      r=mflsss_par_integerized(len,mv,mtarget,mme,precision_level,solution_need,tlimit,dl,du, &
        return_before_mining)
    end if
  end function mflsssparintegerized

  function mflsssparimposeboundsintegerized(max_core,len,mv,mtarget,mme,lb,ub,precision_level, &
      solution_need,tlimit,dl,du,return_before_mining,use_bisrch_in_fb,avg_thread_load) result(r)
    integer,intent(in)::max_core,len,lb(:),ub(:)
    real(dp),intent(in)::mv(:,:),mtarget(:),mme(:)
    integer,intent(in),optional::precision_level(:),solution_need,dl,du,avg_thread_load
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::return_before_mining,use_bisrch_in_fb
    type(integerized_search_result)::r
    call consume_parallel_controls(max_core,use_bisrch_in_fb,avg_thread_load)
    r=mflsss_par_impose_bounds_integerized(len,mv,mtarget,mme,lb,ub,precision_level,solution_need, &
      tlimit,dl,du,return_before_mining)
  end function mflsssparimposeboundsintegerized

  function decomposemflsss(len,mv,mtarget,mme,dl,du,approx_ninstance) result(r)
    integer,intent(in)::len
    real(dp),intent(in)::mv(:,:),mtarget(:),mme(:)
    integer,intent(in),optional::dl,du,approx_ninstance
    type(mflsss_decomposition)::r
    r=decompose_mflsss(len,mv,mtarget,mme,dl,du,approx_ninstance)
  end function decomposemflsss

  function mflsssobjrun(obj,solution_need,tlimit) result(r)
    type(mflsss_object),intent(in)::obj
    integer,intent(in),optional::solution_need
    real(dp),intent(in),optional::tlimit
    type(subset_solutions)::r
    r=mflsss_obj_run(obj,solution_need,tlimit)
  end function mflsssobjrun

  function mmknapsack(max_core,len,items_profits,items_costs,capacities,heuristic,tlimit, &
                      use_bisrch_in_fb,thread_load,verbose) result(r)
    integer,intent(in)::max_core,len
    real(dp),intent(in)::items_profits(:),items_costs(:,:),capacities(:)
    logical,intent(in),optional::heuristic,use_bisrch_in_fb,verbose
    real(dp),intent(in),optional::tlimit
    integer,intent(in),optional::thread_load
    type(knapsack_result)::r
    call consume_parallel_controls(max_core,use_bisrch_in_fb,thread_load)
    if(present(verbose)) then; if(verbose .and. max_core<0) error stop "mmknapsack: invalid max_core"; end if
    r=mm_knapsack(len,items_profits,items_costs,capacities,heuristic,tlimit)
  end function mmknapsack

  function mmknapsackintegerized(max_core,len,items_profits,items_costs,capacities,precision_level, &
      return_before_mining,heuristic,tlimit,use_bisrch_in_fb,thread_load,verbose) result(r)
    integer,intent(in)::max_core,len
    real(dp),intent(in)::items_profits(:),items_costs(:,:),capacities(:)
    integer,intent(in),optional::precision_level(:),thread_load
    logical,intent(in),optional::return_before_mining,heuristic,use_bisrch_in_fb,verbose
    real(dp),intent(in),optional::tlimit
    type(knapsack_result)::r
    call consume_parallel_controls(max_core,use_bisrch_in_fb,thread_load)
    if(present(heuristic)) then; if(heuristic .and. max_core<0) error stop "mmknapsackintegerized"; end if
    if(present(verbose)) then; if(verbose .and. max_core<0) error stop "mmknapsackintegerized"; end if
    r=mm_knapsack_integerized(len,items_profits,items_costs,capacities,precision_level,return_before_mining,tlimit)
  end function mmknapsackintegerized

  function gap(max_core,agents_costs,agents_profits,agents_budgets,heuristic,tlimit,thread_load,verbose) result(r)
    integer,intent(in)::max_core
    real(dp),intent(in)::agents_costs(:,:),agents_profits(:,:),agents_budgets(:)
    logical,intent(in),optional::heuristic,verbose
    real(dp),intent(in),optional::tlimit
    integer,intent(in),optional::thread_load
    type(gap_result)::r
    call consume_parallel_controls(max_core,load=thread_load)
    if(present(verbose)) then; if(verbose .and. max_core<0) error stop "gap: invalid max_core"; end if
    r=gap_solve(agents_costs,agents_profits,agents_budgets,tlimit,heuristic)
  end function gap

  function auxgapbb(cost,profit_or_loss,budget,max_core,tlimit,ub,greedy_branching,optim, &
                    multhread_on,thread_load) result(r)
    real(dp),intent(in)::cost(:,:),profit_or_loss(:,:),budget(:)
    integer,intent(in),optional::max_core,thread_load
    real(dp),intent(in),optional::tlimit
    character(len=*),intent(in),optional::ub,optim,multhread_on
    logical,intent(in),optional::greedy_branching
    type(gap_result)::r
    call consume_text_controls(max_core,thread_load,ub,multhread_on)
    r=aux_gap_bb(cost,profit_or_loss,budget,optim,tlimit,greedy_branching)
  end function auxgapbb

  function auxgapbbdp(cost,profit_or_loss,budget,max_core,tlimit,greedy_branching,optim, &
                      multhread_on,thread_load) result(r)
    real(dp),intent(in)::cost(:,:),profit_or_loss(:,:),budget(:)
    integer,intent(in),optional::max_core,thread_load
    real(dp),intent(in),optional::tlimit
    character(len=*),intent(in),optional::optim,multhread_on
    logical,intent(in),optional::greedy_branching
    type(gap_result)::r
    call consume_text_controls(max_core,thread_load,multhread_on=multhread_on)
    r=aux_gap_bbdp(cost,profit_or_loss,budget,optim,tlimit,greedy_branching)
  end function auxgapbbdp

  function auxgapga(cost,profit_or_loss,budget,trials,population_size,generations,random_seed,max_core,optim) result(r)
    real(dp),intent(in)::cost(:,:),profit_or_loss(:,:),budget(:)
    integer,intent(in)::trials,population_size,generations
    integer(i8),intent(in),optional::random_seed
    integer,intent(in),optional::max_core
    character(len=*),intent(in),optional::optim
    type(gap_result)::r
    if(present(max_core)) then; if(max_core<1) error stop "auxgapga: max_core must be positive"; end if
    r=aux_gap_ga(cost,profit_or_loss,budget,trials,population_size,generations,random_seed,optim)
  end function auxgapga

  function auxknapsack01bb(weight,value,caps,item_ncaps,max_core,tlimit,ub,simplify) result(r)
    real(dp),intent(in)::weight(:),value(:),caps(:)
    integer,intent(in),optional::item_ncaps(:),max_core
    real(dp),intent(in),optional::tlimit
    character(len=*),intent(in),optional::ub
    logical,intent(in),optional::simplify
    type(knapsack_multi_result)::r
    call consume_text_controls(max_core,ub=ub)
    if(present(simplify)) then; if(simplify .and. size(caps)<0) error stop "auxknapsack01bb"; end if
    r=aux_knapsack01bb(weight,value,caps,item_ncaps,tlimit)
  end function auxknapsack01bb

  function auxknapsack01dp(weight,value,caps,max_core,tlimit,simplify) result(r)
    integer,intent(in)::weight(:),caps(:)
    real(dp),intent(in)::value(:)
    integer,intent(in),optional::max_core
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::simplify
    type(knapsack_multi_result)::r
    if(present(max_core)) then; if(max_core<1) error stop "auxknapsack01dp: max_core"; end if
    if(present(simplify)) then; if(simplify .and. size(caps)<0) error stop "auxknapsack01dp"; end if
    r=aux_knapsack01dp(weight,value,caps,tlimit=tlimit)
  end function auxknapsack01dp

  function arbflsss(len,v,target,given_ksum_table,solution_need,max_core,tlimit,approx_ninstance, &
                    ksum_k,ksum_table_size_scaler,verbose) result(r)
    integer,intent(in)::len
    character(len=*),intent(in)::v(:,:),target(:)
    type(ksum_table),intent(in),optional::given_ksum_table
    integer,intent(in),optional::solution_need,max_core,approx_ninstance,ksum_k,ksum_table_size_scaler
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::verbose
    type(subset_solutions)::r
    type(arb_flsss_decomposition)::dec
    integer::cores,ninst
    call consume_arb_controls(max_core,approx_ninstance,ksum_k,ksum_table_size_scaler,verbose)
    cores=1;if(present(max_core))cores=max_core
    ninst=max(1,cores*8);if(present(approx_ninstance))ninst=approx_ninstance
    if(cores>1 .and. openmp_enabled() .and. .not.present(given_ksum_table)) then
      dec=decompose_arb_flsss(len,v,target,ninst)
      r=arb_flsss_decomp_run(dec,solution_need,tlimit,.true.,cores)
    else
      r=arb_flsss(len,v,target,solution_need,tlimit,given_ksum_table)
    end if
  end function arbflsss

  function decomposearbflsss(len,v,target,approx_ninstance,max_core) result(r)
    integer,intent(in)::len
    character(len=*),intent(in)::v(:,:),target(:)
    integer,intent(in),optional::approx_ninstance,max_core
    type(arb_flsss_decomposition)::r
    if(present(max_core)) then; if(max_core<1) error stop "decomposearbflsss: max_core"; end if
    r=decompose_arb_flsss(len,v,target,approx_ninstance)
  end function decomposearbflsss

  function arbflsssobjrun(x,solution_need,tlimit,max_core) result(r)
    type(arb_flsss_object),intent(in)::x
    integer,intent(in),optional::solution_need,max_core
    real(dp),intent(in),optional::tlimit
    type(subset_solutions)::r
    if(present(max_core)) then; if(max_core<1) error stop "arbflsssobjrun: max_core"; end if
    r=arb_flsss_obj_run(x,solution_need,tlimit)
  end function arbflsssobjrun

  function ksumhash(ksum_k,v,ksum_table_size_scaler) result(r)
    integer,intent(in)::ksum_k
    character(len=*),intent(in)::v(:,:)
    integer,intent(in),optional::ksum_table_size_scaler
    type(ksum_table)::r
    integer::cap
    cap=1000000
    if(present(ksum_table_size_scaler)) cap=max(1000,ksum_table_size_scaler*100000)
    r=build_ksum_hash(ksum_k,v,cap)
  end function ksumhash

  function addnumstrings(s) result(r)
    character(len=*),intent(in)::s(:)
    character(len=:),allocatable::r
    r=add_num_strings(s)
  end function addnumstrings

  subroutine consume_parallel_controls(max_core,use_bisrch,load)
    integer,intent(in)::max_core
    logical,intent(in),optional::use_bisrch
    integer,intent(in),optional::load
    if(max_core<1) error stop "FLSSS: max_core must be positive"
    if(present(use_bisrch)) then; if(use_bisrch .and. max_core<0) error stop "FLSSS"; end if
    if(present(load)) then; if(load<1) error stop "FLSSS: thread load must be positive"; end if
  end subroutine consume_parallel_controls

  subroutine consume_text_controls(max_core,load,ub,multhread_on)
    integer,intent(in),optional::max_core,load
    character(len=*),intent(in),optional::ub,multhread_on
    if(present(max_core)) then; if(max_core<1) error stop "FLSSS: max_core"; end if
    if(present(load)) then; if(load<1) error stop "FLSSS: thread load"; end if
    if(present(ub)) then; if(len_trim(ub)==0) error stop "FLSSS: empty ub"; end if
    if(present(multhread_on)) then; if(len_trim(multhread_on)==0) error stop "FLSSS: empty threading mode"; end if
  end subroutine consume_text_controls

  subroutine consume_arb_controls(max_core,approx_ninstance,ksum_k,scaler,verbose)
    integer,intent(in),optional::max_core,approx_ninstance,ksum_k,scaler
    logical,intent(in),optional::verbose
    if(present(max_core)) then; if(max_core<1) error stop "arbFLSSS: max_core"; end if
    if(present(approx_ninstance)) then; if(approx_ninstance<1) error stop "arbFLSSS: approx_ninstance"; end if
    if(present(ksum_k)) then; if(ksum_k<0) error stop "arbFLSSS: ksum_k"; end if
    if(present(scaler)) then; if(scaler<1) error stop "arbFLSSS: ksum table scaler"; end if
    if(present(verbose)) then; if(verbose .and. .false.) error stop "arbFLSSS"; end if
  end subroutine consume_arb_controls

end module flsss_compat
