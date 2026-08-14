module tabu_search
   use tabu_search_kinds, only : dp, i8
   use tabu_search_core, only : tabu_control, tabu_result, tabu_objective, run_tabu_search => tabu_search
   use tabu_search_summary, only : tabu_summary_result, summarize_tabu
   implicit none
   public
end module tabu_search
