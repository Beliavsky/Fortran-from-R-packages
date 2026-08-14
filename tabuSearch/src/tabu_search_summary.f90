module tabu_search_summary
   use tabu_search_kinds, only : dp
   use tabu_search_core, only : tabu_result
   implicit none
   private

   type, public :: tabu_summary_result
      real(dp) :: best_value = -huge(1.0_dp)
      integer :: best_count = 0
      integer :: total_iterations = 0
      integer :: unique_configurations = 0
      integer, allocatable :: optimum_iterations(:)
      integer, allocatable :: optimum_nvars(:)
      integer, allocatable :: optimum_configurations(:,:)
      integer, allocatable :: selected_count(:)
      integer, allocatable :: move_frequency(:)
   end type tabu_summary_result

   public :: summarize_tabu

contains

   subroutine summarize_tabu(result, summary)
      type(tabu_result), intent(in) :: result
      type(tabu_summary_result), intent(out) :: summary
      integer :: i, j, nb, count_best
      logical :: seen

      summary%total_iterations = result%n_records
      if (result%n_records == 0) return

      summary%best_value = maxval(result%utility_keep)
      count_best = 0
      do i = 1, result%n_records
         if (result%utility_keep(i) <= summary%best_value .and. &
             result%utility_keep(i) >= summary%best_value) count_best = count_best + 1
      end do
      summary%best_count = count_best
      allocate(summary%optimum_iterations(count_best))
      allocate(summary%optimum_nvars(count_best))
      allocate(summary%optimum_configurations(count_best, result%size))
      allocate(summary%selected_count(result%size), summary%move_frequency(result%size))

      summary%selected_count = sum(result%config_keep, dim=1)
      summary%move_frequency = 0
      do i = 2, result%n_records
         where (result%config_keep(i,:) /= result%config_keep(i-1,:))
            summary%move_frequency = summary%move_frequency + 1
         end where
      end do

      nb = 0
      do i = 1, result%n_records
         if (result%utility_keep(i) <= summary%best_value .and. &
             result%utility_keep(i) >= summary%best_value) then
            nb = nb + 1
            summary%optimum_iterations(nb) = i
            summary%optimum_nvars(nb) = sum(result%config_keep(i,:))
            summary%optimum_configurations(nb,:) = result%config_keep(i,:)
         end if
      end do

      summary%unique_configurations = 0
      do i = 1, result%n_records
         seen = .false.
         do j = 1, i - 1
            if (all(result%config_keep(i,:) == result%config_keep(j,:))) then
               seen = .true.
               exit
            end if
         end do
         if (.not. seen) summary%unique_configurations = summary%unique_configurations + 1
      end do
   end subroutine summarize_tabu

end module tabu_search_summary
