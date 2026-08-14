module tabu_search_core
   use tabu_search_kinds, only : dp, i8
   use tabu_search_rng, only : tabu_rng
   implicit none
   private

   abstract interface
      function tabu_objective(config) result(value)
         import :: dp
         integer, intent(in) :: config(:)
         real(dp) :: value
      end function tabu_objective
   end interface

   type, public :: tabu_control
      integer :: iters = 100
      integer :: neigh = 0
      integer :: list_size = 9
      integer :: n_restarts = 10
      integer :: repeat_all = 1
      integer(i8) :: seed = 123456789_i8
      logical :: verbose = .false.
   end type tabu_control

   type, public :: tabu_result
      integer :: size = 0
      integer :: n_records = 0
      integer :: iters = 0
      integer :: neigh = 0
      integer :: list_size = 0
      integer :: repeat_all = 0
      integer, allocatable :: config_keep(:,:)
      real(dp), allocatable :: utility_keep(:)
   contains
      procedure :: best_value => tabu_best_value
      procedure :: best_configuration => tabu_best_configuration
   end type tabu_result

   public :: tabu_search, tabu_objective

contains

   subroutine tabu_search(size_config, objective, result, control, initial_config)
      integer, intent(in) :: size_config
      procedure(tabu_objective) :: objective
      type(tabu_result), intent(out) :: result
      type(tabu_control), intent(in), optional :: control
      integer, intent(in), optional :: initial_config(:)

      type(tabu_control) :: ctl
      type(tabu_rng) :: rng
      integer :: neigh, max_records, nrec, repeat_index, restarts
      integer :: best_index
      integer, allocatable :: history(:,:), config(:), tabu_list(:), list_order(:)
      integer, allocatable :: random_neighbours(:), frequency(:), selected_tabu(:)
      integer, allocatable :: candidate(:), ties(:), permutation(:)
      real(dp), allocatable :: utility_history(:), neighbour_utility(:)
      logical, allocatable :: evaluated(:)
      real(dp) :: eutility, aspiration, tempo

      ctl = tabu_control()
      if (present(control)) ctl = control
      neigh = ctl%neigh
      if (neigh == 0) neigh = size_config
      call validate_inputs(size_config, neigh, ctl, initial_config)

      max_records = ctl%repeat_all * ctl%iters * (ctl%n_restarts + 2)
      allocate(history(max_records, size_config), utility_history(max_records))
      allocate(config(size_config), tabu_list(size_config), list_order(ctl%list_size))
      allocate(random_neighbours(neigh), frequency(size_config), selected_tabu(ctl%list_size))
      allocate(candidate(size_config), ties(size_config), permutation(size_config))
      allocate(neighbour_utility(size_config), evaluated(size_config))
      history = 0
      utility_history = 0.0_dp
      nrec = 0
      call rng%seed(ctl%seed)

      if (present(initial_config)) then
         config = initial_config
      else
         call random_configuration(rng, config)
      end if

      do repeat_index = 1, ctl%repeat_all
         if (repeat_index > 1) call random_configuration(rng, config)

         tabu_list = 0
         list_order = 0
         eutility = objective(config)
         aspiration = eutility

         if (ctl%verbose) write(*,'(a)') "Preliminary search stage..."
         call preliminary_search()

         tempo = -huge(1.0_dp)
         restarts = 0
         do while (tempo < aspiration .and. restarts < ctl%n_restarts)
            if (ctl%verbose) write(*,'(a)') "Intensification stage..."
            best_index = last_best_index(utility_history, nrec)
            eutility = utility_history(best_index)
            tempo = aspiration
            config = history(best_index,:)
            call preliminary_search()
            restarts = restarts + 1
         end do

         if (ctl%verbose) write(*,'(a)') "Diversification stage..."
         call random_configuration(rng, config)
         eutility = objective(config)
         call move_frequency(history, nrec, frequency)
         call choose_frequent_tabu(rng, frequency, ctl%list_size, selected_tabu)
         tabu_list = 0
         tabu_list(selected_tabu) = 1
         list_order = selected_tabu
         call rng%shuffle(list_order)
         call preliminary_search()
      end do

      result%size = size_config
      result%n_records = nrec
      result%iters = ctl%iters
      result%neigh = neigh
      result%list_size = ctl%list_size
      result%repeat_all = ctl%repeat_all
      allocate(result%config_keep(nrec, size_config), result%utility_keep(nrec))
      result%config_keep = history(1:nrec,:)
      result%utility_keep = utility_history(1:nrec)

   contains

      subroutine preliminary_search()
         real(dp) :: max_non_tabu, max_tabu, new_utility
         integer :: step, j, move, ntie, fallback, old_move, insert_pos
         logical :: have_non_tabu, have_tabu

         call record_state()
         do step = 2, ctl%iters
            neighbour_utility = -huge(1.0_dp)
            evaluated = .false.
            permutation = [(j, j=1,size_config)]
            call rng%shuffle(permutation)
            random_neighbours = permutation(1:neigh)

            do j = 1, neigh
               move = random_neighbours(j)
               candidate = config
               candidate(move) = 1 - candidate(move)
               neighbour_utility(move) = objective(candidate)
               evaluated(move) = .true.
            end do

            call best_by_tabu(.false., max_non_tabu, have_non_tabu)
            call best_by_tabu(.true., max_tabu, have_tabu)

            if (.not. have_non_tabu) then
               fallback = choose_random_non_tabu()
               candidate = config
               candidate(fallback) = 1 - candidate(fallback)
               neighbour_utility(fallback) = objective(candidate)
               evaluated(fallback) = .true.
               call best_by_tabu(.false., max_non_tabu, have_non_tabu)
            end if

            if (have_tabu .and. max_tabu > max_non_tabu .and. max_tabu > aspiration) then
               call collect_ties(.true., max_tabu, ties, ntie)
            else
               call collect_ties(.false., max_non_tabu, ties, ntie)
            end if
            move = ties(rng%randint(1, ntie))
            new_utility = neighbour_utility(move)

            if (eutility >= new_utility) then
               if (tabu_list(move) == 0) then
                  tabu_list(move) = 1
                  if (count(tabu_list == 1) > ctl%list_size) then
                     old_move = list_order(1)
                     if (old_move > 0) tabu_list(old_move) = 0
                     if (ctl%list_size > 1) then
                        list_order(1:ctl%list_size-1) = list_order(2:ctl%list_size)
                     end if
                     list_order(ctl%list_size) = 0
                  end if
                  insert_pos = first_zero(list_order)
                  if (insert_pos > 0) list_order(insert_pos) = move
               end if
            else if (new_utility > aspiration) then
               aspiration = new_utility
            end if

            eutility = new_utility
            config(move) = 1 - config(move)
            call record_state()
         end do
      end subroutine preliminary_search

      subroutine best_by_tabu(want_tabu, best, have_value)
         logical, intent(in) :: want_tabu
         real(dp), intent(out) :: best
         logical, intent(out) :: have_value
         integer :: idx

         best = -huge(1.0_dp)
         have_value = .false.
         do idx = 1, size_config
            if (.not. evaluated(idx)) cycle
            if ((tabu_list(idx) == 1) .neqv. want_tabu) cycle
            if (.not. have_value .or. neighbour_utility(idx) > best) then
               best = neighbour_utility(idx)
               have_value = .true.
            end if
         end do
      end subroutine best_by_tabu

      subroutine collect_ties(want_tabu, best, tie_index, n_ties)
         logical, intent(in) :: want_tabu
         real(dp), intent(in) :: best
         integer, intent(out) :: tie_index(:)
         integer, intent(out) :: n_ties
         integer :: idx

         n_ties = 0
         do idx = 1, size_config
            if (.not. evaluated(idx)) cycle
            if ((tabu_list(idx) == 1) .neqv. want_tabu) cycle
            if (neighbour_utility(idx) <= best .and. neighbour_utility(idx) >= best) then
               n_ties = n_ties + 1
               tie_index(n_ties) = idx
            end if
         end do
         if (n_ties == 0) error stop "tabu_search: no admissible neighbour"
      end subroutine collect_ties

      function choose_random_non_tabu() result(index_value)
         integer :: index_value
         integer :: idx, nfree, pick

         nfree = count(tabu_list == 0)
         if (nfree == 0) error stop "tabu_search: tabu list covers all moves"
         pick = rng%randint(1, nfree)
         nfree = 0
         index_value = 0
         do idx = 1, size_config
            if (tabu_list(idx) /= 0) cycle
            nfree = nfree + 1
            if (nfree == pick) then
               index_value = idx
               exit
            end if
         end do
      end function choose_random_non_tabu

      subroutine record_state()
         nrec = nrec + 1
         if (nrec > max_records) error stop "tabu_search: history capacity exceeded"
         history(nrec,:) = config
         utility_history(nrec) = eutility
      end subroutine record_state

   end subroutine tabu_search

   subroutine validate_inputs(size_config, neigh, ctl, initial_config)
      integer, intent(in) :: size_config, neigh
      type(tabu_control), intent(in) :: ctl
      integer, intent(in), optional :: initial_config(:)

      if (size_config < 2) error stop "tabu_search: size must be at least 2"
      if (ctl%iters < 2) error stop "tabu_search: iters must be at least 2"
      if (ctl%list_size < 1 .or. ctl%list_size >= size_config) then
         error stop "tabu_search: list_size must be in [1,size-1]"
      end if
      if (neigh < 1 .or. neigh > size_config) then
         error stop "tabu_search: neigh must be in [1,size]"
      end if
      if (ctl%n_restarts < 0) error stop "tabu_search: n_restarts must be nonnegative"
      if (ctl%repeat_all < 1) error stop "tabu_search: repeat_all must be positive"
      if (present(initial_config)) then
         if (size(initial_config) /= size_config) then
            error stop "tabu_search: initial_config has wrong size"
         end if
         if (any(initial_config /= 0 .and. initial_config /= 1)) then
            error stop "tabu_search: initial_config must contain only 0 and 1"
         end if
      end if
   end subroutine validate_inputs

   subroutine random_configuration(rng, config)
      type(tabu_rng), intent(inout) :: rng
      integer, intent(out) :: config(:)
      integer, allocatable :: index(:)
      integer :: i, k

      allocate(index(size(config)))
      index = [(i, i=1,size(config))]
      call rng%shuffle(index)
      k = rng%randint(1, size(config))
      config = 0
      config(index(1:k)) = 1
   end subroutine random_configuration

   subroutine move_frequency(history, nrec, frequency)
      integer, intent(in) :: history(:,:)
      integer, intent(in) :: nrec
      integer, intent(out) :: frequency(:)
      integer :: i

      frequency = 0
      if (nrec < 2) return
      do i = 2, nrec
         where (history(i,:) /= history(i-1,:))
            frequency = frequency + 1
         end where
      end do
   end subroutine move_frequency

   subroutine choose_frequent_tabu(rng, frequency, list_size, selected)
      type(tabu_rng), intent(inout) :: rng
      integer, intent(in) :: frequency(:), list_size
      integer, intent(out) :: selected(:)
      integer, allocatable :: order(:)
      integer :: i, j, tmp

      allocate(order(size(frequency)))
      order = [(i, i=1,size(frequency))]
      call rng%shuffle(order)
      do i = 2, size(order)
         tmp = order(i)
         j = i - 1
         do while (j >= 1)
            if (frequency(order(j)) >= frequency(tmp)) exit
            order(j+1) = order(j)
            j = j - 1
         end do
         order(j+1) = tmp
      end do
      selected = order(1:list_size)
   end subroutine choose_frequent_tabu

   function first_zero(x) result(index_value)
      integer, intent(in) :: x(:)
      integer :: index_value, i

      index_value = 0
      do i = 1, size(x)
         if (x(i) == 0) then
            index_value = i
            return
         end if
      end do
   end function first_zero

   function last_best_index(values, n) result(index_value)
      real(dp), intent(in) :: values(:)
      integer, intent(in) :: n
      integer :: index_value, i
      real(dp) :: best

      if (n < 1) error stop "last_best_index: empty history"
      best = maxval(values(1:n))
      index_value = 1
      do i = 1, n
         if (values(i) >= best) index_value = i
      end do
   end function last_best_index

   function tabu_best_value(self) result(value)
      class(tabu_result), intent(in) :: self
      real(dp) :: value

      if (.not. allocated(self%utility_keep) .or. self%n_records == 0) then
         value = -huge(1.0_dp)
      else
         value = maxval(self%utility_keep)
      end if
   end function tabu_best_value

   function tabu_best_configuration(self) result(config)
      class(tabu_result), intent(in) :: self
      integer, allocatable :: config(:)
      integer :: i, idx
      real(dp) :: best

      if (.not. allocated(self%config_keep) .or. self%n_records == 0) then
         allocate(config(0))
         return
      end if
      best = maxval(self%utility_keep)
      idx = 1
      do i = 1, self%n_records
         if (self%utility_keep(i) >= best) idx = i
      end do
      allocate(config(self%size))
      config = self%config_keep(idx,:)
   end function tabu_best_configuration

end module tabu_search_core
