! Upstream OOR license declaration: LGPL (version unspecified).
module oor_poo
   use oor_kinds, only : dp
   use oor_interfaces, only : scalar_objective
   use oor_random, only : random_uniform
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
   implicit none
   private
   public :: poo, poo_result, poo_leaf, poo_tree

   type :: poo_leaf
      real(dp) :: noisyvalue = 0.0_dp
      real(dp) :: evaluated = 0.0_dp
      logical :: sampled = .false.
      integer :: father = 0
      integer :: depth = 0
      integer :: children(2) = 0
      real(dp) :: support(2) = [0.0_dp, 1.0_dp]
      real(dp), allocatable :: rhos(:)
      integer, allocatable :: visited(:)
      real(dp), allocatable :: rewards(:)
      real(dp), allocatable :: bvalue(:)
      real(dp), allocatable :: uvalue(:)
      real(dp), allocatable :: mbvalue(:)
      real(dp), allocatable :: muvalue(:)
   end type poo_leaf

   type :: poo_tree
      type(poo_leaf), allocatable :: leaves(:)
   end type poo_tree

   type :: poo_result
      real(dp) :: par = 0.0_dp
      real(dp) :: value = 0.0_dp
      real(dp) :: best_rho = 0.0_dp
      integer :: evaluations = 0
      integer :: samples = 0
      type(poo_tree) :: tree
   end type poo_result

contains
   subroutine poo(f, horizon, noise_level, result, rhomax, nu)
      procedure(scalar_objective) :: f
      integer, intent(in) :: horizon
      real(dp), intent(in) :: noise_level
      type(poo_result), intent(out) :: result
      integer, intent(in), optional :: rhomax
      real(dp), intent(in), optional :: nu
      integer :: nrho, k, compt, best_k
      real(dp) :: alpha, nu_use
      real(dp), allocatable :: rhos(:), empp(:)
      integer, allocatable :: nsam(:)
      real(dp) :: x, val
      logical :: existed

      if (horizon < 2) error stop "poo: horizon must be at least 2"
      if (noise_level < 0.0_dp) error stop "poo: noise_level must be nonnegative"
      nrho = 20
      if (present(rhomax)) nrho = rhomax
      if (nrho < 1) error stop "poo: rhomax must be positive"
      nu_use = 1.0_dp
      if (present(nu)) nu_use = nu
      if (nu_use <= 0.0_dp) error stop "poo: nu must be positive"

      allocate(rhos(nrho), empp(nrho), nsam(nrho))
      if (nrho == 1) then
         rhos(1) = 0.0_dp
      else
         do k = 1, nrho
            rhos(k) = real(k - 1, dp) / real(nrho - 1, dp)
         end do
      end if
      empp = 0.0_dp
      nsam = 0
      alpha = log(real(horizon, dp)) * noise_level**2
      call init_tree(result%tree, rhos)

      compt = 0
      best_k = 1
      do while (compt <= horizon)
         do k = 1, nrho
            call tree_sample(result%tree, f, alpha, nu_use, k, x, val, existed)
            empp(k) = empp(k) + val
            if (existed) compt = compt + 1
            nsam(k) = nsam(k) + 1
            if (existed .and. compt <= horizon) best_k = best_empirical_index(empp, nsam)
         end do
      end do

      call tree_sample(result%tree, f, alpha, nu_use, best_k, x, val, existed)
      if (existed) compt = compt + 1
      result%par = x
      result%value = val
      result%best_rho = rhos(best_k)
      result%evaluations = compt
      result%samples = sum(nsam) + 1
   end subroutine poo

   subroutine init_tree(tree, rhos)
      type(poo_tree), intent(out) :: tree
      real(dp), intent(in) :: rhos(:)
      allocate(tree%leaves(1))
      call init_leaf(tree%leaves(1), [0.0_dp, 1.0_dp], 0, 0, rhos)
   end subroutine init_tree

   subroutine init_leaf(leaf, support, father, depth, rhos)
      type(poo_leaf), intent(out) :: leaf
      real(dp), intent(in) :: support(2), rhos(:)
      integer, intent(in) :: father, depth
      real(dp) :: pinf, ninf
      integer :: n

      pinf = ieee_value(0.0_dp, ieee_positive_inf)
      ninf = ieee_value(0.0_dp, ieee_negative_inf)
      n = size(rhos)
      leaf%support = support
      leaf%father = father
      leaf%depth = depth
      leaf%children = 0
      leaf%sampled = .false.
      allocate(leaf%rhos(n), leaf%visited(n), leaf%rewards(n), leaf%bvalue(n), &
               leaf%uvalue(n), leaf%mbvalue(n), leaf%muvalue(n))
      leaf%rhos = rhos
      leaf%visited = 0
      leaf%rewards = 0.0_dp
      leaf%bvalue = pinf
      leaf%uvalue = pinf
      leaf%mbvalue = ninf
      leaf%muvalue = ninf
   end subroutine init_leaf

   subroutine add_children(tree, id)
      type(poo_tree), intent(inout) :: tree
      integer, intent(in) :: id
      type(poo_leaf), allocatable :: tmp(:)
      integer :: nold
      real(dp) :: m

      nold = size(tree%leaves)
      m = 0.5_dp * sum(tree%leaves(id)%support)
      allocate(tmp(nold + 2))
      tmp(1:nold) = tree%leaves
      call init_leaf(tmp(nold + 1), [tree%leaves(id)%support(1), m], id, &
                     tree%leaves(id)%depth + 1, tree%leaves(id)%rhos)
      call init_leaf(tmp(nold + 2), [m, tree%leaves(id)%support(2)], id, &
                     tree%leaves(id)%depth + 1, tree%leaves(id)%rhos)
      call move_alloc(tmp, tree%leaves)
      tree%leaves(id)%children = [nold + 1, nold + 2]
   end subroutine add_children

   recursive function explore(tree, id, k) result(idleaf)
      type(poo_tree), intent(inout) :: tree
      integer, intent(in) :: id, k
      integer :: idleaf, pick, c1, c2
      real(dp) :: u

      if (tree%leaves(id)%visited(k) == 0) then
         idleaf = id
         return
      end if
      if (tree%leaves(id)%children(1) == 0) then
         call add_children(tree, id)
         call random_number(u)
         if (u < 0.5_dp) then
            pick = 1
         else
            pick = 2
         end if
         idleaf = tree%leaves(id)%children(pick)
         return
      end if

      c1 = tree%leaves(id)%children(1)
      c2 = tree%leaves(id)%children(2)
      if (tree%leaves(c1)%bvalue(k) >= tree%leaves(c2)%bvalue(k)) then
         idleaf = explore(tree, c1, k)
      else
         idleaf = explore(tree, c2, k)
      end if
   end function explore

   subroutine update_node(tree, id, alpha, nu, k)
      type(poo_tree), intent(inout) :: tree
      integer, intent(in) :: id, k
      real(dp), intent(in) :: alpha, nu
      real(dp) :: meanv, ucb, metric

      meanv = tree%leaves(id)%rewards(k) / real(tree%leaves(id)%visited(k), dp)
      ucb = sqrt(2.0_dp * alpha / real(tree%leaves(id)%visited(k), dp))
      metric = nu * source_metric_power(tree%leaves(id)%rhos(k), tree%leaves(id)%depth, k)
      tree%leaves(id)%uvalue(k) = meanv + ucb + metric
      tree%leaves(id)%muvalue(k) = meanv - ucb - metric
   end subroutine update_node

   recursive subroutine update_path(tree, id, reward, alpha, nu, k)
      type(poo_tree), intent(inout) :: tree
      integer, intent(in) :: id, k
      real(dp), intent(in) :: reward, alpha, nu
      integer :: c1, c2

      tree%leaves(id)%rewards(k) = tree%leaves(id)%rewards(k) + reward
      tree%leaves(id)%visited(k) = tree%leaves(id)%visited(k) + 1
      call update_node(tree, id, alpha, nu, k)
      c1 = tree%leaves(id)%children(1)
      c2 = tree%leaves(id)%children(2)
      if (c1 == 0) then
         tree%leaves(id)%bvalue(k) = tree%leaves(id)%uvalue(k)
         tree%leaves(id)%mbvalue(k) = tree%leaves(id)%muvalue(k)
      else
         tree%leaves(id)%bvalue(k) = min(tree%leaves(id)%uvalue(k), &
              max(tree%leaves(c1)%bvalue(k), tree%leaves(c2)%bvalue(k)))
         tree%leaves(id)%mbvalue(k) = max(tree%leaves(id)%muvalue(k), &
              max(tree%leaves(c1)%mbvalue(k), tree%leaves(c2)%mbvalue(k)))
      end if
      if (tree%leaves(id)%father /= 0) then
         call update_path(tree, tree%leaves(id)%father, reward, alpha, nu, k)
      end if
   end subroutine update_path

   subroutine tree_sample(tree, f, alpha, nu, k, evaluated, noisyvalue, existed)
      type(poo_tree), intent(inout) :: tree
      procedure(scalar_objective) :: f
      real(dp), intent(in) :: alpha, nu
      integer, intent(in) :: k
      real(dp), intent(out) :: evaluated, noisyvalue
      logical, intent(out) :: existed
      integer :: idleaf

      idleaf = explore(tree, 1, k)
      existed = .false.
      if (.not. tree%leaves(idleaf)%sampled) then
         tree%leaves(idleaf)%evaluated = random_uniform(tree%leaves(idleaf)%support(1), &
                                                        tree%leaves(idleaf)%support(2))
         tree%leaves(idleaf)%noisyvalue = f(tree%leaves(idleaf)%evaluated)
         tree%leaves(idleaf)%sampled = .true.
         existed = .true.
      end if
      call update_path(tree, idleaf, tree%leaves(idleaf)%noisyvalue, alpha, nu, k)
      evaluated = tree%leaves(idleaf)%evaluated
      noisyvalue = tree%leaves(idleaf)%noisyvalue
   end subroutine tree_sample

   function source_metric_power(rho, depth, k) result(v)
      real(dp), intent(in) :: rho
      integer, intent(in) :: depth, k
      real(dp) :: v, exponent

      if (depth == 0) then
         v = 1.0_dp
         return
      end if
      if (rho <= 0.0_dp) then
         v = 0.0_dp
      else if (rho >= 1.0_dp) then
         v = 1.0_dp
      else
         exponent = real(depth, dp)**real(k, dp)
         if (exponent > 1.0e6_dp) then
            v = 0.0_dp
         else
            v = exp(log(rho) * exponent)
         end if
      end if
   end function source_metric_power

   function best_empirical_index(empp, nsam) result(idx)
      real(dp), intent(in) :: empp(:)
      integer, intent(in) :: nsam(:)
      integer :: idx, i
      real(dp) :: best, avg
      real(dp) :: ninf

      ninf = ieee_value(0.0_dp, ieee_negative_inf)
      best = ninf
      idx = 1
      do i = 1, size(empp)
         if (nsam(i) > 0) then
            avg = empp(i) / real(nsam(i), dp)
            if (avg > best) then
               best = avg
               idx = i
            end if
         end if
      end do
   end function best_empirical_index
end module oor_poo
