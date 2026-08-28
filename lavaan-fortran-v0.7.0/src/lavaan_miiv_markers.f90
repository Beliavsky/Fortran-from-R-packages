module lavaan_miiv_markers
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model
   use lavaan_linalg, only : inverse_general
   implicit none
   private

   type, public :: miiv_marker
      integer :: latent_node = 0
      integer :: marker_node = 0
   end type miiv_marker

   type, public :: miiv_marker_equation
      integer :: outcome_node = 0
      integer :: proxy_outcome_node = 0
      integer, allocatable :: predictor_nodes(:), proxy_predictor_nodes(:)
      integer, allocatable :: instrument_nodes(:)
      real(dp), allocatable :: proxy_coefficients(:)
      logical :: identified = .false.
   end type miiv_marker_equation

   public :: ram_miiv_marker_equations, miiv_proxy_node

contains

   subroutine ram_miiv_marker_equations(model, markers, equations, status, tol)
      type(ram_model), intent(in) :: model
      type(miiv_marker), intent(in) :: markers(:)
      type(miiv_marker_equation), allocatable, intent(out) :: equations(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: ia(:, :), inv(:, :), covall(:, :), resid_coef(:)
      integer, allocatable :: outcomes(:), pred(:), ppred(:), cand(:), inst(:)
      real(dp), allocatable :: bproxy(:)
      real(dp) :: eps, ly, lp, ce
      integer :: n, i, j, m, neq, info, out, pout, node, ninst

      n = size(model%a, 1)
      eps = 1.0e-8_dp
      if (present(tol)) eps = tol
      if (size(model%a, 2) /= n .or. any(shape(model%s) /= [n, n])) then
         status = -1
         allocate(equations(0))
         return
      end if
      call validate_markers(model, markers, status, eps)
      if (status /= 0) then
      allocate(equations(0))
      return
      end if
      allocate(ia(n, n))
      ia = -model%a
      do i = 1, n
      ia(i, i) = ia(i, i) + 1.0_dp
      end do
      call inverse_general(ia, inv, info)
      if (info /= 0) then
      status = 100 + info
      allocate(equations(0))
      return
      end if
      covall = matmul(inv, matmul(model%s, transpose(inv)))
      if (allocated(model%observed)) then
         allocate(cand(size(model%observed)))
         cand = model%observed
      else
         allocate(cand(n))
         do i = 1, n
         cand(i) = i
         end do
      end if

      allocate(outcomes(n))
      neq = 0
      do i = 1, n
         pred = pack([(j, j=1,n)], abs(model%a(i, :)) > eps)
         if (size(pred) > 0) then
            pout = miiv_proxy_node(i, model, markers, ly)
            if (pout == 0) cycle
            if (all([(miiv_proxy_node(pred(j), model, markers, lp) > 0, j=1,size(pred))])) then
               neq = neq + 1
               outcomes(neq) = i
            end if
         end if
      end do
      allocate(equations(neq))
      do m = 1, neq
         out = outcomes(m)
         pred = pack([(j, j=1,n)], abs(model%a(out, :)) > eps)
         pout = miiv_proxy_node(out, model, markers, ly)
         allocate(ppred(size(pred)), bproxy(size(pred)))
         do j = 1, size(pred)
            ppred(j) = miiv_proxy_node(pred(j), model, markers, lp)
            bproxy(j) = ly * model%a(out, pred(j)) / lp
         end do
         allocate(resid_coef(n))
         resid_coef = 0.0_dp
         resid_coef(pout) = 1.0_dp
         do j = 1, size(ppred)
         resid_coef(ppred(j)) = resid_coef(ppred(j)) - bproxy(j)
         end do
         allocate(inst(size(cand)))
         ninst = 0
         do j = 1, size(cand)
            node = cand(j)
            if (node == pout .or. any(node == ppred)) cycle
            ! Descendants of the structural outcome are not valid model-implied instruments.
            if (abs(inv(node, out)) > eps) cycle
            ce = dot_product(covall(node, :), resid_coef)
            if (abs(ce) > 20.0_dp * eps) cycle
            if (maxval(abs(covall(node, ppred))) <= eps) cycle
            ninst = ninst + 1
            inst(ninst) = node
         end do
         equations(m)%outcome_node = out
         equations(m)%proxy_outcome_node = pout
         equations(m)%predictor_nodes = pred
         equations(m)%proxy_predictor_nodes = ppred
         equations(m)%proxy_coefficients = bproxy
         allocate(equations(m)%instrument_nodes(ninst))
         if (ninst > 0) equations(m)%instrument_nodes = inst(1:ninst)
         equations(m)%identified = ninst >= size(pred)
         deallocate(ppred, bproxy, resid_coef, inst)
      end do
      status = 0
   end subroutine ram_miiv_marker_equations

   integer function miiv_proxy_node(node, model, markers, loading) result(proxy)
      integer, intent(in) :: node
      type(ram_model), intent(in) :: model
      type(miiv_marker), intent(in) :: markers(:)
      real(dp), intent(out) :: loading
      integer :: j
      proxy = 0
      loading = 0.0_dp
      if (allocated(model%observed)) then
         if (any(model%observed == node)) then
            proxy = node
            loading = 1.0_dp
            return
         end if
      end if
      do j = 1, size(markers)
         if (markers(j)%latent_node == node) then
            proxy = markers(j)%marker_node
            loading = model%a(proxy, node)
            return
         end if
      end do
   end function miiv_proxy_node

   subroutine validate_markers(model, markers, status, tol)
      type(ram_model), intent(in) :: model
      type(miiv_marker), intent(in) :: markers(:)
      integer, intent(out) :: status
      real(dp), intent(in) :: tol
      integer :: n, j, k, l, m
      n = size(model%a, 1)
      status = 0
      do j = 1, size(markers)
         l = markers(j)%latent_node
         m = markers(j)%marker_node
         if (l < 1 .or. l > n .or. m < 1 .or. m > n .or. l == m) then
            status = -10 - j
            return
         end if
         if (.not.allocated(model%observed) .or. .not.any(model%observed == m)) then
            status = -100 - j
            return
         end if
         if (abs(model%a(m, l)) <= tol) then
            status = -200 - j
            return
         end if
         if (count([(markers(k)%latent_node == l, k=1,size(markers))]) > 1) then
            status = -300 - j
            return
         end if
      end do
   end subroutine validate_markers

end module lavaan_miiv_markers
