module lavaan_miiv_partable
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model
   use lavaan_miiv_markers, only : miiv_marker, miiv_marker_equation, ram_miiv_marker_equations
   implicit none
   private

   integer, parameter, public :: miiv_op_regression=1, miiv_op_covariance=2, miiv_op_mean=3, miiv_op_loading=4

   type, public :: miiv_partable_entry
      integer :: lhs_node=0, rhs_node=0, op=0, free_id=0
      real(dp) :: value=0.0_dp
      character(len=64) :: lhs_name='', rhs_name='', label=''
      logical :: marker=.false.
   end type miiv_partable_entry

   type, public :: miiv_named_equation
      integer :: outcome_node=0, proxy_outcome_node=0
      character(len=64) :: outcome_name='', proxy_outcome_name=''
      integer, allocatable :: predictor_nodes(:), proxy_predictor_nodes(:), instrument_nodes(:)
      character(len=64), allocatable :: predictor_names(:), proxy_predictor_names(:), instrument_names(:)
      real(dp), allocatable :: proxy_coefficients(:)
      logical :: identified=.false.
   end type miiv_named_equation

   public :: miiv_auto_markers, ram_miiv_named_equations
   public :: miiv_partable_markers

contains

   subroutine miiv_auto_markers(model, latent_nodes, markers, status, tol)
      type(ram_model), intent(in) :: model
      integer, intent(in) :: latent_nodes(:)
      type(miiv_marker), allocatable, intent(out) :: markers(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: tol
      real(dp) :: eps, val, best
      integer :: i,j,node,obs,bestobs
      eps=1.0e-8_dp
      if(present(tol)) eps=tol
      allocate(markers(size(latent_nodes)))
      if(.not.allocated(model%observed)) then
      status=-1
      return
      end if
      do i=1,size(latent_nodes)
         node=latent_nodes(i)
         bestobs=0
         best=-1.0_dp
         if(node<1 .or. node>size(model%a,1)) then
         status=100+i
         return
         end if
         ! Prefer a unit-loading indicator, matching lavaan's marker-variable convention.
         do j=1,size(model%observed)
            obs=model%observed(j)
            val=model%a(obs,node)
            if(abs(val-1.0_dp)<=100.0_dp*eps) then
            bestobs=obs
            exit
            end if
            if(abs(val)>best .and. abs(val)>eps) then
            best=abs(val)
            bestobs=obs
            end if
         end do
         if(bestobs==0) then
         status=200+i
         return
         end if
         markers(i)%latent_node=node
         markers(i)%marker_node=bestobs
      end do
      status=0
   end subroutine miiv_auto_markers

   subroutine miiv_partable_markers(entries, latent_nodes, markers, status)
      type(miiv_partable_entry), intent(in) :: entries(:)
      integer, intent(in) :: latent_nodes(:)
      type(miiv_marker), allocatable, intent(out) :: markers(:)
      integer, intent(out) :: status
      integer :: i,j,bestj
      real(dp) :: best
      allocate(markers(size(latent_nodes)))
      do i=1,size(latent_nodes)
         bestj=0
         best=-1.0_dp
         do j=1,size(entries)
            if(entries(j)%op/=miiv_op_loading .or. entries(j)%rhs_node/=latent_nodes(i)) cycle
            if(entries(j)%marker) then
            bestj=j
            exit
            end if
            if(abs(entries(j)%value-1.0_dp)<1.0e-10_dp) then
            bestj=j
            exit
            end if
            if(abs(entries(j)%value)>best) then
            best=abs(entries(j)%value)
            bestj=j
            end if
         end do
         if(bestj==0) then
         status=100+i
         return
         end if
         markers(i)%latent_node=latent_nodes(i)
         markers(i)%marker_node=entries(bestj)%lhs_node
      end do
      status=0
   end subroutine miiv_partable_markers

   subroutine ram_miiv_named_equations(model, latent_nodes, node_names, equations, status, markers_in)
      type(ram_model), intent(in) :: model
      integer, intent(in) :: latent_nodes(:)
      character(len=*), intent(in) :: node_names(:)
      type(miiv_named_equation), allocatable, intent(out) :: equations(:)
      integer, intent(out) :: status
      type(miiv_marker), intent(in), optional :: markers_in(:)
      type(miiv_marker), allocatable :: markers(:)
      type(miiv_marker_equation), allocatable :: raw(:)
      integer :: i,j,n
      n=size(model%a,1)
      if(size(node_names)<n) then
      status=-1
      allocate(equations(0))
      return
      end if
      if(present(markers_in)) then
         markers=markers_in
      else
         call miiv_auto_markers(model,latent_nodes,markers,status)
         if(status/=0) then
         allocate(equations(0))
         return
         end if
      end if
      call ram_miiv_marker_equations(model,markers,raw,status)
      if(status/=0) then
      allocate(equations(0))
      return
      end if
      allocate(equations(size(raw)))
      do i=1,size(raw)
         equations(i)%outcome_node=raw(i)%outcome_node
         equations(i)%proxy_outcome_node=raw(i)%proxy_outcome_node
         equations(i)%outcome_name=node_names(raw(i)%outcome_node)
         equations(i)%proxy_outcome_name=node_names(raw(i)%proxy_outcome_node)
         equations(i)%predictor_nodes=raw(i)%predictor_nodes
         equations(i)%proxy_predictor_nodes=raw(i)%proxy_predictor_nodes
         equations(i)%instrument_nodes=raw(i)%instrument_nodes
         equations(i)%proxy_coefficients=raw(i)%proxy_coefficients
         equations(i)%identified=raw(i)%identified
         allocate(equations(i)%predictor_names(size(raw(i)%predictor_nodes)))
         allocate(equations(i)%proxy_predictor_names(size(raw(i)%proxy_predictor_nodes)))
         allocate(equations(i)%instrument_names(size(raw(i)%instrument_nodes)))
         do j=1,size(raw(i)%predictor_nodes)
            equations(i)%predictor_names(j)=node_names(raw(i)%predictor_nodes(j))
         end do
         do j=1,size(raw(i)%proxy_predictor_nodes)
            equations(i)%proxy_predictor_names(j)=node_names(raw(i)%proxy_predictor_nodes(j))
         end do
         do j=1,size(raw(i)%instrument_nodes)
            equations(i)%instrument_names(j)=node_names(raw(i)%instrument_nodes(j))
         end do
      end do
      status=0
   end subroutine ram_miiv_named_equations

end module lavaan_miiv_partable
