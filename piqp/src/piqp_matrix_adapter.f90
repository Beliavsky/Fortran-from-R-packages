! SPDX-License-Identifier: BSD-2-Clause
module piqp_matrix_adapter
   use piqp_kinds, only : dp
   use piqp_types, only : piqp_settings_type, piqp_result_type
   use piqp_solver, only : solve_piqp_dense
   use matrix_sparse, only : csc_matrix, csr_matrix, csr_from_csc, csr_to_dense
   implicit none
   private
   public :: solve_piqp_sparse, csc_to_dense_piqp
contains

   function csc_to_dense_piqp(a) result(x)
      type(csc_matrix), intent(in) :: a
      real(dp), allocatable :: x(:,:)
      type(csr_matrix) :: r
      call csr_from_csc(a,r)
      x=csr_to_dense(r)
   end function csc_to_dense_piqp

   subroutine solve_piqp_sparse(pmat,c,result,amat,b,gmat,h_l,h_u,x_l,x_u,settings)
      type(csc_matrix),intent(in)::pmat
      real(dp),intent(in)::c(:)
      type(piqp_result_type),intent(out)::result
      type(csc_matrix),intent(in),optional::amat,gmat
      real(dp),intent(in),optional::b(:),h_l(:),h_u(:),x_l(:),x_u(:)
      type(piqp_settings_type),intent(in),optional::settings
      real(dp),allocatable::pd(:,:),ad(:,:),gd(:,:)
      pd=csc_to_dense_piqp(pmat)
      if(present(amat)) then
         ad=csc_to_dense_piqp(amat)
      else
         allocate(ad(0,size(c)))
      end if
      if(present(gmat)) then
         gd=csc_to_dense_piqp(gmat)
      else
         allocate(gd(0,size(c)))
      end if
      call solve_piqp_dense(pd,c,result,ad,b,gd,h_l,h_u,x_l,x_u,settings)
   end subroutine solve_piqp_sparse
end module piqp_matrix_adapter
