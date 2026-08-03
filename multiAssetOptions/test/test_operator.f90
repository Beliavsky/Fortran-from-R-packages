program test_operator
   use multi_asset_options
   implicit none

   type(grid_set) :: grid
   type(csr_matrix) :: operator
   type(status_type) :: status
   real(dp), allocatable :: dense(:,:), ones(:), product(:)
   real(dp), parameter :: tol = 1.0e-12_dp
   real(dp) :: expected(3,3)

   allocate(grid%asset(1),grid%dims(1),grid%strides(1))
   grid%asset(1)%x = [0.0_dp,1.0_dp,2.0_dp]
   grid%dims = [3]
   grid%strides = [1]
   grid%n_nodes = 3

   call build_fdm_operator(grid,0.1_dp,[0.02_dp],[0.2_dp], &
      reshape([1.0_dp],[1,1]),operator,status)
   call require(status%ok(),status_message(status))
   call csr_to_dense(operator,dense,status)
   call require(status%ok(),status_message(status))

   expected = 0.0_dp
   expected(1,:) = [-0.1_dp,0.0_dp,0.0_dp]
   expected(2,:) = [-0.02_dp,-0.14_dp,0.06_dp]
   expected(3,:) = [0.0_dp,-0.16_dp,0.06_dp]
   call require(maxval(abs(dense-expected)) < tol, &
      'one-dimensional FDM operator does not match the reference stencil')

   allocate(ones(3),product(3))
   ones = 1.0_dp
   call csr_matvec(operator,ones,product)
   call require(maxval(abs(product+0.1_dp)) < tol, &
      'FDM operator does not discount a constant state correctly')


   deallocate(grid%asset,grid%dims,grid%strides,ones,product)
   allocate(grid%asset(2),grid%dims(2),grid%strides(2))
   grid%asset(1)%x = [0.0_dp,1.0_dp,2.0_dp]
   grid%asset(2)%x = [0.0_dp,1.0_dp,2.0_dp]
   grid%dims = [3,3]
   grid%strides = [1,3]
   grid%n_nodes = 9
   call build_fdm_operator(grid,0.0_dp,[0.0_dp,0.0_dp], &
      [0.2_dp,0.3_dp],reshape([1.0_dp,0.5_dp,0.5_dp,1.0_dp],[2,2]), &
      operator,status)
   call require(status%ok(),status_message(status))
   call csr_to_dense(operator,dense,status)
   call require(status%ok(),status_message(status))
   call require(abs(dense(5,1)-0.0075_dp) < tol .and. &
      abs(dense(5,3)+0.0075_dp) < tol .and. &
      abs(dense(5,7)+0.0075_dp) < tol .and. &
      abs(dense(5,9)-0.0075_dp) < tol, &
      'two-dimensional mixed-derivative stencil is incorrect')

   print '(a)', 'test_operator: PASS'

contains

   subroutine require(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*,'(a)') 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine require

   function status_message(status_value) result(message)
      type(status_type), intent(in) :: status_value
      character(len=:), allocatable :: message
      if (allocated(status_value%message)) then
         message = status_value%message
      else
         message = 'unexpected status failure'
      end if
   end function status_message

end program test_operator
