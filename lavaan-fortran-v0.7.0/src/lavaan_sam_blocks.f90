module lavaan_sam_blocks
   use lavaan_kinds, only : dp
   implicit none
   private

   type, public :: sam_block_covariance_result
      real(dp), allocatable :: measurement_vcov(:,:)
      real(dp), allocatable :: structural_conditional_vcov(:,:)
      real(dp), allocatable :: cross_vcov(:,:)
      real(dp), allocatable :: structural_marginal_vcov(:,:)
      real(dp), allocatable :: joint_vcov(:,:)
      integer :: status=0
   end type sam_block_covariance_result

   public :: sam_block_covariance
   public :: sam_second_order_bias

contains

   subroutine sam_block_covariance(v_measurement,v_structural_conditional,jac,result)
      real(dp),intent(in)::v_measurement(:,:),v_structural_conditional(:,:),jac(:,:)
      type(sam_block_covariance_result),intent(out)::result
      integer::k1,k2
      k1=size(v_measurement,1)
      k2=size(v_structural_conditional,1)
      if(size(v_measurement,2)/=k1 .or. size(v_structural_conditional,2)/=k2 .or. &
         any(shape(jac)/=[k2,k1])) then
         result%status=-1
         return
      end if
      result%measurement_vcov=0.5_dp*(v_measurement+transpose(v_measurement))
      result%structural_conditional_vcov=0.5_dp*(v_structural_conditional+transpose(v_structural_conditional))
      result%cross_vcov=matmul(result%measurement_vcov,transpose(jac))
      result%structural_marginal_vcov=result%structural_conditional_vcov + &
         matmul(jac,matmul(result%measurement_vcov,transpose(jac)))
      allocate(result%joint_vcov(k1+k2,k1+k2))
      result%joint_vcov=0.0_dp
      result%joint_vcov(1:k1,1:k1)=result%measurement_vcov
      result%joint_vcov(1:k1,k1+1:k1+k2)=result%cross_vcov
      result%joint_vcov(k1+1:k1+k2,1:k1)=transpose(result%cross_vcov)
      result%joint_vcov(k1+1:k1+k2,k1+1:k1+k2)=result%structural_marginal_vcov
      result%joint_vcov=0.5_dp*(result%joint_vcov+transpose(result%joint_vcov))
      result%status=0
   end subroutine sam_block_covariance

   subroutine sam_second_order_bias(stage1_bias,stage2_conditional_bias,jac,curvature,v_stage1,stage2_bias,status)
      real(dp),intent(in)::stage1_bias(:),stage2_conditional_bias(:),jac(:,:),curvature(:,:,:),v_stage1(:,:)
      real(dp),allocatable,intent(out)::stage2_bias(:)
      integer,intent(out)::status
      integer::k1,k2,a,i,j
      real(dp)::quad
      k1=size(stage1_bias)
      k2=size(stage2_conditional_bias)
      if(any(shape(jac)/=[k2,k1]) .or. any(shape(curvature)/=[k2,k1,k1]) .or. &
         any(shape(v_stage1)/=[k1,k1])) then
         status=-1
         allocate(stage2_bias(0))
         return
      end if
      stage2_bias=stage2_conditional_bias+matmul(jac,stage1_bias)
      do a=1,k2
         quad=0.0_dp
         do j=1,k1
            do i=1,k1
               quad=quad+curvature(a,i,j)*v_stage1(i,j)
            end do
         end do
         stage2_bias(a)=stage2_bias(a)+0.5_dp*quad
      end do
      status=0
   end subroutine sam_second_order_bias

end module lavaan_sam_blocks
