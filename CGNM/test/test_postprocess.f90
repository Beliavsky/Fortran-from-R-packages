program test_postprocess
   use cgnm, only : dp,cgnm_result,column_quantiles,accepted_indices,best_approximate_minimizers
   implicit none
   type(cgnm_result) :: res
   real(dp) :: q(2)
   real(dp),allocatable :: theta(:,:)
   integer,allocatable :: idx(:)
   allocate(res%residual_history(4,1),res%theta(4,2),res%y(4,2))
   res%residual_history(:,1)=[1._dp,1.1_dp,1.2_dp,100._dp]
   res%theta=reshape([1._dp,2._dp,3._dp,4._dp, 5._dp,6._dp,7._dp,8._dp],[4,2])
   res%y=0._dp
   call column_quantiles(res%theta,0.5_dp,q)
   if(maxval(abs(q-[2.5_dp,6.5_dp]))>1.e-12_dp) error stop 'quantile'
   call best_approximate_minimizers(res,2,theta,idx)
   if(any(idx/=[1,2])) error stop 'best'
   call accepted_indices(res,idx,use_accepted=.false.,num_parameters=2)
   if(size(idx)/=2) error stop 'accepted fixed count'
end program test_postprocess
