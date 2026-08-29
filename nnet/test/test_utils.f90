program test_utils
use nnet, only: dp,class_ind,summarize_rows
implicit none
integer::lab(5)
real(dp)::x(5,2),y(5,1)
real(dp),allocatable::ci(:,:),xo(:,:),yo(:,:)
lab=[1,3,2,3,1]
ci=class_ind(lab)
if(any(abs(sum(ci,dim=2)-1.0_dp)>0.0_dp)) error stop 'class_ind'
x=reshape([1._dp,2._dp,1._dp,2._dp,1._dp, 0._dp,0._dp,0._dp,1._dp,0._dp],[5,2])
y(:,1)=[1._dp,2._dp,3._dp,4._dp,5._dp]
call summarize_rows(x,y,xo,yo)
if(size(xo,1)/=3) error stop 'summarize row count'
if(abs(sum(yo)-15.0_dp)>1e-14_dp) error stop 'summarize total'
print *,'test_utils passed'
end program
