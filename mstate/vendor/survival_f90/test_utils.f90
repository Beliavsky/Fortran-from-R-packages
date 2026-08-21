program test_utils
  use survival
  implicit none
  real(dp)::s(1)=[0._dp],e(1)=[2.5_dp],ct(4)=[1._dp,2._dp,3._dp,4._dp]
  real(dp)::cp(4)=[0.9_dp,0.8_dp,0.7_dp,0.6_dp]
  logical::ext(1)=[.true.],keep(4)=[.true.,.true.,.true.,.true.]
  integer::stat(1)=[2]
  integer,allocatable::row(:),add(:),os(:),rid(:)
  real(dp),allocatable::a(:),b(:),wt(:),ss(:),ee(:)
  call finegray_expand(s,e,ct,cp,ext,keep,row,a,b,wt,add)
  if(size(row)/=2) error stop 'finegray size'
  if(abs(wt(2)-cp(4)/cp(3))>1e-12_dp) error stop 'finegray weight'
  call surv_split(s,e,stat,[1._dp,2._dp],ss,ee,os,rid)
  if(size(ss)/=3) error stop 'split size'
  print *, 'test_utils PASS'
end program
