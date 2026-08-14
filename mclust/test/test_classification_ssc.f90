program test_classification_ssc
  use mclust
  implicit none
  integer,parameter::n=120,d=2
  real(dp)::x(n,d)
  integer::truth(n),semi(n),st,i
  type(mclust_da_fit)::da
  type(mclust_fit)::ssc
  real(dp),allocatable::post(:,:)
  integer,allocatable::pred(:)
  integer,parameter::gs(1)=[1]
  character(len=3),parameter::mods(2)=[character(len=3)::'EII','EEE']
  do i=1,n/2
    truth(i)=1; x(i,1)=-3.0_dp+0.3_dp*sin(0.4_dp*i); x(i,2)=-2.0_dp+0.2_dp*cos(0.7_dp*i)
  end do
  do i=n/2+1,n
    truth(i)=2; x(i,1)=3.5_dp+0.3_dp*sin(0.5_dp*i); x(i,2)=2.5_dp+0.3_dp*cos(0.8_dp*i)
  end do
  call fit_mclust_da(x,truth,da,g_values=gs,model_names=mods,status=st)
  if(st/=0)error stop 'DA fit'
  call predict_mclust_da(da,x,post,pred,status=st); if(st/=0)error stop 'DA predict'
  if(class_error_rate(pred,truth)>0.01_dp)error stop 'DA error'
  semi=0; semi(1:10)=1; semi(n/2+1:n/2+10)=2
  call mclust_ssc_select(x,semi,2,ssc,model_names=mods,status=st)
  if(st/=0)error stop 'SSC fit'
  if(class_error_rate(ssc%classification,truth)>0.01_dp) then
    print *, 'ssc error',class_error_rate(ssc%classification,truth)
    error stop 'SSC error'
  end if
  if(any(ssc%classification(1:10)/=1) .or. any(ssc%classification(n/2+1:n/2+10)/=2)) error stop 'known classes changed'
  print *, 'test_classification_ssc PASS '
end program test_classification_ssc
