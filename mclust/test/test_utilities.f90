program test_utilities
  use mclust
  implicit none
  integer,parameter::n=90,d=3
  real(dp)::x(n,d),w(n)
  logical::missing(n,d)
  real(dp),allocatable::ximp(:,:),zcomb(:,:)
  type(mclust_fit)::fit,wfit
  type(mclust_dr_fit)::dr
  type(cluster_combination)::cc
  integer::i,st
  do i=1,n/3
    x(i,:)=[-3.0_dp+0.2_dp*sin(0.3_dp*i),-2.0_dp+0.2_dp*cos(0.4_dp*i),-1.0_dp+0.15_dp*sin(0.7_dp*i)]
  end do
  do i=n/3+1,2*n/3
    x(i,:)=[0.0_dp+0.2_dp*sin(0.5_dp*i),3.0_dp+0.2_dp*cos(0.6_dp*i),2.0_dp+0.15_dp*sin(0.8_dp*i)]
  end do
  do i=2*n/3+1,n
    x(i,:)=[4.0_dp+0.2_dp*sin(0.2_dp*i),-1.0_dp+0.2_dp*cos(0.3_dp*i),3.0_dp+0.15_dp*sin(0.9_dp*i)]
  end do
  call fit_model(x,3,'VVV',fit); if(fit%status<0)error stop 'fit'
  w=1.0_dp; w(1:n:3)=0.5_dp
  call fit_model_weighted(x,3,'VVV',w,wfit); if(wfit%status<0)error stop 'weighted fit'
  call fit_mclust_dr(fit,x,dr,lambda=0.5_dp,status=st); if(st/=0 .or. dr%numdir<1)error stop 'DR'
  missing=.false.; missing(2,1)=.true.; missing(17,2:3)=.true.; missing(41,:)=.true.
  call impute_data(fit,x,missing,ximp,st); if(st/=0)error stop 'impute'
  if(any(.not.(ximp<huge(1.0_dp))))error stop 'impute finite'
  call clust_combi(fit%z,cc); if(size(cc%merge,1)/=2)error stop 'combine tree'
  call apply_combination(fit%z,cc%merge,2,zcomb); if(size(zcomb,2)/=2)error stop 'combine result'
  if(maxval(abs(sum(zcomb,dim=2)-1.0_dp))>1e-10_dp)error stop 'combine row sums'
  if(adjusted_rand_index(fit%classification,fit%classification)<0.999999_dp)error stop 'ARI'
  print *, 'test_utilities PASS ',dr%evalues(1)
end program test_utilities
