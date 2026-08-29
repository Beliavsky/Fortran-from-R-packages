program test_conditional_full
   use SpatialExtremes
   implicit none
   real(dp)::cov(3,3),y(2),coord(3,2),p0
   integer,allocatable::parts(:,:),chain(:,:),hit(:,:)
   real(dp),allocatable::weights(:),x(:,:)
   type(conditional_maxstable_result_t)::fit
   integer::info

   cov=reshape([1.0_dp,0.4_dp,0.3_dp, 0.4_dp,1.0_dp,0.2_dp, 0.3_dp,0.2_dp,1.0_dp],[3,3])
   y=[1.0_dp,1.2_dp]
   coord=reshape([0.0_dp,0.0_dp, 1.0_dp,0.3_dp, 0.2_dp,1.4_dp],[3,2],order=[2,1])

   p0=mvnorm_cdf_qmc([0.0_dp],reshape([1.0_dp],[1,1]),n_per_dim=100,info=info)
   if(info/=0 .or. abs(p0-0.5_dp)>0.03_dp)error stop 'mvnorm QMC regression'
   p0=mvstudent_cdf_qmc([0.0_dp],3.0_dp,[0.0_dp],reshape([1.0_dp],[1,1]),n_per_dim=100,info=info)
   if(info/=0 .or. abs(p0-0.5_dp)>0.03_dp)error stop 'mvstudent QMC regression'

   call schlather_partition_weights(cov(1:2,1:2),y,parts,weights,info,n_per_dim=75)
   if(info/=0 .or. size(weights)/=2 .or. abs(sum(weights)-1.0_dp)>1.0e-10_dp) &
      error stop 'Schlather partition weights'
   call extremalt_partition_weights(cov(1:2,1:2),y,2.0_dp,parts,weights,info,n_per_dim=75)
   if(info/=0 .or. abs(sum(weights)-1.0_dp)>1.0e-10_dp)error stop 'extremal-t partition weights'
   call brownresnick_partition_weights(coord(1:2,:),y,2.0_dp,1.2_dp,parts,weights,info,n_per_dim=75)
   if(info/=0 .or. abs(sum(weights)-1.0_dp)>1.0e-10_dp)error stop 'Brown-Resnick partition weights'

   call gibbs_partitions_schlather(cov(1:2,1:2),y,3,1,1,chain,info,n_per_dim=40)
   if(info/=0 .or. any(shape(chain)/=[3,2]))error stop 'Schlather partition Gibbs'
   call gibbs_partitions_extremalt(cov(1:2,1:2),y,2.0_dp,3,1,1,chain,info,n_per_dim=40)
   if(info/=0 .or. any(shape(chain)/=[3,2]))error stop 'extremal-t partition Gibbs'
   call gibbs_partitions_brownresnick(coord(1:2,:),y,2.0_dp,1.2_dp,3,1,1,chain,info,n_per_dim=40)
   if(info/=0 .or. any(shape(chain)/=[3,2]))error stop 'Brown-Resnick partition Gibbs'

   call sample_conditional_schlather(cov,y,2,fit,n_per_dim=60)
   if(fit%info/=0 .or. maxval(abs(fit%sim(:,1:2)-spread(y,1,2)))>1.0e-10_dp) &
      error stop 'automatic conditional Schlather'
   call sample_conditional_extremalt(cov,y,2.0_dp,2,fit,n_per_dim=60)
   if(fit%info/=0 .or. maxval(abs(fit%sim(:,1:2)-spread(y,1,2)))>1.0e-10_dp) &
      error stop 'automatic conditional extremal-t'
   call sample_conditional_brownresnick(coord,y,2.0_dp,1.2_dp,2,fit,n_per_dim=60,n_subextremal=20)
   if(fit%info/=0 .or. maxval(abs(fit%sim(:,1:2)-spread(y,1,2)))>1.0e-10_dp) &
      error stop 'automatic conditional Brown-Resnick'

   call starting_partitions_schlather(cov,3,parts,info)
   if(info/=0 .or. any(shape(parts)/=[3,3]))error stop 'Schlather starting partitions'
   call starting_partitions_extremalt(cov,2.0_dp,3,parts,info)
   if(info/=0 .or. any(shape(parts)/=[3,3]))error stop 'extremal-t starting partitions'
   call starting_partitions_brownresnick(coord,2.0_dp,1.2_dp,3,parts,info)
   if(info/=0 .or. any(shape(parts)/=[3,3]))error stop 'Brown-Resnick starting partitions'

   call simulate_brownresnick_exact_hitting(2,coord,2.0_dp,1.2_dp,x,hit,info)
   if(info/=0 .or. any(shape(x)/=[2,3]) .or. any(hit<1))error stop 'Brown-Resnick hitting scenario'
   print '(a)','full conditional parity tests passed'
end program test_conditional_full
