program test_mapping_crimcoords
  use mclust
  implicit none
  integer,parameter::n=12,d=2
  real(dp)::x(n,d)
  integer::class(n),perm(n),mapped(n),vote(n,3),maj(n),st,i
  integer,allocatable::mapping(:,:)
  real(dp),allocatable::z(:,:)
  type(crimcoords_fit)::cc
  do i=1,4
    x(i,:)=[-2.0_dp+0.1_dp*i,0.2_dp*i]; class(i)=10
    x(i+4,:)=[0.1_dp*i,2.5_dp+0.1_dp*i]; class(i+4)=20
    x(i+8,:)=[2.5_dp+0.1_dp*i,-1.5_dp+0.05_dp*i]; class(i+8)=30
  end do
  where(class==10) perm=2
  where(class==20) perm=3
  where(class==30) perm=1
  call match_clusters(class,perm,mapped,mapping,st)
  if(st/=0 .or. any(mapped/=class)) error stop 'match_clusters'
  call unmap_classes(class,z)
  if(size(z,2)/=3 .or. any(abs(sum(z,dim=2)-1.0_dp)>epsilon(1.0_dp))) error stop 'unmap'
  vote(:,1)=class; vote(:,2)=mapped; vote(:,3)=class
  call majority_vote(vote,maj)
  if(any(maj/=class)) error stop 'majority_vote'
  call fit_crimcoords(x,class,cc,status=st)
  if(st/=0 .or. cc%numdir<1) error stop 'crimcoords status'
  if(cc%evalues(1)<=0.0_dp) error stop 'crimcoords eigenvalue'
  if(size(cc%projection,1)/=n) error stop 'crimcoords projection'
  print *, 'test_mapping_crimcoords PASS',cc%evalues(1)
end program test_mapping_crimcoords
