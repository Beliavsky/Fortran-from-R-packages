program test_random_algorithms
   use, intrinsic :: iso_fortran_env, only : int64
   use rfast
   implicit none
   real(dp) :: mu(2), sigma(2,2), sample(20000,2), m(2)
   real(dp) :: train(4,1), test(2,1), w(4,4), fw(4,4), pts(4,2), sm(2)
   integer :: cls(4), pred(2), perm(3), ord(4)
   logical :: ok, adj(4,4)

   call set_seed(12345_int64)
   mu=[1.0_dp,-2.0_dp]
   sigma=reshape([1.0_dp,0.3_dp,0.3_dp,2.0_dp],[2,2])
   sample=rmvnorm(size(sample,1),mu,sigma)
   m=sum(sample,dim=1)/real(size(sample,1),dp)
   call assert_true(maxval(abs(m-mu))<0.05_dp,'mvnormal mean')

   train(:,1)=[0.0_dp,1.0_dp,10.0_dp,11.0_dp]
   cls=[1,1,2,2]
   test(:,1)=[0.2_dp,10.2_dp]
   pred=knn_classify(train,cls,test,1)
   call assert_true(all(pred==[1,2]),'knn')

   w=huge(1.0_dp)/8.0_dp
   w(1,1)=0.0_dp;w(2,2)=0.0_dp;w(3,3)=0.0_dp;w(4,4)=0.0_dp
   w(1,2)=2.0_dp;w(2,3)=3.0_dp;w(1,4)=10.0_dp;w(3,4)=1.0_dp
   fw=floyd_warshall(w)
   call assert_close(fw(1,4),6.0_dp,1e-12_dp,'floyd')

   adj=.false.;adj(1,2)=.true.;adj(1,3)=.true.;adj(2,4)=.true.;adj(3,4)=.true.
   ord=topological_sort(adj,ok)
   call assert_true(ok .and. ord(1)==1 .and. ord(4)==4,'topological')
   perm=[1,2,3];ok=next_permutation(perm)
   call assert_true(ok .and. all(perm==[1,3,2]),'next permutation')

   pts=reshape([0.0_dp,2.0_dp,0.0_dp,2.0_dp,0.0_dp,0.0_dp,2.0_dp,2.0_dp],[4,2])
   sm=spatial_median(pts)
   call assert_true(maxval(abs(sm-[1.0_dp,1.0_dp]))<1e-8_dp,'spatial median')

   print *, 'test_random_algorithms: PASS'
contains
   subroutine assert_close(got,want,tol,msg)
      real(dp),intent(in)::got,want,tol
      character(*),intent(in)::msg
      if(abs(got-want)>tol)then
         print *, 'FAIL ',trim(msg),got,want
         error stop 1
      end if
   end subroutine
   subroutine assert_true(good,msg)
      logical,intent(in)::good
      character(*),intent(in)::msg
      if(.not.good)then
         print *, 'FAIL ',trim(msg)
         error stop 1
      end if
   end subroutine
end program test_random_algorithms
