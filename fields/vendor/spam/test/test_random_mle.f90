program test_random_mle
use spam
use r_mod,only:set_seed_int
implicit none
integer,parameter::n=5000
real(dp)::sig(2,2),mu(2),m(2),s(2,2),loc(5,1),y(5),xx(5,1)
real(dp),allocatable::z(:,:)
type(csr_matrix)::a,dist,q;type(mle_result)::fit
integer::i
sig=reshape([1d0,0.4d0,0.4d0,2d0],[2,2]);a=csr_from_dense(sig);mu=[1d0,-2d0]
call set_seed_int(12345);z=rmvnorm_cov(n,mu,a);m=sum(z,dim=1)/real(n,dp)
s=0;do i=1,n;s=s+outer(z(i,:)-m,z(i,:)-m);end do;s=s/real(n,dp)
call check(maxval(abs(m-mu))<0.07d0,'rmvnorm mean');call check(maxval(abs(s-sig))<0.1d0,'rmvnorm covariance')
q=csr_from_dense(reshape([2d0,0.5d0,0.5d0,1.5d0],[2,2]))
call set_seed_int(777);z=rmvnorm_prec(n,[0d0,0d0],q);m=sum(z,dim=1)/real(n,dp)
s=0;do i=1,n;s=s+outer(z(i,:)-m,z(i,:)-m);end do;s=s/real(n,dp)
call check(maxval(abs(matmul(csr_to_dense(q),s)-reshape([1d0,0d0,0d0,1d0],[2,2])))<0.12d0,'rmvnorm.prec')
z=rmvnorm_prec_const(100,[0d0,0d0],q,reshape([1d0,1d0],[1,2]),[2d0])
call check(maxval(abs(sum(z,dim=2)-2d0))<1d-11,'rmvnorm.prec.const')
loc(:,1)=[0d0,1d0,2d0,3d0,4d0];dist=nearest_dist(loc,delta=10d0,full=.true.);y=[1d0,2d0,3d0,4d0,5d0];xx(:,1)=1d0
fit=mle_spam(y,xx,dist,[2d0],[1d0],[0.01d0],[20d0],'nugget',maxit=100,reltol=1d-8)
call check(fit%convergence==0,'mle convergence');call check(abs(fit%par(1)-3d0)<2d-3,'mle mean')
call check(abs(fit%par(2)-2d0)<3d-3,'mle variance')
print *,'test_random_mle: PASS'
contains
subroutine check(ok,msg);logical,intent(in)::ok;character(*),intent(in)::msg;if(.not.ok)then;print *,msg
error stop;end if;end
pure function outer(a,b) result(c);real(dp),intent(in)::a(:),b(:);real(dp)::c(size(a),size(b));integer::j
do j=1,size(a);c(j,:)=a(j)*b;end do;end
end program
