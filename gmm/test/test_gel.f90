program test_gel
use gmm, only: dp,gel_fit_result_t,gel_fit,GEL_EL,GEL_ET,GEL_CUE,GEL_HD,GEL_ETEL,GEL_ETHD
implicit none
real(dp)::data(8,2),theta0(1)
type(gel_fit_result_t)::r
integer::i
data(:,1)=[(real(i,dp),i=1,8)]
data(:,2)=[0.5_dp,1.2_dp,0.7_dp,1.5_dp,0.9_dp,1.8_dp,1.1_dp,2.0_dp]
theta0=4.5_dp
call chk(GEL_EL,5.07303379_dp,2e-4_dp)
call chk(GEL_ET,4.57388202_dp,2e-4_dp)
call chk(GEL_HD,4.82242889_dp,3e-4_dp)
call chk(GEL_CUE,3.69096154_dp,3e-4_dp)
call chk(GEL_ETEL,4.83733867_dp,4e-4_dp)
call chk(GEL_ETHD,4.75781925_dp,4e-4_dp)
print '(a)','test_gel: ok'
contains
subroutine chk(tp,expected,tol)
integer,intent(in)::tp
real(dp),intent(in)::expected,tol
call gel_fit(mom,data,theta0,tp,r,maxit=800,maxiterlam=200,tol=1e-10_dp)
if(abs(r%coefficients(1)-expected)>tol)then
 print *,'GEL type',tp,'got',r%coefficients(1),'expected',expected,'obj',r%objective,'conv',r%conv_par,r%conv_lambda
 error stop 1
end if
if(abs(sum(r%prob)-1.0_dp)>1e-10_dp .or. minval(r%prob)<-1e-12_dp)error stop 'bad implied probabilities'
end subroutine
pure function mom(theta,d) result(gt)
real(dp),intent(in)::theta(:),d(:,:)
real(dp),allocatable::gt(:,:)
allocate(gt(size(d,1),2))
gt(:,1)=d(:,1)-theta(1)
gt(:,2)=(d(:,1)-theta(1))*d(:,2)
end function
end program
