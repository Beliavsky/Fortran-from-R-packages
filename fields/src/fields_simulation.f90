! GPL-2.0-or-later. Random-field/data simulation corresponding to fields sim.* routines.
module fields_simulation
use fields_kinds, only: dp
use fields_covariance, only: stationary_covariance, paciorek_covariance
use fields_linalg, only: mvn_sample
use fields_kriging, only: krig_fit, krig_simulate
use r_mod, only: rnorm1
implicit none
private
public :: simulate_random_field, simulate_paciorek_field, simulate_spatial_data, conditional_field_simulation

contains

function simulate_random_field(x,nsim,model,a_range,smoothness,power,sigma2,nugget,info) result(draws)
real(dp),intent(in)::x(:,:)
integer,intent(in)::nsim
character(len=*),intent(in),optional::model
real(dp),intent(in),optional::a_range,smoothness,power,sigma2,nugget
integer,intent(out),optional::info
real(dp),allocatable::draws(:,:),k(:,:),mu(:)
real(dp)::s2,t2
integer::i,ierr
s2=1.0_dp;if(present(sigma2))s2=sigma2
t2=0.0_dp;if(present(nugget))t2=nugget
k=stationary_covariance(x,x,model,a_range,smoothness,power,phi=s2)
do i=1,size(k,1);k(i,i)=k(i,i)+t2;end do
allocate(mu(size(x,1)));mu=0.0_dp;draws=mvn_sample(mu,k,nsim,ierr);if(present(info))info=ierr
end function simulate_random_field

function simulate_paciorek_field(x,local_range,nsim,smoothness,sigma2,info) result(draws)
real(dp),intent(in)::x(:,:),local_range(:)
integer,intent(in)::nsim
real(dp),intent(in),optional::smoothness,sigma2
integer,intent(out),optional::info
real(dp),allocatable::draws(:,:),k(:,:),mu(:)
real(dp)::s2
integer::ierr
s2=1.0_dp;if(present(sigma2))s2=sigma2
k=s2*paciorek_covariance(x,x,local_range,local_range,smoothness)
allocate(mu(size(x,1)));mu=0.0_dp;draws=mvn_sample(mu,k,nsim,ierr);if(present(info))info=ierr
end function simulate_paciorek_field

function simulate_spatial_data(x,nsim,model,a_range,smoothness,power,sigma2,tau2,trend,info) result(y)
real(dp),intent(in)::x(:,:)
integer,intent(in)::nsim
character(len=*),intent(in),optional::model
real(dp),intent(in),optional::a_range,smoothness,power,sigma2,tau2,trend(:)
integer,intent(out),optional::info
real(dp),allocatable::y(:,:),f(:,:)
real(dp)::t2
integer::i,j,ierr
f=simulate_random_field(x,nsim,model,a_range,smoothness,power,sigma2,info=ierr)
y=f;t2=0.0_dp;if(present(tau2))t2=tau2
if(present(trend))then;if(size(trend)/=size(x,1))error stop 'simulate_spatial_data: trend mismatch';do j=1,nsim;y(:,j)=y(:,j)+trend;end do;end if
if(t2>0.0_dp)then;do j=1,nsim;do i=1,size(x,1);y(i,j)=y(i,j)+sqrt(t2)*rnorm1();end do;end do;end if
if(present(info))info=ierr
end function simulate_spatial_data

function conditional_field_simulation(fit,k_new,k_newnew,nsim,t_new,include_nugget,new_weights,info) result(draws)
type(krig_fit),intent(in)::fit
real(dp),intent(in)::k_new(:,:),k_newnew(:,:)
integer,intent(in)::nsim
real(dp),intent(in),optional::t_new(:,:),new_weights(:)
logical,intent(in),optional::include_nugget
integer,intent(out),optional::info
real(dp),allocatable::draws(:,:)
draws=krig_simulate(fit,k_new,k_newnew,nsim,t_new,include_nugget,new_weights,info)
end function conditional_field_simulation

end module fields_simulation
