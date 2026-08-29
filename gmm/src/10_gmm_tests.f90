! SPDX-License-Identifier: GPL-2.0-or-later
module gmm_tests
use r_compat, only: dp,pchisq
use gmm_linalg, only: invert_matrix, center_columns
use gmm_covariance, only: moment_covariance, hac_covariance
implicit none
private
public :: kleibergen_result_t, kleibergen_k_from_blocks, kleibergen_k, j_test
public :: gmm_moment_jacobian_function

abstract interface
   pure function gmm_moment_jacobian_function(theta,data) result(jac)
      import :: dp
      real(dp),intent(in)::theta(:),data(:,:)
      real(dp),allocatable::jac(:,:,:) ! n x q x k
   end function gmm_moment_jacobian_function
end interface

type :: kleibergen_result_t
   real(dp)::k_stat=0,j_stat=0,s_stat=0,k_pvalue=1,j_pvalue=1,s_pvalue=1
   integer::df_k=0,df_j=0,df_s=0
   real(dp),allocatable::d(:,:),bread(:),meat(:,:)
end type
contains
subroutine j_test(objective,n,df,stat,pvalue)
real(dp),intent(in)::objective
integer,intent(in)::n,df
real(dp),intent(out)::stat,pvalue
stat=objective*real(n,dp)
if(df>0)then
pvalue=1-pchisq(stat,real(df,dp))
else
pvalue=1.0_dp
end if
end subroutine


subroutine kleibergen_k(moment,jacobian,theta0,data,df_k,res,use_hac,bandwidth,kernel)
use gmm_estimation, only: gmm_moment_function
procedure(gmm_moment_function)::moment
procedure(gmm_moment_jacobian_function)::jacobian
real(dp),intent(in)::theta0(:),data(:,:)
integer,intent(in)::df_k
type(kleibergen_result_t),intent(out)::res
logical,intent(in),optional::use_hac
real(dp),intent(in),optional::bandwidth
character(len=*),intent(in),optional::kernel
real(dp),allocatable::gt(:,:),jac(:,:,:),allm(:,:),v(:,:),vff(:,:),vtf(:,:),qt(:,:)
real(dp)::bw
character(len=32)::kern
logical::hac
integer::n,q,k,j,c1,c2
gt=moment(theta0,data)
jac=jacobian(theta0,data)
n=size(gt,1)
q=size(gt,2)
k=size(theta0)
allocate(allm(n,q+q*k))
allm(:,1:q)=gt
do j=1,k
 c1=q+(j-1)*q+1
 c2=q+j*q
 allm(:,c1:c2)=jac(:,:,j)
end do
allm=center_columns(allm)
hac=.false.
if(present(use_hac))hac=use_hac
kern='Quadratic Spectral'
if(present(kernel))kern=kernel
bw=0.0_dp
if(present(bandwidth))bw=bandwidth
if(hac)then
v=hac_covariance(allm,bw,kern,.false.)
else
v=moment_covariance(allm,.false.)
end if
vff=v(1:q,1:q)
vtf=v(q+1:,1:q)
allocate(qt(q,k))
do j=1,k
qt(:,j)=sum(jac(:,:,j),dim=1)
end do
call kleibergen_k_from_blocks(sum(gt,dim=1),qt,vff,vtf,n,df_k,res)
end subroutine kleibergen_k

subroutine kleibergen_k_from_blocks(f_t,q_t,v_ff,v_theta_f,n,df_k,res)
! Direct computational core of KTest after .BigCov has formed its covariance blocks.
real(dp),intent(in)::f_t(:),q_t(:,:),v_ff(:,:),v_theta_f(:,:)
integer,intent(in)::n,df_k
type(kleibergen_result_t),intent(out)::res
real(dp),allocatable::vfi(:,:),adj(:),dvec(:),bread(:),meati(:,:)
integer::info,q,k
q=size(f_t)
k=size(q_t,2)
call invert_matrix(v_ff,vfi,info)
if(info/=0)then
res%k_stat=huge(1.0_dp)
return
end if
adj=matmul(v_theta_f,matmul(vfi,f_t))
allocate(dvec(q*k))
dvec=reshape(q_t,[q*k])-adj
allocate(res%d(q,k))
res%d=reshape(dvec,[q,k])
res%meat=matmul(transpose(res%d),matmul(vfi,res%d))
bread=matmul(transpose(res%d),matmul(vfi,f_t))
res%bread=bread
call invert_matrix(res%meat,meati,info)
if(info/=0)then
res%k_stat=huge(1.0_dp)
return
end if
res%k_stat=dot_product(bread,matmul(meati,bread))/real(n,dp)
res%s_stat=dot_product(f_t,matmul(vfi,f_t))/real(n,dp)
res%j_stat=res%s_stat-res%k_stat
res%df_k=df_k
res%df_j=q-k
res%df_s=res%df_k+res%df_j
if(res%df_k>0)res%k_pvalue=1-pchisq(res%k_stat,real(res%df_k,dp))
if(res%df_j>0)res%j_pvalue=1-pchisq(res%j_stat,real(res%df_j,dp))
if(res%df_s>0)res%s_pvalue=1-pchisq(res%s_stat,real(res%df_s,dp))
end subroutine
end module gmm_tests
