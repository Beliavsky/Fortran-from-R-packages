module spam_grf
use spam_kinds,only:dp
use spam_types,only:csr_matrix
use spam_distance,only:nearest_dist
use spam_likelihood,only:make_covariance
use spam_random,only:rmvnorm_cov
use r_mod,only:runif_vec
implicit none
private
public::rgrf,rgrf_grid
contains
function rgrf(n,locs,theta,model,beta,xdesign,compact) result(samples)
integer,intent(in)::n;real(dp),intent(in)::locs(:,:),theta(:);character(len=*),intent(in)::model
real(dp),intent(in),optional::beta(:),xdesign(:,:);logical,intent(in),optional::compact
real(dp),allocatable::samples(:,:),mu(:);type(csr_matrix)::dist,sigma
real(dp)::del;logical::cp;integer::nl
nl=size(locs,1);allocate(mu(nl));mu=0.0_dp
if(present(beta))then
 if(size(beta)==1)then
  if(present(xdesign))then
   if(size(xdesign,1)/=nl)error stop 'rgrf: X row mismatch';mu=xdesign(:,1)*beta(1)
  else;mu=beta(1);end if
 else
  if(.not.present(xdesign))error stop 'rgrf: X required for vector beta'
  if(size(xdesign,1)/=nl.or.size(xdesign,2)/=size(beta))error stop 'rgrf: X/beta mismatch'
  mu=matmul(xdesign,beta)
 end if
end if
cp=.false.;if(present(compact))cp=compact
if(cp .and. size(theta)>0)then;del=abs(theta(1));else;del=sqrt(huge(1.0_dp))*0.25_dp;end if
dist=nearest_dist(locs,y=locs,method='euclidean',delta=del)
sigma=make_covariance(dist,theta,model);samples=rmvnorm_cov(n,mu,sigma)
end function

function rgrf_grid(nsamp,nx,ny,theta,model,beta,xlim,ylim,tau,compact,locs_out) result(samples)
integer,intent(in)::nsamp,nx,ny;real(dp),intent(in)::theta(:);character(len=*),intent(in)::model
real(dp),intent(in),optional::beta,xlim(2),ylim(2),tau;logical,intent(in),optional::compact
real(dp),allocatable,intent(out),optional::locs_out(:,:);real(dp),allocatable::samples(:,:),locs(:,:),rv(:)
real(dp)::xl(2),yl(2),tt,dx,dy,bv;integer::i,j,k
if(nx<1.or.ny<1)error stop 'rgrf_grid: nx,ny must be positive'
xl=[0.0_dp,1.0_dp];yl=xl;if(present(xlim))xl=xlim;if(present(ylim))yl=ylim
tt=0.0_dp;if(present(tau))tt=tau;if(tt<0.0_dp.or.tt>=0.5_dp)error stop 'rgrf_grid: tau must be in [0,0.5)'
dx=(xl(2)-xl(1))/real(nx,dp);dy=(yl(2)-yl(1))/real(ny,dp);allocate(locs(nx*ny,2));k=0
do j=1,ny;do i=1,nx;k=k+1;locs(k,:)=[xl(1)+(real(i,dp)-0.5_dp)*dx,yl(1)+(real(j,dp)-0.5_dp)*dy];end do;end do
if(tt>0.0_dp)then
 rv=runif_vec(2*nx*ny);do k=1,nx*ny;locs(k,1)=locs(k,1)+(2*rv(2*k-1)-1)*tt*dx;locs(k,2)=locs(k,2)+(2*rv(2*k)-1)*tt*dy;end do
end if
bv=0.0_dp;if(present(beta))bv=beta
samples=rgrf(nsamp,locs,theta,model,beta=[bv],compact=compact)
if(present(locs_out))then;allocate(locs_out(size(locs,1),2));locs_out=locs;end if
end function
end module spam_grf
