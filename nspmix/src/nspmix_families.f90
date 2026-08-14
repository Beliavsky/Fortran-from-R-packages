module nspmix_families
   use nspmix_kinds, only : dp
   use nspmix_types
   implicit none
   private
   public :: make_npnorm_data, make_nppois_data, make_npgeom_data, make_npnbinom_data
   public :: make_cvps_data, make_mlogit_data, nobs_data, data_weights, support_bounds
   public :: gridpoints, logd_eval, default_beta, default_mix_points, valid_parameters
contains
   subroutine make_npnorm_data(v,data,w)
      real(dp), intent(in) :: v(:)
      type(nsp_data), intent(out) :: data
      real(dp), intent(in), optional :: w(:)
      data%family=NSP_NORMAL; data%v=v; allocate(data%w(size(v))); data%w=1.0_dp
      if(present(w)) then
         if(size(w)==1) then; data%w=w(1); else; data%w=w; end if
      end if
   end subroutine

   subroutine make_nppois_data(v,data,w)
      real(dp), intent(in) :: v(:)
      type(nsp_data), intent(out) :: data
      real(dp), intent(in), optional :: w(:)
      call make_simple(v,NSP_POISSON,data,w)
   end subroutine
   subroutine make_npgeom_data(v,data,w)
      real(dp), intent(in) :: v(:)
      type(nsp_data), intent(out) :: data
      real(dp), intent(in), optional :: w(:)
      call make_simple(v,NSP_GEOM,data,w)
   end subroutine
   subroutine make_npnbinom_data(v,size_par,data,w)
      real(dp), intent(in) :: v(:),size_par
      type(nsp_data), intent(out) :: data
      real(dp), intent(in), optional :: w(:)
      call make_simple(v,NSP_NBINOM,data,w); data%size=size_par
   end subroutine
   subroutine make_simple(v,fam,data,w)
      real(dp), intent(in) :: v(:); integer,intent(in)::fam
      type(nsp_data), intent(out)::data; real(dp),intent(in),optional::w(:)
      data%family=fam; data%v=v; allocate(data%w(size(v))); data%w=1.0_dp
      if(present(w)) then
         if(size(w)==1) then; data%w=w(1); else; data%w=w; end if
      end if
   end subroutine

   subroutine make_cvps_data(ni,mi,ri,data)
      real(dp), intent(in) :: ni(:),mi(:),ri(:)
      type(nsp_data), intent(out) :: data
      data%family=NSP_CVPS; data%ni=ni; data%mi=mi; data%ri=ri
   end subroutine

   subroutine make_mlogit_data(group,y,trials,xmat,data)
      integer, intent(in) :: group(:)
      real(dp), intent(in) :: y(:),trials(:),xmat(:,:)
      type(nsp_data), intent(out) :: data
      integer, allocatable :: uniq(:)
      integer :: i,j,k,ng
      if(size(group)/=size(y) .or. size(y)/=size(trials) .or. size(xmat,1)/=size(y)) error stop "mlogit dimensions"
      allocate(uniq(size(group)))
      data%family=NSP_MLOGIT; allocate(data%group(size(group)))
      ng=0; uniq=0
      do i=1,size(group)
         k=0
         do j=1,ng
            if(uniq(j)==group(i)) then; k=j; exit; end if
         end do
         if(k==0) then; ng=ng+1; uniq(ng)=group(i); k=ng; end if
         data%group(i)=k
      end do
      data%y=y; data%trials=trials; data%xmat=xmat
   end subroutine

   integer function nobs_data(data)
      type(nsp_data), intent(in) :: data
      select case(data%family)
      case(NSP_NORMAL,NSP_POISSON,NSP_GEOM,NSP_NBINOM); nobs_data=size(data%v)
      case(NSP_CVPS); nobs_data=size(data%ni)
      case(NSP_MLOGIT); nobs_data=maxval(data%group)
      case default; nobs_data=0
      end select
   end function

   subroutine data_weights(data,beta,w)
      type(nsp_data),intent(in)::data; real(dp),intent(in)::beta(:)
      real(dp),allocatable,intent(out)::w(:)
      integer :: n
      if(size(beta)<0) error stop "unreachable"
      n=nobs_data(data); allocate(w(n)); w=1.0_dp
      if(allocated(data%w)) w=data%w
   end subroutine

   subroutine support_bounds(data,beta,lo,hi)
      type(nsp_data),intent(in)::data; real(dp),intent(in)::beta(:)
      real(dp),intent(out)::lo,hi
      if(size(beta)<0) error stop "unreachable"
      lo=-huge(1.0_dp); hi=huge(1.0_dp)
      select case(data%family)
      case(NSP_POISSON); lo=0.0_dp
      case(NSP_GEOM,NSP_NBINOM); lo=0.0_dp; hi=1.0_dp
      end select
   end subroutine

   subroutine default_beta(data,beta)
      type(nsp_data),intent(in)::data; real(dp),allocatable,intent(out)::beta(:)
      integer :: q
      select case(data%family)
      case(NSP_NORMAL); allocate(beta(1)); beta=1.0_dp
      case(NSP_CVPS)
         allocate(beta(1)); beta(1)=sqrt(max(sum(data%ri)/max(sum(data%ni)-real(size(data%ni),dp),1.0_dp),1.0e-12_dp))
      case(NSP_MLOGIT); q=size(data%xmat,2); allocate(beta(q)); beta=0.0_dp
      case default; allocate(beta(0))
      end select
   end subroutine

   subroutine default_mix_points(data,beta,kmax,pt)
      type(nsp_data),intent(in)::data; real(dp),intent(in)::beta(:); integer,intent(in)::kmax
      real(dp),allocatable,intent(out)::pt(:)
      integer :: k
      k=max(1,min(kmax,20)); call gridpoints(data,beta,k,pt)
   end subroutine

   subroutine gridpoints(data,beta,ngrid,pt)
      type(nsp_data),intent(in)::data; real(dp),intent(in)::beta(:); integer,intent(in)::ngrid
      real(dp),allocatable,intent(out)::pt(:)
      real(dp)::a,b,t,pbar,xb,eta
      integer::i,g,ng
      allocate(pt(max(1,ngrid)))
      select case(data%family)
      case(NSP_NORMAL)
         a=minval(data%v); b=maxval(data%v)
      case(NSP_POISSON)
         a=sqrt(max(minval(data%v),0.0_dp)); b=sqrt(max(maxval(data%v),1.0_dp))
         do i=1,size(pt); t=real(i-1,dp)/real(max(size(pt)-1,1),dp); pt(i)=((1-t)*a+t*b)**2; end do; return
      case(NSP_GEOM)
         a=-log(1.0_dp+maxval(data%v)); b=-log(1.0_dp+minval(data%v))
         do i=1,size(pt); t=real(i-1,dp)/real(max(size(pt)-1,1),dp); pt(i)=exp((1-t)*a+t*b); end do; return
      case(NSP_NBINOM)
         a=log(data%size/(data%size+maxval(data%v))); b=log(data%size/(data%size+minval(data%v)))
         do i=1,size(pt); t=real(i-1,dp)/real(max(size(pt)-1,1),dp); pt(i)=exp((1-t)*a+t*b); end do; return
      case(NSP_CVPS)
         a=minval(data%mi); b=maxval(data%mi)
      case(NSP_MLOGIT)
         ng=nobs_data(data); a=huge(1.0_dp); b=-huge(1.0_dp)
         do g=1,ng
            pbar=sum(data%y,mask=data%group==g)/max(sum(data%trials,mask=data%group==g),1.0e-12_dp)
            pbar=min(max(pbar,1.0e-4_dp),0.9999_dp); eta=log(pbar/(1.0_dp-pbar))
            do i=1,size(data%group)
               if(data%group(i)==g) then
                  xb=dot_product(data%xmat(i,:),beta); a=min(a,eta-xb); b=max(b,eta-xb)
               end if
            end do
         end do
      case default
         a=0.0_dp; b=1.0_dp
      end select
      if(size(pt)==1) then; pt(1)=0.5_dp*(a+b); return; end if
      do i=1,size(pt); t=real(i-1,dp)/real(size(pt)-1,dp); pt(i)=(1-t)*a+t*b; end do
   end subroutine

   logical function valid_parameters(data,beta,pt)
      type(nsp_data),intent(in)::data; real(dp),intent(in)::beta(:),pt(:)
      real(dp)::lo,hi
      call support_bounds(data,beta,lo,hi)
      valid_parameters=all(pt>=lo) .and. all(pt<=hi)
      if(data%family==NSP_NORMAL .or. data%family==NSP_CVPS) then
         if(size(beta)/=1) then
            valid_parameters=.false.
         else
            valid_parameters=valid_parameters .and. beta(1)>0.0_dp
         end if
      end if
   end function

   subroutine logd_eval(data,beta,pt,ld,dt,db)
      type(nsp_data),intent(in)::data; real(dp),intent(in)::beta(:),pt(:)
      real(dp),allocatable,intent(out)::ld(:,:),dt(:,:),db(:,:,:)
      integer::n,k,i,j,g,q,r,ng
      real(dp)::th,z,b2,p,eta,resid,const
      n=nobs_data(data); k=size(pt); q=size(beta)
      allocate(ld(n,k),dt(n,k),db(n,k,q)); ld=0.0_dp; dt=0.0_dp; if(q>0) db=0.0_dp
      select case(data%family)
      case(NSP_NORMAL)
         b2=beta(1)*beta(1); const=-0.5_dp*log(2.0_dp*acos(-1.0_dp)*b2)
         do j=1,k; do i=1,n
            z=data%v(i)-pt(j); ld(i,j)=const-0.5_dp*z*z/b2; dt(i,j)=z/b2
            db(i,j,1)=(z*z/b2-1.0_dp)/beta(1)
         end do; end do
      case(NSP_POISSON)
         do j=1,k; th=max(pt(j),1.0e-100_dp); do i=1,n
            ld(i,j)=data%v(i)*log(th)-th-log_gamma(data%v(i)+1.0_dp); dt(i,j)=data%v(i)/th-1.0_dp
         end do; end do
      case(NSP_GEOM)
         do j=1,k; th=min(max(pt(j),1.0e-100_dp),1.0_dp-1.0e-10_dp); do i=1,n
            ld(i,j)=log(th)+data%v(i)*log(1.0_dp-th); dt(i,j)=1.0_dp/th-data%v(i)/(1.0_dp-th)
         end do; end do
      case(NSP_NBINOM)
         do j=1,k; th=min(max(pt(j),1.0e-100_dp),1.0_dp-1.0e-10_dp); do i=1,n
            ld(i,j)=log_gamma(data%size+data%v(i))-log_gamma(data%size)-log_gamma(data%v(i)+1.0_dp) &
                    +data%size*log(th)+data%v(i)*log(1.0_dp-th)
            dt(i,j)=data%size/th-data%v(i)/(1.0_dp-th)
         end do; end do
      case(NSP_CVPS)
         b2=beta(1)*beta(1)
         do j=1,k; do i=1,n
            z=data%mi(i)-pt(j); resid=data%ri(i)+data%ni(i)*z*z
            ld(i,j)=-0.5_dp*data%ni(i)*log(2.0_dp*acos(-1.0_dp)*b2)-0.5_dp*resid/b2
            dt(i,j)=data%ni(i)*z/b2; db(i,j,1)=-data%ni(i)/beta(1)+resid/(beta(1)**3)
         end do; end do
      case(NSP_MLOGIT)
         ng=nobs_data(data)
         do j=1,k
            do r=1,size(data%group)
               g=data%group(r); eta=pt(j)+dot_product(data%xmat(r,:),beta)
               if(eta>=0.0_dp) then; p=1.0_dp/(1.0_dp+exp(-min(eta,700.0_dp)))
               else; p=exp(max(eta,-700.0_dp))/(1.0_dp+exp(max(eta,-700.0_dp))); end if
               p=min(max(p,1.0e-10_dp),1.0_dp-1.0e-10_dp)
               ld(g,j)=ld(g,j)+data%y(r)*log(p)+(data%trials(r)-data%y(r))*log(1.0_dp-p)
               resid=data%y(r)-data%trials(r)*p; dt(g,j)=dt(g,j)+resid
               do i=1,q; db(g,j,i)=db(g,j,i)+resid*data%xmat(r,i); end do
            end do
         end do
      end select
   end subroutine
end module nspmix_families
