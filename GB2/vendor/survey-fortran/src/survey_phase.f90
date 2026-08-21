! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_phase
  use survey_kinds, only : dp
  use survey_types, only : svystat_t, phase_variance_t
  use survey_pps, only : ht_variance
  use survey_linalg, only : sym_pinv
  implicit none
  private
  public :: dcheck_strat, dcheck_multi, dcheck_multi_subset, combine_dcheck, &
            twophase2_variance, multiphase_variance, multiphase_total, multiphase_mean, project_phase_calibration
contains

  subroutine dcheck_strat(strata,prob,dcheck)
    integer,intent(in) :: strata(:)
    real(dp),intent(in) :: prob(:)
    real(dp),intent(out) :: dcheck(:,:)
    integer :: n,i,j,nh
    n=size(strata);if(size(prob)/=n.or.any(shape(dcheck)/=[n,n]))error stop 'dcheck_strat: shape mismatch'
    if(any(prob<0.0_dp).or.any(prob>1.0_dp))error stop 'dcheck_strat: probabilities outside [0,1]'
    dcheck=0.0_dp
    do i=1,n
      nh=count(strata==strata(i))
      if(nh>1)then
        do j=1,n
          if(i/=j.and.strata(j)==strata(i))dcheck(i,j)=-(1.0_dp-prob(i))/real(nh-1,dp)
        end do
      end if
      dcheck(i,i)=1.0_dp-prob(i)
    end do
  end subroutine dcheck_strat

  subroutine combine_dcheck(d1,d2,dout)
    real(dp),intent(in) :: d1(:,:),d2(:,:)
    real(dp),intent(out) :: dout(:,:)
    if(any(shape(d1)/=shape(d2)).or.any(shape(dout)/=shape(d1)))error stop 'combine_dcheck: shape mismatch'
    dout=-d1*d2+d1+d2
  end subroutine combine_dcheck

  subroutine dcheck_multi(id,strata,prob,dcheck)
    integer,intent(in) :: id(:,:),strata(:,:)
    real(dp),intent(in) :: prob(:,:)
    real(dp),intent(out) :: dcheck(:,:)
    integer,allocatable :: uid(:),map(:),srep(:)
    real(dp),allocatable :: prep(:),dc_u(:,:),stage_dc(:,:),tmp(:,:)
    integer :: n,ns,s,i,k
    n=size(id,1);ns=size(id,2)
    if(any(shape(strata)/=shape(id)).or.any(shape(prob)/=shape(id)).or.any(shape(dcheck)/=[n,n]))error stop 'dcheck_multi: shape mismatch'
    allocate(stage_dc(n,n),tmp(n,n));dcheck=0.0_dp
    do s=1,ns
      call unique_map(id(:,s),uid,map)
      allocate(srep(size(uid)),prep(size(uid)),dc_u(size(uid),size(uid)))
      do k=1,size(uid)
        i=find_first_int(map,k);srep(k)=strata(i,s);prep(k)=prob(i,s)
      end do
      call dcheck_strat(srep,prep,dc_u)
      do i=1,n;do k=1,n;stage_dc(i,k)=dc_u(map(i),map(k));end do;end do
      call combine_dcheck(dcheck,stage_dc,tmp);dcheck=tmp
      deallocate(uid,map,srep,prep,dc_u)
    end do
  end subroutine dcheck_multi

  subroutine dcheck_multi_subset(id,strata,subset,prob,dcheck,with_replacement,phase1_iid)
    integer,intent(in) :: id(:,:),strata(:,:)
    logical,intent(in) :: subset(:)
    real(dp),intent(in) :: prob(:,:)
    real(dp),intent(out) :: dcheck(:,:)
    logical,intent(in),optional :: with_replacement,phase1_iid
    integer,allocatable :: selected_rows(:),uid(:),map_sel(:),srep(:),nh(:),all_uid(:),all_map(:)
    real(dp),allocatable :: prep(:),dc_u(:,:),stage_dc(:,:),tmp(:,:)
    integer :: n,ns,nsel,s,i,j,k,repidx
    logical :: wr,iid
    n=size(id,1);ns=size(id,2);nsel=count(subset);wr=.false.;iid=.false.;if(present(with_replacement))wr=with_replacement;if(present(phase1_iid))iid=phase1_iid
    if(size(subset)/=n.or.any(shape(strata)/=shape(id)).or.any(shape(prob)/=shape(id)).or.any(shape(dcheck)/=[nsel,nsel]))error stop 'dcheck_multi_subset: shape mismatch'
    dcheck=0.0_dp;if(iid)then;do i=1,nsel;dcheck(i,i)=1.0_dp;end do;return;end if
    allocate(selected_rows(nsel));selected_rows=pack([(i,i=1,n)],subset);allocate(stage_dc(nsel,nsel),tmp(nsel,nsel))
    do s=1,ns
      call unique_map(pack(id(:,s),subset),uid,map_sel);call unique_map(id(:,s),all_uid,all_map)
      allocate(srep(size(uid)),prep(size(uid)),nh(size(uid)),dc_u(size(uid),size(uid)));dc_u=0.0_dp
      do k=1,size(uid)
        repidx=selected_rows(find_first_int(map_sel,k));srep(k)=strata(repidx,s);prep(k)=prob(repidx,s)
        nh(k)=count_unique_in_stratum(id(:,s),strata(:,s),srep(k))
      end do
      do i=1,size(uid)
        dc_u(i,i)=1.0_dp-prep(i)
        if(.not.wr.and.nh(i)>1)then
          do j=1,size(uid)
            if(i/=j.and.srep(j)==srep(i))dc_u(i,j)=-(1.0_dp-prep(i))/real(nh(i)-1,dp)
          end do
        end if
      end do
      do i=1,nsel;do j=1,nsel;stage_dc(i,j)=dc_u(map_sel(i),map_sel(j));end do;end do
      call combine_dcheck(dcheck,stage_dc,tmp);dcheck=tmp
      deallocate(uid,map_sel,srep,prep,nh,dc_u,all_uid,all_map)
    end do
  end subroutine dcheck_multi_subset

  subroutine twophase2_variance(x,dcheck_full,dcheck_phase2,result)
    real(dp),intent(in) :: x(:,:),dcheck_full(:,:),dcheck_phase2(:,:)
    type(phase_variance_t),intent(out) :: result
    real(dp),allocatable :: vfull(:,:),v2(:,:)
    integer :: p
    if(size(dcheck_full,1)/=size(x,1).or.size(dcheck_full,2)/=size(x,1)) error stop 'twophase2_variance: full Dcheck shape'
    if(any(shape(dcheck_phase2)/=shape(dcheck_full))) error stop 'twophase2_variance: phase-2 Dcheck shape'
    p=size(x,2);allocate(vfull(p,p),v2(p,p),result%variance(p,p),result%phase(p,p,2))
    vfull=ht_variance(x,dcheck_full);v2=ht_variance(x,dcheck_phase2)
    result%phase(:,:,1)=vfull-v2
    result%phase(:,:,2)=v2
    result%variance=vfull
  end subroutine twophase2_variance

  subroutine multiphase_variance(x,phase_weight,dcheck,result,calibration_x,calibration_weight)
    real(dp),intent(in) :: x(:,:),phase_weight(:,:),dcheck(:,:,:)
    type(phase_variance_t),intent(out) :: result
    real(dp),intent(in),optional :: calibration_x(:,:,:),calibration_weight(:,:)
    real(dp),allocatable :: z(:,:),zp(:,:)
    integer :: n,p,h,np
    n=size(x,1);p=size(x,2);np=size(phase_weight,2)
    if(size(phase_weight,1)/=n) error stop 'multiphase_variance: phase_weight rows'
    if(size(dcheck,1)/=n.or.size(dcheck,2)/=n.or.size(dcheck,3)/=np) error stop 'multiphase_variance: Dcheck shape'
    if(present(calibration_x).neqv.present(calibration_weight)) error stop 'multiphase_variance: calibration arguments must be paired'
    if(present(calibration_x)) then
      if(size(calibration_x,1)/=n.or.size(calibration_x,3)/=np.or.any(shape(calibration_weight)/=[n,np])) &
        error stop 'multiphase_variance: calibration shape'
    end if
    allocate(result%variance(p,p),result%phase(p,p,np),z(n,p),zp(n,p));result%variance=0.0_dp
    do h=1,np
      z=x*spread(phase_weight(:,h),2,p)
      if(present(calibration_x)) then
        call project_phase_calibration(z,calibration_x(:,:,h),calibration_weight(:,h),zp)
      else
        zp=z
      end if
      result%phase(:,:,h)=ht_variance(zp,dcheck(:,:,h))
      result%variance=result%variance+result%phase(:,:,h)
    end do
  end subroutine multiphase_variance

  function multiphase_total(x,final_weight,phase_weight,dcheck,calibration_x,calibration_weight) result(ans)
    real(dp),intent(in) :: x(:,:),final_weight(:),phase_weight(:,:),dcheck(:,:,:)
    real(dp),intent(in),optional :: calibration_x(:,:,:),calibration_weight(:,:)
    type(svystat_t) :: ans
    type(phase_variance_t) :: pv
    integer :: j,p
    if(size(x,1)/=size(final_weight)) error stop 'multiphase_total: final_weight size'
    p=size(x,2);allocate(ans%estimate(p),ans%variance(p,p),ans%influence(size(x,1),p))
    do j=1,p;ans%estimate(j)=dot_product(final_weight,x(:,j));end do
    if(present(calibration_x)) then
      call multiphase_variance(x,phase_weight,dcheck,pv,calibration_x,calibration_weight)
    else
      call multiphase_variance(x,phase_weight,dcheck,pv)
    end if
    ans%variance=pv%variance
    ans%influence=x*spread(final_weight,2,p)
  end function multiphase_total

  function multiphase_mean(x,final_weight,phase_weight,dcheck,calibration_x,calibration_weight) result(ans)
    real(dp),intent(in) :: x(:,:),final_weight(:),phase_weight(:,:),dcheck(:,:,:)
    real(dp),intent(in),optional :: calibration_x(:,:,:),calibration_weight(:,:)
    type(svystat_t) :: ans
    type(phase_variance_t) :: pv
    real(dp),allocatable :: lin(:,:)
    real(dp) :: sw
    integer :: j,p,n
    n=size(x,1);p=size(x,2);if(size(final_weight)/=n) error stop 'multiphase_mean: final_weight size'
    sw=sum(final_weight);if(sw<=0.0_dp) error stop 'multiphase_mean: nonpositive total weight'
    allocate(ans%estimate(p),ans%variance(p,p),ans%influence(n,p),lin(n,p))
    do j=1,p;ans%estimate(j)=dot_product(final_weight,x(:,j))/sw;lin(:,j)=(x(:,j)-ans%estimate(j))/sw;end do
    if(present(calibration_x)) then
      call multiphase_variance(lin,phase_weight,dcheck,pv,calibration_x,calibration_weight)
    else
      call multiphase_variance(lin,phase_weight,dcheck,pv)
    end if
    ans%variance=pv%variance
    ans%influence=lin*spread(final_weight,2,p)
  end function multiphase_mean

  subroutine project_phase_calibration(x,z,w,resid)
    ! Equivalent to upstream project_ps(): qr.resid(calibration_qr_x, x/scale)*scale.
    ! z is the matrix used to construct the upstream QR (typically model_matrix*sqrt(phase_weight)); w is the post-calibration scale.
    real(dp),intent(in) :: x(:,:),z(:,:),w(:)
    real(dp),intent(out) :: resid(:,:)
    real(dp),allocatable :: zw(:,:),a(:,:),ainv(:,:),rhs(:,:),coef(:,:)
    integer :: n,k,p,i,j,l,rank,info
    n=size(x,1);p=size(x,2);k=size(z,2)
    if(size(z,1)/=n.or.size(w)/=n.or.any(shape(resid)/=[n,p])) error stop 'project_phase_calibration: shape mismatch'
    if(any(w<=0.0_dp)) error stop 'project_phase_calibration: weights must be positive'
    allocate(zw(n,k),a(k,k),ainv(k,k),rhs(k,p),coef(k,p));zw=z;a=0.0_dp;rhs=0.0_dp
    do i=1,n
      do j=1,k
        do l=1,k;a(j,l)=a(j,l)+zw(i,j)*zw(i,l);end do
        rhs(j,:)=rhs(j,:)+zw(i,j)*(x(i,:)/w(i))
      end do
    end do
    call sym_pinv(a,ainv,rank,info=info);coef=matmul(ainv,rhs)
    resid=(x/spread(w,2,p)-matmul(zw,coef))*spread(w,2,p)
  end subroutine project_phase_calibration

  subroutine unique_map(x,u,map)
    integer,intent(in) :: x(:)
    integer,allocatable,intent(out) :: u(:),map(:)
    integer,allocatable :: tmp(:)
    integer :: i,j,m
    allocate(tmp(size(x)),map(size(x)));m=0
    do i=1,size(x)
      j=0
      if(m>0)then
        do while(j<m)
          j=j+1;if(tmp(j)==x(i))exit
        end do
        if(tmp(j)/=x(i))j=0
      end if
      if(j==0)then;m=m+1;tmp(m)=x(i);j=m;end if
      map(i)=j
    end do
    allocate(u(m));if(m>0)u=tmp(1:m)
  end subroutine unique_map

  integer function find_first_int(x,v) result(k)
    integer,intent(in)::x(:),v;integer::i;k=0;do i=1,size(x);if(x(i)==v)then;k=i;return;end if;end do
  end function find_first_int

  integer function count_unique_in_stratum(id,strata,svalue) result(nu)
    integer,intent(in)::id(:),strata(:),svalue
    integer,allocatable::u(:),m(:)
    call unique_map(pack(id,strata==svalue),u,m);nu=size(u)
  end function count_unique_in_stratum
end module survey_phase
