module trawl_simulation
  use trawl_kinds, only : dp
  use trawl_types, only : trawl_spec,trawl_ok,trawl_invalid_argument
  use trawl_functions, only : eval_trawl
  use trawl_rng, only : runif_scalar,rpois_scalar,rlogarithmic_scalar
  use trawl_distributions, only : bivariate_lsdsim
  implicit none
  private
  public :: sim_univariate_trawl,sim_bivariate_trawl
contains

  subroutine sim_univariate_trawl(t,x,delta_t,burnin,marginal,spec,v,m,theta,status)
    real(dp),intent(in)::t
    integer,allocatable,intent(out)::x(:)
    real(dp),intent(in),optional::delta_t,burnin,v,m,theta
    character(len=*),intent(in),optional::marginal
    type(trawl_spec),intent(in),optional::spec
    integer,intent(out),optional::status
    real(dp)::dt,bi,vv,mm,th,total,now,tdiff
    character(len=16)::marg
    type(trawl_spec)::sp
    integer::nj,nall,ncut,nout,i,k,s
    real(dp),allocatable::jt(:),jh(:)
    integer,allocatable::mark(:),allx(:)
    dt=1.0_dp;if(present(delta_t))dt=delta_t
    bi=10.0_dp;if(present(burnin))bi=burnin
    vv=0.0_dp;if(present(v))vv=v
    mm=0.0_dp;if(present(m))mm=m
    th=0.0_dp;if(present(theta))th=theta
    marg='Poi';if(present(marginal))marg=trim(marginal)
    sp%kind='Exp';if(present(spec))sp=spec
    s=trawl_ok
    if(t<=0.0_dp .or. dt<=0.0_dp .or. bi<0.0_dp) s=trawl_invalid_argument
    if(trim(marg)=='NegBin') then
      if(mm<0.0_dp .or. th<=0.0_dp .or. th>=1.0_dp) s=trawl_invalid_argument
      vv=mm*abs(log(1.0_dp-th))
    else if(trim(marg)/='Poi') then
      s=trawl_invalid_argument
    end if
    if(vv<0.0_dp) s=trawl_invalid_argument
    if(s/=trawl_ok) then
      allocate(x(0));if(present(status))status=s;return
    end if
    total=bi+t;nall=floor(total/dt);ncut=floor(bi/dt);nout=max(0,nall-ncut)
    nj=rpois_scalar(vv*total)
    allocate(jt(nj),jh(nj),mark(nj),allx(nall));allx=0
    do i=1,nj
      jt(i)=runif_scalar()*total;jh(i)=runif_scalar()
      if(trim(marg)=='NegBin') then; mark(i)=rlogarithmic_scalar(th); else; mark(i)=1; end if
    end do
    do k=1,nall
      now=real(k,dp)*dt
      do i=1,nj
        if(jt(i)<=now) then
          tdiff=now-jt(i)
          if(jh(i)<=eval_trawl(-tdiff,sp)) allx(k)=allx(k)+mark(i)
        end if
      end do
    end do
    allocate(x(nout));if(nout>0)x=allx(ncut+1:nall)
    if(present(status))status=s
  end subroutine

  subroutine sim_bivariate_trawl(t,x,delta_t,burnin,marginal,dependencetype,spec1,spec2, &
      v1,v2,v12,kappa1,kappa2,kappa12,a1,a2,status)
    real(dp),intent(in)::t
    integer,allocatable,intent(out)::x(:,:)
    real(dp),intent(in),optional::delta_t,burnin,v1,v2,v12,kappa1,kappa2,kappa12,a1,a2
    character(len=*),intent(in),optional::marginal,dependencetype
    type(trawl_spec),intent(in),optional::spec1,spec2
    integer,intent(out),optional::status
    real(dp)::dt,bi,total,vv1,vv2,vv12,k1,k2,k12,aa1,aa2,p1,p2,now,tdiff
    character(len=16)::marg,dep
    type(trawl_spec)::s1,s2
    integer::nall,ncut,nout,n12,n1,n2,i,k,s
    real(dp),allocatable::jt12(:),jh12(:),jt1(:),jh1(:),jt2(:),jh2(:)
    integer,allocatable::c12(:,:),c1(:),c2(:),aout(:,:),joint_marks(:,:)
    dt=1.0_dp;if(present(delta_t))dt=delta_t
    bi=10.0_dp;if(present(burnin))bi=burnin
    vv1=0.0_dp;if(present(v1))vv1=v1
    vv2=0.0_dp;if(present(v2))vv2=v2
    vv12=0.0_dp;if(present(v12))vv12=v12
    k1=0.0_dp;if(present(kappa1))k1=kappa1
    k2=0.0_dp;if(present(kappa2))k2=kappa2
    k12=0.0_dp;if(present(kappa12))k12=kappa12
    aa1=0.0_dp;if(present(a1))aa1=a1
    aa2=0.0_dp;if(present(a2))aa2=a2
    marg='Poi';if(present(marginal))marg=trim(marginal)
    dep='fullydep';if(present(dependencetype))dep=trim(dependencetype)
    s1%kind='Exp';if(present(spec1))s1=spec1
    s2%kind='Exp';if(present(spec2))s2=spec2
    s=trawl_ok
    if(t<=0.0_dp .or. dt<=0.0_dp .or. bi<0.0_dp) s=trawl_invalid_argument
    if(trim(marg)=='NegBin') then
      if(aa1<0.0_dp .or. aa2<0.0_dp .or. min(k1,k2,k12)<0.0_dp) s=trawl_invalid_argument
      vv12=k12*log(1.0_dp+aa1+aa2)
      if(trim(dep)=='dep') then
        vv1=k1*log(1.0_dp+aa1);vv2=k2*log(1.0_dp+aa2)
      end if
    else if(trim(marg)/='Poi') then
      s=trawl_invalid_argument
    end if
    if(trim(dep)/='fullydep' .and. trim(dep)/='dep') s=trawl_invalid_argument
    if(s/=trawl_ok) then
      allocate(x(0,2));if(present(status))status=s;return
    end if
    total=bi+t;nall=floor(total/dt);ncut=floor(bi/dt);nout=max(0,nall-ncut)
    allocate(aout(nall,2));aout=0

    n12=rpois_scalar(vv12*total)
    if(trim(dep)=='dep' .and. n12==0) n12=rpois_scalar(vv12*total)
    allocate(jt12(n12),jh12(n12),c12(n12,2))
    do i=1,n12;jt12(i)=runif_scalar()*total;jh12(i)=runif_scalar();end do
    if(trim(marg)=='NegBin' .and. n12>0) then
      p1=aa1/(aa1+aa2+1.0_dp);p2=aa2/(aa1+aa2+1.0_dp)
      call bivariate_lsdsim(n12,p1,p2,joint_marks,s);c12=joint_marks
    else
      c12=1
    end if

    if(trim(dep)=='fullydep') then
      do k=1,nall
        now=real(k,dp)*dt
        do i=1,n12
          if(jt12(i)<=now) then
            tdiff=now-jt12(i)
            if(jh12(i)<=eval_trawl(-tdiff,s1)) aout(k,1)=aout(k,1)+c12(i,1)
            if(jh12(i)<=eval_trawl(-tdiff,s2)) aout(k,2)=aout(k,2)+c12(i,2)
          end if
        end do
      end do
    else
      if(vv1==0.0_dp)then;n1=0;else;n1=rpois_scalar(vv1*total);end if
      if(vv2==0.0_dp)then;n2=0;else;n2=rpois_scalar(vv2*total);end if
      allocate(jt1(n1),jh1(n1),c1(n1),jt2(n2),jh2(n2),c2(n2))
      do i=1,n1
        jt1(i)=runif_scalar()*total;jh1(i)=runif_scalar()
        if(trim(marg)=='NegBin')then;c1(i)=rlogarithmic_scalar(aa1/(1.0_dp+aa1));else;c1(i)=1;end if
      end do
      do i=1,n2
        jt2(i)=runif_scalar()*total;jh2(i)=runif_scalar()
        if(trim(marg)=='NegBin')then;c2(i)=rlogarithmic_scalar(aa2/(1.0_dp+aa2));else;c2(i)=1;end if
      end do
      do k=1,nall
        now=real(k,dp)*dt
        do i=1,n12
          if(jt12(i)<=now) then
            tdiff=now-jt12(i)
            if(jh12(i)<=eval_trawl(-tdiff,s1))aout(k,1)=aout(k,1)+c12(i,1)
            if(jh12(i)<=eval_trawl(-tdiff,s2))aout(k,2)=aout(k,2)+c12(i,2)
          end if
        end do
        do i=1,n1
          if(jt1(i)<=now) then;tdiff=now-jt1(i);if(jh1(i)<=eval_trawl(-tdiff,s1))aout(k,1)=aout(k,1)+c1(i);end if
        end do
        do i=1,n2
          if(jt2(i)<=now) then;tdiff=now-jt2(i);if(jh2(i)<=eval_trawl(-tdiff,s2))aout(k,2)=aout(k,2)+c2(i);end if
        end do
      end do
    end if
    allocate(x(nout,2));if(nout>0)x=aout(ncut+1:nall,:)
    if(present(status))status=s
  end subroutine
end module trawl_simulation
