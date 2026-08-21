! SPDX-License-Identifier: GPL-2.0-or-later
module mc2d_stats
  use mc2d_kinds, only : dp, nan_dp
  use mc2d_utils, only : mean_dp, sd_dp, quantile_dp, correlation_dp
  use mc2d_node, only : mcnode, mc, mc_type_0, mc_type_v, mc_type_u, mc_type_vu, mcdata, extractvar
  implicit none
  private
  public :: node_summary, node_summary_each, node_quantiles, node_quantiles_each, mcratio, tornado, tornadounc
  public :: mcapply_reduce, mcapply_elemental, running_convergence

  type,public :: summary_result
    real(dp),allocatable::value(:,:)
    character(len=32),allocatable::row_name(:),col_name(:)
  end type
  type,public :: tornado_result
    real(dp),allocatable::value(:,:,:)
    character(len=64),allocatable::input_name(:),output_name(:),stat_name(:)
  end type

  abstract interface
    real(dp) function reduce_proc(x)
      import dp
      real(dp),intent(in)::x(:)
    end function reduce_proc
    real(dp) function scalar_proc(x)
      import dp
      real(dp),intent(in)::x
    end function scalar_proc
  end interface
contains
  function node_summary(x,probs,lim) result(s)
    type(mcnode),intent(in)::x
    real(dp),intent(in),optional::probs(:),lim(:)
    type(summary_result)::s
    real(dp),allocatable::pp(:),ll(:),tmp(:,:)
    integer::i,j,k,np,nl
    pp=[0.0_dp,0.025_dp,0.25_dp,0.5_dp,0.75_dp,0.975_dp,1.0_dp];if(present(probs))pp=probs
    ll=[0.025_dp,0.975_dp];if(present(lim))ll=lim;np=size(pp);nl=size(ll)
    select case(x%node_type)
    case(mc_type_0)
      allocate(s%value(1,1),s%row_name(1),s%col_name(1));s%value(1,1)=x%value(1,1,1);s%row_name='NoVar';s%col_name='NoUnc'
    case(mc_type_v)
      allocate(s%value(1,4+np),s%row_name(1),s%col_name(4+np));s%row_name='NoUnc'
      s%value(1,1)=mean_dp(x%value(:,1,1));s%value(1,2)=sd_dp(x%value(:,1,1))
      do i=1,np;s%value(1,2+i)=quantile_dp(x%value(:,1,1),pp(i));end do
      s%value(1,3+np)=real(x%nsv(),dp);s%value(1,4+np)=0
      s%col_name(1:2)=['mean                            ','sd                              ']
      do i=1,np;write(s%col_name(2+i),'(f6.2,"%")')100*pp(i);end do;s%col_name(3+np)='nsv';s%col_name(4+np)="Na's"
    case(mc_type_u)
      allocate(s%value(2+nl,1),s%row_name(2+nl),s%col_name(1));s%col_name='NoVar'
      s%value(1,1)=quantile_dp(x%value(1,:,1),0.5_dp);s%value(2,1)=mean_dp(x%value(1,:,1))
      do i=1,nl;s%value(2+i,1)=quantile_dp(x%value(1,:,1),ll(i));end do
      s%row_name(1)='median';s%row_name(2)='mean';do i=1,nl;write(s%row_name(2+i),'(f6.2,"%")')100*ll(i);end do
    case(mc_type_vu)
      allocate(tmp(x%nsu(),4+np))
      do j=1,x%nsu()
        tmp(j,1)=mean_dp(x%value(:,j,1));tmp(j,2)=sd_dp(x%value(:,j,1))
        do i=1,np;tmp(j,2+i)=quantile_dp(x%value(:,j,1),pp(i));end do
        tmp(j,3+np)=real(x%nsv(),dp);tmp(j,4+np)=0
      end do
      allocate(s%value(2+nl,4+np),s%row_name(2+nl),s%col_name(4+np))
      do k=1,4+np
        s%value(1,k)=quantile_dp(tmp(:,k),0.5_dp);s%value(2,k)=mean_dp(tmp(:,k))
        do i=1,nl;s%value(2+i,k)=quantile_dp(tmp(:,k),ll(i));end do
      end do
      s%row_name(1)='median';s%row_name(2)='mean';do i=1,nl;write(s%row_name(2+i),'(f6.2,"%")')100*ll(i);end do
      s%col_name(1:2)=['mean                            ','sd                              ']
      do i=1,np;write(s%col_name(2+i),'(f6.2,"%")')100*pp(i);end do;s%col_name(3+np)='nsv';s%col_name(4+np)="Na's"
    end select
  end function node_summary

  function node_summary_each(x,probs,lim) result(s)
    type(mcnode),intent(in)::x
    real(dp),intent(in),optional::probs(:),lim(:)
    type(summary_result),allocatable::s(:)
    type(mcnode)::one
    integer::k
    allocate(s(x%nvariates()))
    do k=1,x%nvariates()
      one=extractvar(x,[k])
      if(present(probs).and.present(lim))then
        s(k)=node_summary(one,probs=probs,lim=lim)
      else if(present(probs))then
        s(k)=node_summary(one,probs=probs)
      else if(present(lim))then
        s(k)=node_summary(one,lim=lim)
      else
        s(k)=node_summary(one)
      end if
    end do
  end function node_summary_each

  function node_quantiles(x,probs,lim) result(s)
    type(mcnode),intent(in)::x;real(dp),intent(in)::probs(:);real(dp),intent(in),optional::lim(:);type(summary_result)::s
    real(dp),allocatable::ll(:),prem(:,:);integer::i,j
    ll=[0.025_dp,0.975_dp];if(present(lim))ll=lim
    select case(x%node_type)
    case(mc_type_0)
      allocate(s%value(1,1));s%value=x%value(1,1,1)
    case(mc_type_v)
      allocate(s%value(1,size(probs)));do i=1,size(probs);s%value(1,i)=quantile_dp(x%value(:,1,1),probs(i));end do
    case(mc_type_u)
      allocate(s%value(2+size(ll),1));s%value(1,1)=quantile_dp(x%value(1,:,1),0.5_dp);s%value(2,1)=mean_dp(x%value(1,:,1))
      do i=1,size(ll);s%value(2+i,1)=quantile_dp(x%value(1,:,1),ll(i));end do
    case(mc_type_vu)
      allocate(prem(size(probs),x%nsu()))
      do j=1,x%nsu()
      do i=1,size(probs)
      prem(i,j)=quantile_dp(x%value(:,j,1),probs(i))
      end do
      end do
      allocate(s%value(2+size(ll),size(probs)))
      do i=1,size(probs)
      s%value(1,i)=quantile_dp(prem(i,:),0.5_dp)
      s%value(2,i)=mean_dp(prem(i,:))
      do j=1,size(ll)
      s%value(2+j,i)=quantile_dp(prem(i,:),ll(j))
      end do
      end do
    end select
  end function node_quantiles

  function node_quantiles_each(x,probs,lim) result(s)
    type(mcnode),intent(in)::x
    real(dp),intent(in)::probs(:)
    real(dp),intent(in),optional::lim(:)
    type(summary_result),allocatable::s(:)
    type(mcnode)::one
    integer::k
    allocate(s(x%nvariates()))
    do k=1,x%nvariates()
      one=extractvar(x,[k])
      if(present(lim))then
        s(k)=node_quantiles(one,probs,lim)
      else
        s(k)=node_quantiles(one,probs)
      end if
    end do
  end function node_quantiles_each

  function mcratio(x,pcentral,pvar,punc) result(r)
    type(mcnode),intent(in)::x;real(dp),intent(in),optional::pcentral,pvar,punc;real(dp)::r(7)
    real(dp)::pc,pv,pu,a,b,c,d;real(dp),allocatable::qv(:),qv2(:);integer::j
    pc=.5_dp;if(present(pcentral))pc=pcentral;pv=.975_dp;if(present(pvar))pv=pvar;pu=.975_dp;if(present(punc))pu=punc
    select case(x%node_type)
    case(mc_type_vu)
      allocate(qv(x%nsu()),qv2(x%nsu()))
      do j=1,x%nsu()
      qv(j)=quantile_dp(x%value(:,j,1),pc)
      qv2(j)=quantile_dp(x%value(:,j,1),pv)
      end do
      a=quantile_dp(qv,pc);b=quantile_dp(qv2,pc);c=quantile_dp(qv,pu);d=quantile_dp(qv2,pu)
    case(mc_type_v)
      a=quantile_dp(x%value(:,1,1),pc);b=quantile_dp(x%value(:,1,1),pv);c=a;d=b
    case(mc_type_u)
      a=quantile_dp(x%value(1,:,1),pc);b=a;c=quantile_dp(x%value(1,:,1),pu);d=c
    case default;a=x%value(1,1,1);b=a;c=a;d=a
    end select
    r=[a,b,c,d,b/a,c/a,d/a]
  end function mcratio

  logical function node_output_none(x) result(is_none)
    type(mcnode),intent(in)::x
    is_none=.false.
    if(allocated(x%outm))is_none=trim(x%outm)=='none'
  end function node_output_none

  subroutine require_supported_outm(x,where)
    type(mcnode),intent(in)::x
    character(len=*),intent(in)::where
    if(x%nvariates()<=1)return
    if(.not.allocated(x%outm))return
    if(trim(x%outm)=='each')return
    error stop trim(where)//': multivariate custom outm reducers are not supported; use outm="each"'
  end subroutine require_supported_outm

  function variate_name(base,k,nv) result(name)
    character(len=*),intent(in)::base
    integer,intent(in)::k,nv
    character(len=64)::name
    if(nv==1)then
      name=trim(base)
    else
      write(name,'(a,".",i0)')trim(base),k
    end if
  end function variate_name

  function variability_stat_name(k,q) result(name)
    integer,intent(in)::k
    real(dp),intent(in)::q(:)
    character(len=32)::name
    if(k==1)then
      name='mean'
    else if(k==2)then
      name='sd'
    else if(q(k-2)==0.0_dp)then
      name='Min'
    else if(q(k-2)==1.0_dp)then
      name='Max'
    else
      write(name,'(f6.2,"%")')100.0_dp*q(k-2)
      name=adjustl(name)
    end if
  end function variability_stat_name

  subroutine variability_stats(x,j,k,q,stat)
    type(mcnode),intent(in)::x
    integer,intent(in)::j,k
    real(dp),intent(in)::q(:)
    real(dp),intent(out)::stat(:)
    integer::iq,jj
    jj=min(j,x%nsu())
    if(size(stat)/=2+size(q))error stop 'variability_stats: invalid output size'
    stat(1)=mean_dp(x%value(:,jj,k))
    stat(2)=sd_dp(x%value(:,jj,k))
    do iq=1,size(q)
      stat(2+iq)=quantile_dp(x%value(:,jj,k),q(iq))
    end do
  end subroutine variability_stats

  function tornado(model,output,method,lim) result(tr)
    type(mc),intent(in)::model
    integer,intent(in)::output
    character(len=*),intent(in),optional::method
    real(dp),intent(in),optional::lim(2)
    type(tornado_result)::tr
    character(len=16)::meth
    real(dp)::ll(2)
    integer::i,j,ji,ki,ko,nin,nout,nco,icol,nrow
    real(dp),allocatable::corr(:)
    logical::eligible

    meth='spearman'
    if(present(method))meth=method
    ll=[0.025_dp,0.975_dp]
    if(present(lim))ll=lim
    if(output<1.or.output>model%size())error stop 'tornado: invalid output'
    if(model%node(output)%node_type/=mc_type_v .and. model%node(output)%node_type/=mc_type_vu) &
      error stop 'tornado: output must be V or VU'
    if(node_output_none(model%node(output)))error stop 'tornado: output has outm="none"'
    call require_supported_outm(model%node(output),'tornado')

    nin=0
    do i=1,model%size()
      if(i==output .or. node_output_none(model%node(i)))cycle
      eligible=model%node(i)%node_type==mc_type_v
      if(model%node(output)%node_type==mc_type_vu) &
        eligible=eligible.or.model%node(i)%node_type==mc_type_vu
      if(.not.eligible)cycle
      call require_supported_outm(model%node(i),'tornado')
      nin=nin+model%node(i)%nvariates()
    end do
    if(nin==0)error stop 'tornado: no valid input nodes'

    nout=model%node(output)%nvariates()
    nco=model%node(output)%nsu()
    nrow=merge(1,4,model%node(output)%node_type==mc_type_v)
    allocate(tr%value(nrow,nin,nout),tr%stat_name(nrow))
    allocate(tr%input_name(nin),tr%output_name(nout),corr(nco))
    if(nrow==1)then
      tr%stat_name(1)='corr'
    else
      tr%stat_name=['median                          ', &
                    'mean                            ', &
                    'lower                           ', &
                    'upper                           ']
    end if
    do ko=1,nout
      tr%output_name(ko)=variate_name(model%name(output),ko,nout)
    end do

    icol=0
    do i=1,model%size()
      if(i==output .or. node_output_none(model%node(i)))cycle
      eligible=model%node(i)%node_type==mc_type_v
      if(model%node(output)%node_type==mc_type_vu) &
        eligible=eligible.or.model%node(i)%node_type==mc_type_vu
      if(.not.eligible)cycle
      do ki=1,model%node(i)%nvariates()
        icol=icol+1
        tr%input_name(icol)=variate_name(model%name(i),ki,model%node(i)%nvariates())
        do ko=1,nout
          do j=1,nco
            ji=min(j,size(model%node(i)%value,2))
            corr(j)=correlation_dp(model%node(output)%value(:,j,ko), &
              model%node(i)%value(:,ji,ki),meth)
          end do
          if(nrow==1)then
            tr%value(1,icol,ko)=corr(1)
          else
            tr%value(1,icol,ko)=quantile_dp(corr,0.5_dp)
            tr%value(2,icol,ko)=mean_dp(corr)
            tr%value(3,icol,ko)=quantile_dp(corr,ll(1))
            tr%value(4,icol,ko)=quantile_dp(corr,ll(2))
          end if
        end do
      end do
    end do
  end function tornado

  function tornadounc(model,output,method,quant) result(tr)
    type(mc),intent(in)::model
    integer,intent(in)::output
    character(len=*),intent(in),optional::method
    real(dp),intent(in),optional::quant(:)
    type(tornado_result)::tr
    character(len=16)::meth
    real(dp),allocatable::q(:),outstats(:,:),instats(:,:),tmp(:)
    integer::i,j,ji,ki,ko,si,so,nin,nout,nco,nstat,icol,noutstat
    logical::eligible
    character(len=64)::vname

    meth='spearman'
    if(present(method))meth=method
    q=[0.5_dp,0.75_dp,0.975_dp]
    if(present(quant))q=quant
    if(output<1.or.output>model%size())error stop 'tornadounc: invalid output'
    if(model%node(output)%node_type/=mc_type_u .and. model%node(output)%node_type/=mc_type_vu) &
      error stop 'tornadounc: output must be U or VU'
    if(node_output_none(model%node(output)))error stop 'tornadounc: output has outm="none"'
    call require_supported_outm(model%node(output),'tornadounc')

    nstat=2+size(q)
    noutstat=merge(1,nstat,model%node(output)%node_type==mc_type_u)
    nin=0
    do i=1,model%size()
      if(i==output .or. node_output_none(model%node(i)))cycle
      eligible=model%node(i)%node_type==mc_type_u
      if(model%node(output)%node_type==mc_type_vu) &
        eligible=eligible.or.model%node(i)%node_type==mc_type_vu
      if(.not.eligible)cycle
      call require_supported_outm(model%node(i),'tornadounc')
      if(model%node(i)%node_type==mc_type_u)then
        nin=nin+model%node(i)%nvariates()
      else
        nin=nin+model%node(i)%nvariates()*nstat
      end if
    end do
    if(nin==0)error stop 'tornadounc: no valid input nodes'

    nout=model%node(output)%nvariates()
    nco=model%node(output)%nsu()
    allocate(tr%value(noutstat,nin,nout),tr%stat_name(noutstat))
    allocate(tr%input_name(nin),tr%output_name(nout))
    if(noutstat==1)then
      tr%stat_name(1)='corr'
    else
      do so=1,nstat
        tr%stat_name(so)=variability_stat_name(so,q)
      end do
    end if
    do ko=1,nout
      tr%output_name(ko)=variate_name(model%name(output),ko,nout)
    end do

    allocate(outstats(nco,noutstat))
    icol=0
    do ko=1,nout
      if(model%node(output)%node_type==mc_type_u)then
        do j=1,nco
          outstats(j,1)=model%node(output)%value(1,j,ko)
        end do
      else
        allocate(tmp(nstat))
        do j=1,nco
          call variability_stats(model%node(output),j,ko,q,tmp)
          outstats(j,:)=tmp
        end do
        deallocate(tmp)
      end if

      icol=0
      do i=1,model%size()
        if(i==output .or. node_output_none(model%node(i)))cycle
        eligible=model%node(i)%node_type==mc_type_u
        if(model%node(output)%node_type==mc_type_vu) &
          eligible=eligible.or.model%node(i)%node_type==mc_type_vu
        if(.not.eligible)cycle
        do ki=1,model%node(i)%nvariates()
          vname=variate_name(model%name(i),ki,model%node(i)%nvariates())
          if(model%node(i)%node_type==mc_type_u)then
            allocate(instats(nco,1))
            do j=1,nco
              ji=min(j,size(model%node(i)%value,2))
              instats(j,1)=model%node(i)%value(1,ji,ki)
            end do
            icol=icol+1
            if(ko==1)tr%input_name(icol)=vname
            do so=1,noutstat
              tr%value(so,icol,ko)=correlation_dp(outstats(:,so),instats(:,1),meth)
            end do
            deallocate(instats)
          else
            allocate(instats(nco,nstat),tmp(nstat))
            do j=1,nco
              call variability_stats(model%node(i),j,ki,q,tmp)
              instats(j,:)=tmp
            end do
            do si=1,nstat
              icol=icol+1
              if(ko==1)tr%input_name(icol)=trim(variability_stat_name(si,q))//' '//trim(vname)
              do so=1,noutstat
                tr%value(so,icol,ko)=correlation_dp(outstats(:,so),instats(:,si),meth)
              end do
            end do
            deallocate(instats,tmp)
          end if
        end do
      end do
    end do
  end function tornadounc

  function mcapply_reduce(x,margin,fun) result(r)
    type(mcnode),intent(in)::x;character(len=*),intent(in)::margin;procedure(reduce_proc)::fun;type(mcnode)::r
    integer::i,j,k
    select case(trim(margin))
    case('all');r=mcdata(fun(reshape(x%value,[size(x%value)])),type='0')
    case('var')
    allocate(r%value(1,x%nsu(),x%nvariates()))
    do k=1,x%nvariates()
    do j=1,x%nsu()
    r%value(1,j,k)=fun(x%value(:,j,k))
    end do
    end do
    r%node_type=merge(mc_type_0,mc_type_u,x%nsu()==1)
    r%outm=x%outm
    case('unc')
    allocate(r%value(x%nsv(),1,x%nvariates()))
    do k=1,x%nvariates()
    do i=1,x%nsv()
    r%value(i,1,k)=fun(x%value(i,:,k))
    end do
    end do
    r%node_type=merge(mc_type_0,mc_type_v,x%nsv()==1)
    r%outm=x%outm
    case('variates')
    allocate(r%value(x%nsv(),x%nsu(),1))
    do j=1,x%nsu()
    do i=1,x%nsv()
    r%value(i,j,1)=fun(x%value(i,j,:))
    end do
    end do
    r%node_type=x%node_type
    r%outm=x%outm
    case default;error stop 'mcapply_reduce: invalid margin'
    end select
  end function mcapply_reduce

  function mcapply_elemental(x,fun) result(r)
    type(mcnode),intent(in)::x
    procedure(scalar_proc)::fun
    type(mcnode)::r
    integer::i,j,k
    r=x
    do k=1,size(x%value,3)
      do j=1,size(x%value,2)
        do i=1,size(x%value,1)
          r%value(i,j,k)=fun(x%value(i,j,k))
        end do
      end do
    end do
  end function mcapply_elemental

  subroutine running_convergence(x,probs,out)
    real(dp),intent(in)::x(:),probs(:);real(dp),intent(out)::out(1+size(probs),size(x));integer::i,j
    do i=1,size(x);out(1,i)=mean_dp(x(:i));do j=1,size(probs);out(1+j,i)=quantile_dp(x(:i),probs(j));end do;end do
  end subroutine running_convergence
end module mc2d_stats
