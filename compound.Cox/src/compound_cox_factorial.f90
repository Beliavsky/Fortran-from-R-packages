! SPDX-License-Identifier: GPL-2.0-only
module compound_cox_factorial
  use compound_cox_kinds, only : dp
  use compound_cox_types, only : factorial_result, cg_result
  use compound_cox_cg, only : cg_clayton, cg_frank, cg_gumbel
  use compound_cox_math, only : covariance_matrix, pseudoinverse_sym, symmetric_eigenvalues, &
       normal_cdf, normal_quantile, chisq_cdf, chisq_quantile, randn
  implicit none
  private
  public :: surv_factorial
contains
  subroutine cg_dispatch(time,status,alpha,copula,res)
    real(dp),intent(in)::time(:),alpha
    integer,intent(in)::status(:)
    character(len=*),intent(in)::copula
    type(cg_result),intent(out)::res
    select case(copula(1:1))
    case('c','C')
    call cg_clayton(time,status,alpha,res)
    case('f','F')
    call cg_frank(time,status,alpha,res)
    case default
    call cg_gumbel(time,status,alpha,res)
    end select
  end subroutine cg_dispatch

  subroutine p_func(time,status,group,d,alpha,copula,tupper,pvec)
    real(dp),intent(in)::time(:),alpha,tupper
    integer,intent(in)::status(:),group(:),d
    character(len=*),intent(in)::copula
    real(dp),intent(out)::pvec(d)
    real(dp),allocatable::w(:,:),ti(:),tl(:),sip(:),sin(:),slp(:),sln(:),dsl(:)
    integer,allocatable::di(:),dl(:)
    type(cg_result)::ci,cl
    integer::i,l,ni,nl,j,countle,ku
    real(dp)::siu,slu
    allocate(w(d,d))
    w=0
    do i=1,d
      call get_group(time,status,group,i,ti,di)
      ni=size(ti)
      call cg_dispatch(ti,di,alpha,copula,ci)
      allocate(sip(ni),sin(ni))
      sip=ci%surv
      sip(ni)=0
      sin(1)=1
      if(ni>1)sin(2:)=ci%surv(:ni-1)
      ku=count(ci%time<=tupper)
      if(ku>0)then
      siu=sip(ku)
      else
      siu=1
      end if
      do l=1,d
        call get_group(time,status,group,l,tl,dl)
        nl=size(tl)
        call cg_dispatch(tl,dl,alpha,copula,cl)
        allocate(slp(nl),sln(nl),dsl(nl))
        slp=cl%surv
        slp(nl)=0
        sln(1)=1
        if(nl>1)sln(2:)=cl%surv(:nl-1)
        dsl=sln-slp
        do j=1,nl
        if(cl%time(j)>tupper)dsl(j)=0
        end do
        ku=count(cl%time<=tupper)
        if(ku>0)then
        slu=slp(ku)
        else
        slu=1
        end if
        do j=1,nl
        countle=count(ci%time<=cl%time(j))
        countle=max(countle,1)
        w(i,l)=w(i,l)+(sip(countle)+sin(countle))*dsl(j)/2
        end do
        w(i,l)=w(i,l)+siu*slu/2
        deallocate(tl,dl,slp,sln,dsl)
      end do
      deallocate(ti,di,sip,sin)
    end do
    do i=1,d
    pvec(i)=sum(w(i,:))/real(d,dp)
    end do
  end subroutine p_func

  subroutine get_group(time,status,group,g,t,dv)
    real(dp),intent(in)::time(:)
    integer,intent(in)::status(:),group(:),g
    real(dp),allocatable,intent(out)::t(:)
    integer,allocatable,intent(out)::dv(:)
    integer::m,i,k
    m=count(group==g)
    allocate(t(m),dv(m))
    k=0
    do i=1,size(time)
    if(group(i)==g)then
    k=k+1
    t(k)=time(i)
    dv(k)=status(i)
    end if
    end do
  end subroutine get_group

  subroutine surv_factorial(time,status,group,alpha,res,copula,nsim,t_upper,contrast)
    real(dp),intent(in)::time(:),alpha
    integer,intent(in)::status(:),group(:)
    type(factorial_result),intent(out)::res
    character(len=*),intent(in),optional::copula
    integer,intent(in),optional::nsim
    real(dp),intent(in),optional::t_upper
    real(dp),intent(in),optional::contrast(:,:)
    integer::d,n,r,i,ns
    real(dp)::tu,denom
    character(len=16)::cop
    real(dp),allocatable::p(:),ptmp(:),del(:,:),covd(:,:),jvar(:,:),c(:,:),cc(:,:),pinv(:,:),tmat(:,:),tv(:,:),eval(:),qsim(:),tt(:)
    integer,allocatable::st(:),gr(:)
    d=maxval(group)
    n=size(time)
    cop='clayton'
    if(present(copula))cop=copula
    ns=1000
    if(present(nsim))ns=nsim
    if(present(t_upper))then
    tu=t_upper
    else
    tu=huge(1.0_dp)
    do i=1,d
    tu=min(tu,maxval(time,mask=group==i))
    end do
    end if
    allocate(p(d),ptmp(d))
    call p_func(time,status,group,d,alpha,cop,tu,p)
    allocate(del(n,d))
    do i=1,n
    allocate(tt(n-1),st(n-1),gr(n-1))
    tt=[time(:i-1),time(i+1:)]
    st=[status(:i-1),status(i+1:)]
    gr=[group(:i-1),group(i+1:)]
    call p_func(tt,st,gr,d,alpha,cop,tu,ptmp)
    del(i,:)=ptmp
    deallocate(tt,st,gr)
    end do
    call covariance_matrix(del,covd)
    allocate(jvar(d,d))
    jvar=real((n-1)*(n-1),dp)/real(n,dp)*covd
    allocate(res%estimate(d),res%se(d),res%lower(d),res%upper(d),res%p(d),res%variance(d,d))
    res%estimate=p
    res%variance=jvar
    res%se=sqrt(max(diagonal(jvar),0.0_dp))
    do i=1,d
    if(res%se(i)>0)then
    res%p(i)=1-normal_cdf((p(i)-0.5_dp)/res%se(i))
    else
    res%p(i)=1
    end if
    res%lower(i)=max(0.0_dp,p(i)-normal_quantile(0.975_dp)*res%se(i))
    res%upper(i)=min(1.0_dp,p(i)+normal_quantile(0.975_dp)*res%se(i))
    end do
    if(present(contrast))then
    allocate(c(size(contrast,1),size(contrast,2)))
    c=contrast
    else
    allocate(c(d,d))
    c=-1.0_dp/real(d,dp)
    do i=1,d
    c(i,i)=c(i,i)+1
    end do
    end if
    cc=matmul(c,transpose(c))
    call pseudoinverse_sym(cc,pinv)
    tmat=matmul(transpose(c),matmul(pinv,c))
    tv=real(n,dp)*matmul(tmat,jvar)
    denom=trace(tv)
    res%f_stat=real(n,dp)*dot_product(p,matmul(tmat,p))/denom
    call symmetric_eigenvalues(tv,eval)
    allocate(qsim(ns))
    do r=1,ns
    qsim(r)=0
    do i=1,d
    qsim(r)=qsim(r)+randn()**2*eval(i)/denom
    end do
    end do
    call sort_real(qsim)
    res%c_simu=[qsim(max(1,int(0.90_dp*ns))),qsim(max(1,int(0.95_dp*ns))),qsim(max(1,int(0.99_dp*ns)))]
    res%p_simu=real(count(qsim>res%f_stat),dp)/real(ns,dp)
    res%df=denom**2/trace(matmul(tv,tv))
    res%c_anal=[chisq_quantile(0.90_dp,res%df)/res%df,chisq_quantile(0.95_dp,res%df)/res%df,chisq_quantile(0.99_dp,res%df)/res%df]
    res%p_anal=1-chisq_cdf(res%f_stat*res%df,res%df)
  end subroutine surv_factorial

  pure function diagonal(a) result(d)
    real(dp),intent(in)::a(:,:)
    real(dp)::d(min(size(a,1),size(a,2)))
    integer::i
    do i=1,size(d)
    d(i)=a(i,i)
    end do
  end function diagonal
  pure real(dp) function trace(a) result(v)
    real(dp),intent(in)::a(:,:)
    integer::i
    v=0
    do i=1,min(size(a,1),size(a,2))
    v=v+a(i,i)
    end do
  end function trace
  subroutine sort_real(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::v
    do i=2,size(x)
    v=x(i)
    j=i-1
    do while(j>=1)
    if(x(j)<=v)exit
    x(j+1)=x(j)
    j=j-1
    end do
    x(j+1)=v
    end do
  end subroutine sort_real
end module compound_cox_factorial
