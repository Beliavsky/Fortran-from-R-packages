module relsurv_ratetable
  use relsurv_kinds, only : dp
  implicit none
  private

  type, public :: ratetable_type
    integer :: ndim = 0
    integer, allocatable :: dims(:)
    integer, allocatable :: factor(:)
    integer, allocatable :: ncuts(:)
    real(dp), allocatable :: cuts(:,:)
    real(dp), allocatable :: rate(:)
  end type ratetable_type

  type, public :: net_summary_type
    real(dp), allocatable :: yidli(:), yidsi(:), dnisi(:), yisi(:)
    real(dp), allocatable :: yidlisi(:), sidli(:), yi(:), dnisisq(:)
    real(dp), allocatable :: yisisq(:), dni(:), sis(:), yisidli(:)
    real(dp), allocatable :: yisis(:), sit(:), yisitt(:), yidlisitt(:)
    real(dp), allocatable :: yidlisiw(:), dci(:)
  end type net_summary_type

  public :: make_ratetable, pystep, pystep2, expected_survival
  public :: population_survival_curve, netwei_summary, netfast_summary
  public :: population_hazard_increment
  public :: expprep2_expected, expprep2_summary

contains
  function make_ratetable(dims, factor, cuts, ncuts, rate) result(tab)
    integer, intent(in) :: dims(:), factor(:), ncuts(:)
    real(dp), intent(in) :: cuts(:,:), rate(:)
    type(ratetable_type) :: tab
    tab%ndim=size(dims)
    allocate(tab%dims(tab%ndim),tab%factor(tab%ndim),tab%ncuts(tab%ndim))
    tab%dims=dims; tab%factor=factor; tab%ncuts=ncuts
    allocate(tab%cuts(size(cuts,1),size(cuts,2))); tab%cuts=cuts
    allocate(tab%rate(size(rate))); tab%rate=rate
  end function make_ratetable

  subroutine pystep(tab, data, step, edge, amount, index, index2, wt)
    type(ratetable_type), intent(in) :: tab
    real(dp), intent(in) :: data(:), step
    integer, intent(in) :: edge
    real(dp), intent(out) :: amount, wt
    integer, intent(out) :: index, index2
    integer :: i,j,kk,dtemp,juse,lin0,lin2
    real(dp) :: maxtime,shortfall,temp
    kk=1; lin0=0; lin2=0; wt=1.0_dp; shortfall=0.0_dp; maxtime=step
    do i=1,tab%ndim
      if (tab%factor(i)==1) then
        juse=max(0,min(tab%dims(i)-1,nint(data(i))-1))
        lin0=lin0+juse*kk
      else
        dtemp=tab%ncuts(i)
        j=1
        do while(j<=dtemp)
          if (data(i)<tab%cuts(j,i)) exit
          j=j+1
        end do
        if (j==1) then
          temp=tab%cuts(1,i)-data(i)
          if (edge==0 .and. temp>shortfall) shortfall=min(temp,step)
          maxtime=min(maxtime,temp)
          juse=0
        else if (j>dtemp) then
          if (edge==0) then
            temp=tab%cuts(dtemp,i)-data(i)
            if (temp<=0.0_dp) then
              shortfall=step
            else
              maxtime=min(maxtime,temp)
            end if
          end if
          if (tab%factor(i)>1) then
            juse=tab%dims(i)-1
          else
            juse=tab%dims(i)-1
          end if
        else
          temp=tab%cuts(j,i)-data(i)
          maxtime=min(maxtime,temp)
          juse=j-2
          if (tab%factor(i)>1) then
            wt=1.0_dp-real(mod(juse,tab%factor(i)),dp)/real(tab%factor(i),dp)
            juse=juse/tab%factor(i)
            lin2=kk
          end if
        end if
        lin0=lin0+juse*kk
      end if
      kk=kk*tab%dims(i)
    end do
    index=lin0+1; index2=lin0+lin2+1
    if (shortfall==0.0_dp) then
      amount=max(0.0_dp,maxtime)
    else
      index=0; index2=0; amount=shortfall
    end if
  end subroutine pystep

  subroutine pystep2(tab,data,index,index2,wt)
    type(ratetable_type), intent(in) :: tab
    real(dp), intent(in) :: data(:)
    integer, intent(out) :: index,index2
    real(dp), intent(out) :: wt
    integer :: i,j,kk,juse,lin0
    kk=1; lin0=0; wt=1.0_dp
    do i=1,tab%ndim
      if(tab%factor(i)==1) then
        juse=max(0,min(tab%dims(i)-1,nint(data(i))-1))
      else
        j=1
        do while(j<=tab%ncuts(i))
          if(data(i)<tab%cuts(j,i)) exit
          j=j+1
        end do
        if(j/=1) j=j-1
        juse=max(0,min(tab%dims(i)-1,j-1))
      end if
      lin0=lin0+juse*kk; kk=kk*tab%dims(i)
    end do
    index=lin0+1; index2=index
  end subroutine pystep2

  function expected_survival(tab,x,time) result(surv)
    type(ratetable_type), intent(in) :: tab
    real(dp), intent(in) :: x(:,:), time(:)
    real(dp) :: surv(size(time))
    real(dp) :: data(tab%ndim),remain,dt,wt,haz
    integer :: i,k,idx,idx2
    if(size(x,1)/=size(time) .or. size(x,2)/=tab%ndim) then
      surv=0.0_dp; return
    end if
    do i=1,size(time)
      data=x(i,:); remain=max(0.0_dp,time(i)); haz=0.0_dp
      do while(remain>100.0_dp*epsilon(1.0_dp))
        call pystep(tab,data,remain,1,dt,idx,idx2,wt)
        if(dt<=0.0_dp .or. idx<1) exit
        if(wt<1.0_dp) then
          haz=haz+dt*(wt*tab%rate(idx)+(1.0_dp-wt)*tab%rate(idx2))
        else
          haz=haz+dt*tab%rate(idx)
        end if
        do k=1,tab%ndim
          if(tab%factor(k)/=1) data(k)=data(k)+dt
        end do
        remain=remain-dt
      end do
      surv(i)=exp(-haz)
    end do
  end function expected_survival

  subroutine population_survival_curve(tab,x,times,surv)
    type(ratetable_type), intent(in) :: tab
    real(dp), intent(in) :: x(:,:),times(:)
    real(dp), intent(out) :: surv(:,:)
    real(dp), allocatable :: tvec(:)
    integer :: j
    allocate(tvec(size(x,1)))
    do j=1,size(times)
      tvec=times(j)
      surv(:,j)=expected_survival(tab,x,tvec)
    end do
  end subroutine population_survival_curve

  subroutine init_summary(out,ntime,n)
    type(net_summary_type), intent(out) :: out
    integer,intent(in)::ntime,n
    allocate(out%yidli(ntime),out%yidsi(ntime),out%dnisi(ntime),out%yisi(ntime), &
      out%yidlisi(ntime),out%sidli(ntime),out%yi(ntime),out%dnisisq(ntime), &
      out%yisisq(ntime),out%dni(ntime),out%sis(ntime),out%yisidli(ntime), &
      out%yisis(ntime),out%sit(n),out%yisitt(ntime),out%yidlisitt(ntime), &
      out%yidlisiw(ntime),out%dci(ntime))
    out%yidli=0.0_dp; out%yidsi=0.0_dp; out%dnisi=0.0_dp; out%yisi=0.0_dp
    out%yidlisi=0.0_dp; out%sidli=0.0_dp; out%yi=0.0_dp; out%dnisisq=0.0_dp
    out%yisisq=0.0_dp; out%dni=0.0_dp; out%sis=0.0_dp; out%yisidli=0.0_dp
    out%yisis=0.0_dp; out%sit=0.0_dp; out%yisitt=0.0_dp
    out%yidlisitt=0.0_dp; out%yidlisiw=0.0_dp; out%dci=0.0_dp
  end subroutine init_summary

  subroutine netwei_summary(tab,x,y,status,times,out)
    type(ratetable_type), intent(in) :: tab
    real(dp), intent(in) :: x(:,:),y(:),times(:)
    integer,intent(in)::status(:)
    type(net_summary_type), intent(out)::out
    real(dp),allocatable::si(:)
    real(dp)::data(tab%ndim),hazard,remain,dt,wt,thiscell,time0,rate0,den
    integer::i,j,k,idx,idx2,n
    n=size(y); call init_summary(out,size(times),n); allocate(si(n)); si=1.0_dp
    time0=0.0_dp
    do j=1,size(times)
      thiscell=times(j)-time0
      do i=1,n
        data=x(i,:)
        do k=1,tab%ndim
          if(tab%factor(k)/=1) data(k)=data(k)+time0
        end do
        remain=max(0.0_dp,thiscell); hazard=0.0_dp
        do while(remain>100.0_dp*epsilon(1.0_dp))
          call pystep(tab,data,remain,1,dt,idx,idx2,wt)
          if(dt<=0.0_dp .or. idx<1) exit
          rate0=max(tab%rate(idx),1.0e-12_dp)
          if(wt<1.0_dp) then
            hazard=hazard+dt*(wt*rate0+(1.0_dp-wt)*tab%rate(idx2))
          else
            hazard=hazard+dt*rate0
          end if
          do k=1,tab%ndim
            if(tab%factor(k)/=1) data(k)=data(k)+dt
          end do
          remain=remain-dt
        end do
        if(thiscell>0.0_dp) then
          den=hazard/thiscell
          if(abs(den)>tiny(1.0_dp)) out%sit(i)=out%sit(i)+si(i)*(1.0_dp-exp(-hazard))/den
        end if
        si(i)=si(i)*exp(-hazard)
        out%sis(j)=out%sis(j)+si(i); out%sidli(j)=out%sidli(j)+hazard*si(i)
        if(y(i)>=times(j)) then
          out%yidsi(j)=out%yidsi(j)+exp(-hazard); out%yidli(j)=out%yidli(j)+hazard
          out%yisidli(j)=out%yisidli(j)+hazard*si(i); out%yi(j)=out%yi(j)+1.0_dp
          out%yisi(j)=out%yisi(j)+1.0_dp/si(i); out%yisisq(j)=out%yisisq(j)+1.0_dp/(si(i)*si(i))
          out%yisis(j)=out%yisis(j)+si(i); out%yidlisi(j)=out%yidlisi(j)+hazard/si(i)
          if(y(i)==times(j)) then
            out%dnisi(j)=out%dnisi(j)+real(status(i),dp)/si(i)
            out%dni(j)=out%dni(j)+real(status(i),dp)
            out%dci(j)=out%dci(j)+real(1-status(i),dp)
            out%dnisisq(j)=out%dnisisq(j)+real(status(i),dp)/(si(i)*si(i))
          end if
        end if
      end do
      time0=times(j)
    end do
  end subroutine netwei_summary

  subroutine netfast_summary(tab,x,y,ystart,status,times,precision,out)
    type(ratetable_type), intent(in) :: tab
    real(dp), intent(in) :: x(:,:),y(:),ystart(:),times(:),precision
    integer,intent(in)::status(:)
    type(net_summary_type),intent(out)::out
    real(dp),allocatable::si(:),sitt(:)
    real(dp)::data(tab%ndim),time0,thiscell,ftime,dt,tstart,wt,lambda1,lambda2
    real(dp)::fyisi,fyisi2,fyidlisi,fyidlisi2,fint,sisum,sisumtt,hazard
    integer::i,j,k,idx,idx2,n,jfine
    n=size(y); call init_summary(out,size(times),n); allocate(si(n),sitt(n)); si=1.0_dp
    time0=0.0_dp
    do j=1,size(times)
      thiscell=times(j)-time0; ftime=0.0_dp; fyisi=0.0_dp; fyisi2=0.0_dp
      fyidlisi=0.0_dp; fyidlisi2=0.0_dp; fint=0.0_dp; hazard=0.0_dp
      sitt=si; jfine=0; sisum=0.0_dp; sisumtt=0.0_dp
      do while(ftime<thiscell-100.0_dp*epsilon(1.0_dp))
        jfine=jfine+1; out%dnisi(j)=0.0_dp; out%dni(j)=0.0_dp; out%dci(j)=0.0_dp; out%dnisisq(j)=0.0_dp
        tstart=time0+ftime; dt=min(max(precision,epsilon(1.0_dp)),thiscell-ftime)
        sisum=0.0_dp; sisumtt=0.0_dp
        do i=1,n
          if(y(i)>=times(j)) then
            data=x(i,:)
            do k=1,tab%ndim
              if(tab%factor(k)/=1) data(k)=data(k)+tstart
            end do
            call pystep2(tab,data,idx,idx2,wt)
            lambda1=tab%rate(idx); lambda2=tab%rate(idx2)
            if(ystart(i)<times(j)) then
              fyidlisi=fyidlisi+lambda1/si(i)
              fyidlisi2=fyidlisi2+lambda1/(si(i)*exp(-dt*lambda1))
              fyisi=fyisi+1.0_dp/si(i)
              fyisi2=fyisi2+1.0_dp/(si(i)*exp(-dt*lambda1))
              if(wt<1.0_dp) then
                hazard=hazard+dt*(wt*lambda1+(1.0_dp-wt)*lambda2)
              else
                hazard=hazard+dt*lambda1
              end if
            end if
            si(i)=si(i)*exp(-dt*lambda1)
            if(ystart(i)<=times(j)) then
              sisum=sisum+1.0_dp/si(i); sisumtt=sisumtt+1.0_dp/sitt(i)
              if(jfine==1) out%yi(j)=out%yi(j)+1.0_dp
            end if
            if(y(i)==times(j)) then
              out%dnisi(j)=out%dnisi(j)+real(status(i),dp)/si(i)
              out%dni(j)=out%dni(j)+real(status(i),dp)
              out%dci(j)=out%dci(j)+real(1-status(i),dp)
              out%dnisisq(j)=out%dnisisq(j)+real(status(i),dp)/(si(i)*si(i))
            end if
          end if
        end do
        if(fyisi>0.0_dp .and. fyisi2>0.0_dp) fint=fint+0.5_dp*(fyidlisi/fyisi+fyidlisi2/fyisi2)*dt
        ftime=ftime+dt
      end do
      out%yisi(j)=sisum; out%yisitt(j)=sisumtt
      if(sisum>0.0_dp) out%yidlisi(j)=hazard/sisum
      if(sisumtt>0.0_dp) out%yidlisitt(j)=hazard/sisumtt
      out%yidlisiw(j)=fint; out%yidli(j)=hazard
      time0=times(j)
    end do
  end subroutine netfast_summary

  subroutine population_hazard_increment(tab,x,t0,t1,dh)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::x(:,:),t0,t1
    real(dp),intent(out)::dh(:)
    real(dp)::data(tab%ndim),remain,dt,wt,haz
    integer::i,k,idx,idx2
    do i=1,size(x,1)
      data=x(i,:)
      do k=1,tab%ndim
        if(tab%factor(k)/=1) data(k)=data(k)+t0
      end do
      remain=max(0.0_dp,t1-t0); haz=0.0_dp
      do while(remain>100.0_dp*epsilon(1.0_dp))
        call pystep(tab,data,remain,1,dt,idx,idx2,wt)
        if(dt<=0.0_dp .or. idx<1) exit
        if(wt<1.0_dp) haz=haz+dt*(wt*tab%rate(idx)+(1.0_dp-wt)*tab%rate(idx2))
        if(wt>=1.0_dp) haz=haz+dt*tab%rate(idx)
        do k=1,tab%ndim
          if(tab%factor(k)/=1) data(k)=data(k)+dt
        end do
        remain=remain-dt
      end do
      dh(i)=haz
    end do
  end subroutine population_hazard_increment
  function expprep2_expected(tab,x,time) result(surv)
    type(ratetable_type), intent(in) :: tab
    real(dp), intent(in) :: x(:,:), time(:)
    real(dp) :: surv(size(time))
    surv=expected_survival(tab,x,time)
  end function expprep2_expected

  subroutine expprep2_summary(tab,x,y,status,times,out,fast,ystart,precision)
    type(ratetable_type), intent(in) :: tab
    real(dp), intent(in) :: x(:,:),y(:),times(:)
    integer, intent(in) :: status(:)
    type(net_summary_type), intent(out) :: out
    logical, intent(in), optional :: fast
    real(dp), intent(in), optional :: ystart(:),precision
    logical :: use_fast
    real(dp), allocatable :: ys(:)
    real(dp) :: prec
    use_fast=.false.
    if(present(fast)) use_fast=fast
    if(use_fast) then
      allocate(ys(size(y)))
      ys=0.0_dp
      if(present(ystart)) ys=ystart
      prec=1.0_dp
      if(present(precision)) prec=precision
      call netfast_summary(tab,x,y,ys,status,times,prec,out)
    else
      call netwei_summary(tab,x,y,status,times,out)
    end if
  end subroutine expprep2_summary

end module relsurv_ratetable
