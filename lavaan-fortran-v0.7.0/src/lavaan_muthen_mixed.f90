module lavaan_muthen_mixed
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : inverse_general
   use lavaan_optimizer, only : bfgs_minimize
   use lavaan_ordinal, only : ordinal_thresholds, bvn_rectangle
   implicit none
   private

   integer, parameter, public :: lavaan_numeric=0, lavaan_ordered=1

   type, public :: muthen_mixed_result
      real(dp), allocatable :: intercept(:), thresholds(:,:), slopes(:,:), variance(:)
      real(dp), allocatable :: correlation(:,:), covariance(:,:), stats(:), gamma(:,:), weight(:,:)
      integer, allocatable :: ncat(:), stat_offset(:)
      integer :: nobs=0, nvar=0, nexo=0, status=0
      logical :: gamma_computed=.false.
   end type muthen_mixed_result

   public :: muthen1984_mixed

contains

   subroutine muthen1984_mixed(data,ov_type,ncat,result,exo,compute_gamma)
      real(dp),intent(in)::data(:,:)
      integer,intent(in)::ov_type(:),ncat(:)
      type(muthen_mixed_result),intent(out)::result
      real(dp),intent(in),optional::exo(:,:)
      logical,intent(in),optional::compute_gamma
      logical::dogamma
      integer::n,k,i,info
      real(dp),allocatable::jk(:,:),meanjk(:),d(:),g(:,:),ginv(:,:)
      type(muthen_mixed_result)::ri

      dogamma=.true.
      if(present(compute_gamma)) dogamma=compute_gamma
      call mixed_core(data,ov_type,ncat,result,exo)
      if(result%status/=0 .or. .not.dogamma) return
      n=size(data,1)
      k=size(result%stats)
      if(n<8 .or. k==0) then
      result%status=9001
      return
      end if
      allocate(jk(n,k))
      do i=1,n
         if(present(exo)) then
            call mixed_core(pack_rows(data,i),ov_type,ncat,ri,pack_rows(exo,i))
         else
            call mixed_core(pack_rows(data,i),ov_type,ncat,ri)
         end if
         if(ri%status/=0 .or. size(ri%stats)/=k) then
         result%status=9100+i
         return
         end if
         jk(i,:)=ri%stats
      end do
      meanjk=sum(jk,dim=1)/real(n,dp)
      allocate(g(k,k))
      g=0.0_dp
      do i=1,n
         d=jk(i,:)-meanjk
         g=g+spread(d,2,k)*spread(d,1,k)
      end do
      ! lavaan Gamma is N times the asymptotic covariance of the statistic vector.
      result%gamma=real(n-1,dp)*g
      g=result%gamma
      do i=1,k
      g(i,i)=g(i,i)+1.0e-10_dp*max(1.0_dp,abs(g(i,i)))
      end do
      call inverse_general(g,ginv,info)
      allocate(result%weight(k,k))
      result%weight=0.0_dp
      if(info==0) then
         result%weight=ginv
      else
         do i=1,k
            if(result%gamma(i,i)>tiny(1.0_dp)) result%weight(i,i)=1.0_dp/result%gamma(i,i)
         end do
      end if
      result%gamma_computed=.true.
   end subroutine muthen1984_mixed

   subroutine mixed_core(data,ov_type,ncat,result,exo)
      real(dp),intent(in)::data(:,:)
      integer,intent(in)::ov_type(:),ncat(:)
      type(muthen_mixed_result),intent(out)::result
      real(dp),intent(in),optional::exo(:,:)
      real(dp),allocatable::x(:,:),resid(:,:),sd(:),th(:),b(:),rho(:,:)
      integer::n,p,nexo,maxth,j,i,k,info,pos,nstat
      n=size(data,1)
      p=size(data,2)
      nexo=0
      if(present(exo)) nexo=size(exo,2)
      result%nobs=n
      result%nvar=p
      result%nexo=nexo
      result%ncat=ncat
      if(size(ov_type)/=p .or. size(ncat)/=p .or. n<4) then
      result%status=-1
      return
      end if
      if(present(exo)) then
         if(size(exo,1)/=n) then
         result%status=-2
         return
         end if
         allocate(x(n,nexo+1))
         x(:,1)=1.0_dp
         if(nexo>0) x(:,2:)=exo
      else
         allocate(x(n,1))
         x=1.0_dp
      end if
      maxth=max(1,maxval(ncat)-1)
      allocate(result%intercept(p),result%thresholds(maxth,p),result%slopes(p,nexo),result%variance(p))
      allocate(resid(n,p),sd(p))
      result%intercept=0.0_dp
      result%thresholds=0.0_dp
      result%slopes=0.0_dp
      result%variance=1.0_dp
      resid=0.0_dp
      sd=1.0_dp
      do j=1,p
         select case(ov_type(j))
         case(lavaan_numeric)
            call fit_linear(data(:,j),x,b,result%variance(j),info)
            if(info/=0) then
            result%status=100+j
            return
            end if
            result%intercept(j)=b(1)
            if(nexo>0) result%slopes(j,:)=b(2:)
            resid(:,j)=data(:,j)-matmul(x,b)
            sd(j)=sqrt(max(result%variance(j),1.0e-12_dp))
         case(lavaan_ordered)
            if(ncat(j)<2) then
            result%status=200+j
            return
            end if
            call fit_ordinal_probit(nint(data(:,j)),x,ncat(j),th,b,info)
            if(info/=0) then
            result%status=300+j
            return
            end if
            result%thresholds(1:ncat(j)-1,j)=th
            if(nexo>0) result%slopes(j,:)=b(2:)
            result%intercept(j)=0.0_dp
            result%variance(j)=1.0_dp
            sd(j)=1.0_dp
         case default
            result%status=400+j
            return
         end select
      end do
      allocate(rho(p,p))
      rho=0.0_dp
      do j=1,p
      rho(j,j)=1.0_dp
      end do
      do j=1,p-1
         do i=j+1,p
            call pair_rho(i,j,data,ov_type,ncat,result,x,resid,sd,rho(i,j),info)
            if(info/=0) then
            result%status=500+10*j+i
            return
            end if
            rho(j,i)=rho(i,j)
         end do
      end do
      result%correlation=rho
      allocate(result%covariance(p,p))
      result%covariance=0.0_dp
      do j=1,p
         do i=1,p
            result%covariance(i,j)=rho(i,j)*sd(i)*sd(j)
         end do
      end do
      allocate(result%stat_offset(p+1))
      result%stat_offset(1)=1
      nstat=0
      do j=1,p
         if(ov_type(j)==lavaan_numeric) then
            nstat=nstat+1+nexo+1
         else
            nstat=nstat+(ncat(j)-1)+nexo
         end if
         result%stat_offset(j+1)=nstat+1
      end do
      nstat=nstat+p*(p-1)/2
      allocate(result%stats(nstat))
      pos=0
      do j=1,p
         if(ov_type(j)==lavaan_numeric) then
            pos=pos+1
            result%stats(pos)=result%intercept(j)
            if(nexo>0) then
            result%stats(pos+1:pos+nexo)=result%slopes(j,:)
            pos=pos+nexo
            end if
            pos=pos+1
            result%stats(pos)=result%variance(j)
         else
            k=ncat(j)-1
            result%stats(pos+1:pos+k)=result%thresholds(1:k,j)
            pos=pos+k
            if(nexo>0) then
            result%stats(pos+1:pos+nexo)=result%slopes(j,:)
            pos=pos+nexo
            end if
         end if
      end do
      do j=1,p-1
      do i=j+1,p
      pos=pos+1
      result%stats(pos)=rho(i,j)
      end do
      end do
      result%status=0
   end subroutine mixed_core

   subroutine fit_linear(y,x,b,var,info)
      real(dp),intent(in)::y(:),x(:,:)
      real(dp),allocatable,intent(out)::b(:)
      real(dp),intent(out)::var
      integer,intent(out)::info
      real(dp),allocatable::xtx(:,:),inv(:,:),r(:)
      xtx=matmul(transpose(x),x)
      call inverse_general(xtx,inv,info)
      if(info/=0) then
      allocate(b(size(x,2)))
      b=0.0_dp
      var=huge(1.0_dp)
      return
      end if
      b=matmul(inv,matmul(transpose(x),y))
      r=y-matmul(x,b)
      var=sum(r*r)/real(max(1,size(y)-size(x,2)),dp)
   end subroutine fit_linear

   subroutine fit_ordinal_probit(y,x,nc,threshold,b,info)
      integer,intent(in)::y(:),nc
      real(dp),intent(in)::x(:,:)
      real(dp),allocatable,intent(out)::threshold(:),b(:)
      integer,intent(out)::info
      integer,allocatable::cnt(:)
      real(dp),allocatable::t0(:),par(:)
      real(dp)::fval
      integer::i,k,it
      logical::conv
      allocate(cnt(nc))
      cnt=0
      do i=1,size(y)
         if(y(i)<1 .or. y(i)>nc) then
         info=-1
         return
         end if
         cnt(y(i))=cnt(y(i))+1
      end do
      if(any(cnt==0)) then
      info=-2
      return
      end if
      t0=ordinal_thresholds(cnt)
      k=nc-1
      allocate(par(k+size(x,2)-1))
      par=0.0_dp
      par(1)=t0(1)
      do i=2,k
      par(i)=log(max(t0(i)-t0(i-1),1.0e-3_dp))
      end do
      call bfgs_minimize(nll,par,fval,conv,it,maxiter=500,tol=2.0e-7_dp)
      call unpack_cut(par(1:k),threshold)
      allocate(b(size(x,2)))
      b=0.0_dp
      if(size(x,2)>1) b(2:)=par(k+1:)
      info=merge(0,1,conv .or. fval<huge(1.0_dp)/100.0_dp)
   contains
      function nll(v) result(f)
         real(dp),intent(in)::v(:)
         real(dp)::f,lo,hi,eta,pr
         real(dp),allocatable::tt(:)
         integer::r,c
         call unpack_cut(v(1:k),tt)
         f=0.0_dp
         do r=1,size(y)
            c=y(r)
            eta=0.0_dp
            if(size(x,2)>1) eta=dot_product(x(r,2:),v(k+1:))
            if(c==1) then
            lo=-huge(1.0_dp)
            else
            lo=tt(c-1)-eta
            end if
            if(c==nc) then
            hi=huge(1.0_dp)
            else
            hi=tt(c)-eta
            end if
            pr=max(normal_cdf(hi)-normal_cdf(lo),1.0e-300_dp)
            f=f-log(pr)
         end do
      end function nll
   end subroutine fit_ordinal_probit

   subroutine unpack_cut(v,t)
      real(dp),intent(in)::v(:)
      real(dp),allocatable,intent(out)::t(:)
      integer::i
      allocate(t(size(v)))
      t(1)=v(1)
      do i=2,size(v)
      t(i)=t(i-1)+exp(v(i))
      end do
   end subroutine unpack_cut

   subroutine pair_rho(i,j,data,typ,ncat,result,x,resid,sd,rho,info)
      integer,intent(in)::i,j,typ(:),ncat(:)
      real(dp),intent(in)::data(:,:),x(:,:),resid(:,:),sd(:)
      type(muthen_mixed_result),intent(in)::result
      real(dp),intent(out)::rho
      integer,intent(out)::info
      real(dp)::a,b,x1,x2,f1,f2,gr
      integer::iter
      a=-0.995_dp
      b=0.995_dp
      gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
      x1=b-gr*(b-a)
      x2=a+gr*(b-a)
      f1=negll(x1)
      f2=negll(x2)
      do iter=1,120
         if(abs(b-a)<2.0e-7_dp) exit
         if(f1>f2) then
         a=x1
         x1=x2
         f1=f2
         x2=a+gr*(b-a)
         f2=negll(x2)
         else
         b=x2
         x2=x1
         f2=f1
         x1=b-gr*(b-a)
         f1=negll(x1)
         end if
      end do
      rho=0.5_dp*(a+b)
      info=0
   contains
      function negll(rr) result(v)
         real(dp),intent(in)::rr
         real(dp)::v,eta1,eta2,l1,u1,l2,u2,pr,zv,den
         integer::r,c1,c2
         v=0.0_dp
         den=sqrt(max(1.0e-10_dp,1.0_dp-rr*rr))
         if(typ(i)==lavaan_numeric .and. typ(j)==lavaan_numeric) then
            v=0.5_dp*real(size(data,1),dp)*log(max(1.0e-12_dp,1.0_dp-rr*rr))
            do r=1,size(data,1)
               v=v+0.5_dp*( (resid(r,i)/sd(i))**2 - 2.0_dp*rr*(resid(r,i)/sd(i))*(resid(r,j)/sd(j)) + &
                  (resid(r,j)/sd(j))**2 )/max(1.0e-12_dp,1.0_dp-rr*rr)
            end do
            return
         end if
         do r=1,size(data,1)
            if(typ(i)==lavaan_ordered) then
               eta1=0.0_dp
               if(size(x,2)>1) eta1=dot_product(x(r,2:),result%slopes(i,:))
               c1=nint(data(r,i))
               call bounds(c1,ncat(i),result%thresholds(:,i),eta1,l1,u1)
            end if
            if(typ(j)==lavaan_ordered) then
               eta2=0.0_dp
               if(size(x,2)>1) eta2=dot_product(x(r,2:),result%slopes(j,:))
               c2=nint(data(r,j))
               call bounds(c2,ncat(j),result%thresholds(:,j),eta2,l2,u2)
            end if
            if(typ(i)==lavaan_ordered .and. typ(j)==lavaan_ordered) then
               pr=bvn_rectangle(l1,u1,l2,u2,rr)
            else if(typ(i)==lavaan_numeric) then
               zv=resid(r,i)/sd(i)
               pr=max(normal_cdf((u2-rr*zv)/den)-normal_cdf((l2-rr*zv)/den),1.0e-300_dp)
            else
               zv=resid(r,j)/sd(j)
               pr=max(normal_cdf((u1-rr*zv)/den)-normal_cdf((l1-rr*zv)/den),1.0e-300_dp)
            end if
            v=v-log(max(pr,1.0e-300_dp))
         end do
      end function negll
   end subroutine pair_rho

   subroutine bounds(c,nc,t,eta,lo,hi)
      integer,intent(in)::c,nc
      real(dp),intent(in)::t(:),eta
      real(dp),intent(out)::lo,hi
      if(c==1) then
      lo=-huge(1.0_dp)
      else
      lo=t(c-1)-eta
      end if
      if(c==nc) then
      hi=huge(1.0_dp)
      else
      hi=t(c)-eta
      end if
   end subroutine bounds

   pure function normal_cdf(z) result(p)
      real(dp),intent(in)::z
      real(dp)::p
      if(z>8.0_dp) then
      p=1.0_dp
      else if(z< -8.0_dp) then
      p=0.0_dp
      else
      p=0.5_dp*erfc(-z/sqrt(2.0_dp))
      end if
   end function normal_cdf

   function pack_rows(a,drop) result(b)
      real(dp),intent(in)::a(:,:)
      integer,intent(in)::drop
      real(dp),allocatable::b(:,:)
      integer::n,i,pos
      n=size(a,1)
      allocate(b(n-1,size(a,2)))
      pos=0
      do i=1,n
         if(i==drop) cycle
         pos=pos+1
         b(pos,:)=a(i,:)
      end do
   end function pack_rows

end module lavaan_muthen_mixed
