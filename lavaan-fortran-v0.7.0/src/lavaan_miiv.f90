module lavaan_miiv
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : inverse_general, sample_mean_cov
   implicit none
   private

   type, public :: miiv_result
      real(dp), allocatable :: beta(:), vcov(:, :), fitted(:), residual(:)
      real(dp), allocatable :: first_stage_f(:)
      real(dp) :: sigma2=huge(1.0_dp), sargan=0.0_dp
      integer :: sargan_df=0, rank=0, status=0
   end type miiv_result

   type, public :: miiv_equation
      integer :: outcome_node=0
      integer, allocatable :: predictor_nodes(:), instrument_nodes(:)
      logical :: identified=.false.
   end type miiv_equation

   public :: miiv_2sls, miiv_2sls_cov, ram_miiv_candidates, ram_miiv_equations

contains

   subroutine miiv_2sls(y,x,z,result,intercept)
      real(dp),intent(in)::y(:),x(:,:),z(:,:)
      type(miiv_result),intent(out)::result
      logical,intent(in),optional::intercept
      real(dp),allocatable::xx(:,:),zz(:,:),ztz(:,:),ztzi(:,:),pxx(:,:),pxy(:),inv(:,:),zhat(:,:)
      real(dp),allocatable::gamma(:),zr(:)
      logical::icpt
      integer::n,k,l,info,j
      real(dp)::sst,sse,r2
      n=size(y)
      k=size(x,2)
      l=size(z,2)
      icpt=.true.
      if(present(intercept)) icpt=intercept
      if(size(x,1)/=n .or. size(z,1)/=n .or. l<k) then
      result%status=-1
      return
      end if
      if(icpt) then
         allocate(xx(n,k+1),zz(n,l+1))
         xx(:,1)=1.0_dp
         xx(:,2:)=x
         zz(:,1)=1.0_dp
         zz(:,2:)=z
      else
         xx=x
         zz=z
      end if
      k=size(xx,2)
      l=size(zz,2)
      ztz=matmul(transpose(zz),zz)
      call inverse_general(ztz,ztzi,info)
      if(info/=0) then
      result%status=info
      return
      end if
      zhat=matmul(zz,matmul(ztzi,matmul(transpose(zz),xx)))
      pxx=matmul(transpose(xx),zhat)
      pxy=matmul(transpose(zhat),y)
      call inverse_general(pxx,inv,info)
      if(info/=0) then
      result%status=info
      return
      end if
      result%beta=matmul(inv,pxy)
      result%fitted=matmul(xx,result%beta)
      result%residual=y-result%fitted
      result%rank=k
      result%sigma2=dot_product(result%residual,result%residual)/real(max(1,n-k),dp)
      allocate(result%vcov(k,k))
      result%vcov=result%sigma2*inv
      allocate(result%first_stage_f(k-merge(1,0,icpt)))
      result%first_stage_f=0.0_dp
      do j=1,size(result%first_stage_f)
         call first_stage_stat(xx(:,j+merge(1,0,icpt)),zz,result%first_stage_f(j))
      end do
      ! Sargan over-identification statistic: n*R^2 from residuals on all instruments.
      gamma=matmul(ztzi,matmul(transpose(zz),result%residual))
      zr=matmul(zz,gamma)
      sst=dot_product(result%residual,result%residual)
      sse=dot_product(result%residual-zr,result%residual-zr)
      r2=max(0.0_dp,min(1.0_dp,1.0_dp-sse/max(sst,tiny(1.0_dp))))
      result%sargan=real(n,dp)*r2
      result%sargan_df=max(0,l-k)
      result%status=0
   contains
      subroutine first_stage_stat(v,inst,f)
         real(dp),intent(in)::v(:),inst(:,:)
         real(dp),intent(out)::f
         real(dp),allocatable::ii(:,:),coef(:),fit(:)
         real(dp)::ss0,ss1
         integer::istat,df1,df2
         call inverse_general(matmul(transpose(inst),inst),ii,istat)
         if(istat/=0) then
         f=0.0_dp
         return
         end if
         coef=matmul(ii,matmul(transpose(inst),v))
         fit=matmul(inst,coef)
         ss1=dot_product(v-fit,v-fit)
         ss0=dot_product(v-sum(v)/real(size(v),dp),v-sum(v)/real(size(v),dp))
         df1=max(1,size(inst,2)-1)
         df2=max(1,size(v)-size(inst,2))
         f=max(0.0_dp,(ss0-ss1)/real(df1,dp))/(ss1/real(df2,dp)+tiny(1.0_dp))
      end subroutine first_stage_stat
   end subroutine miiv_2sls

   subroutine miiv_2sls_cov(sxx,sxz,szz,sxy,szy,beta,vcov,nobs,status)
      real(dp),intent(in)::sxx(:,:),sxz(:,:),szz(:,:),sxy(:),szy(:)
      real(dp),allocatable,intent(out)::beta(:),vcov(:,:)
      integer,intent(in)::nobs
      integer,intent(out)::status
      real(dp),allocatable::zzi(:,:),m(:,:),mi(:,:),rhs(:)
      real(dp)::sigma2
      integer::k,info
      k=size(sxx,1)
      if(size(sxx,2)/=k .or. size(sxz,1)/=k .or. size(sxz,2)/=size(szz,1) .or. &
         size(szz,1)/=size(szz,2) .or. size(sxy)/=k .or. size(szy)/=size(szz,1)) then
         status=-1
         return
      end if
      call inverse_general(szz,zzi,info)
      if(info/=0) then
      status=info
      return
      end if
      m=matmul(sxz,matmul(zzi,transpose(sxz)))
      rhs=matmul(sxz,matmul(zzi,szy))
      call inverse_general(m,mi,info)
      if(info/=0) then
      status=info
      return
      end if
      beta=matmul(mi,rhs)
      sigma2=max(1.0e-12_dp,1.0_dp-2.0_dp*dot_product(beta,sxy)+dot_product(beta,matmul(sxx,beta)))
      vcov=sigma2*mi/real(max(1,nobs),dp)
      status=0
   end subroutine miiv_2sls_cov

   subroutine ram_miiv_equations(model,equations,status,tol,observed_only)
      use lavaan_ram, only : ram_model
      use lavaan_linalg, only : inverse_general
      type(ram_model),intent(in)::model
      type(miiv_equation),allocatable,intent(out)::equations(:)
      integer,intent(out)::status
      real(dp),intent(in),optional::tol
      logical,intent(in),optional::observed_only
      real(dp),allocatable::ia(:,:),inv(:,:),covall(:,:),errcov(:,:)
      integer,allocatable::pred(:),cand(:),inst(:),tmpout(:)
      real(dp)::eps
      integer::n,i,j,m,neq,info,node,ninst
      logical::obs_only
      n=size(model%a,1)
      eps=1.0e-8_dp
      if(present(tol)) eps=tol
      obs_only=.true.
      if(present(observed_only)) obs_only=observed_only
      if(size(model%a,2)/=n .or. any(shape(model%s)/=[n,n])) then
         status=-1
         allocate(equations(0))
         return
      end if
      allocate(ia(n,n))
      ia=-model%a
      do i=1,n
      ia(i,i)=ia(i,i)+1.0_dp
      end do
      call inverse_general(ia,inv,info)
      if(info/=0) then
      status=info
      allocate(equations(0))
      return
      end if
      covall=matmul(inv,matmul(model%s,transpose(inv)))
      errcov=matmul(inv,model%s)
      allocate(tmpout(n))
      neq=0
      do i=1,n
         pred=pack([(j,j=1,n)],abs(model%a(i,:))>eps)
         if(size(pred)>0) then
         neq=neq+1
         tmpout(neq)=i
         end if
      end do
      allocate(equations(neq))
      do m=1,neq
         i=tmpout(m)
         pred=pack([(j,j=1,n)],abs(model%a(i,:))>eps)
         if(obs_only .and. allocated(model%observed)) then
            cand=model%observed
         else
            cand=[(j,j=1,n)]
         end if
         allocate(inst(size(cand)))
         ninst=0
         do j=1,size(cand)
            node=cand(j)
            if(node==i .or. any(node==pred)) cycle
            ! Exclude descendants of the dependent variable and variables carrying its disturbance.
            if(abs(inv(node,i))>eps) cycle
            if(abs(errcov(node,i))>eps) cycle
            ! An instrument must be relevant for at least one endogenous predictor.
            if(maxval(abs(covall(node,pred)))<=eps) cycle
            ninst=ninst+1
            inst(ninst)=node
         end do
         equations(m)%outcome_node=i
         equations(m)%predictor_nodes=pred
         allocate(equations(m)%instrument_nodes(ninst))
         if(ninst>0) equations(m)%instrument_nodes=inst(1:ninst)
         equations(m)%identified=ninst>=size(pred)
         deallocate(inst)
      end do
      status=0
   end subroutine ram_miiv_equations

   subroutine ram_miiv_candidates(model,outcome_node,predictor_nodes,candidate_nodes,instruments,status,tol)
      use lavaan_ram, only : ram_model
      use lavaan_linalg, only : inverse_general
      type(ram_model),intent(in)::model
      integer,intent(in)::outcome_node,predictor_nodes(:),candidate_nodes(:)
      integer,allocatable,intent(out)::instruments(:)
      integer,intent(out)::status
      real(dp),intent(in),optional::tol
      real(dp),allocatable::ia(:,:),inv(:,:),xe(:,:),covall(:,:)
      integer,allocatable::tmp(:)
      real(dp)::eps
      integer::n,i,node,m,info
      n=size(model%a,1)
      eps=1.0e-8_dp
      if(present(tol)) eps=tol
      if(outcome_node<1 .or. outcome_node>n .or. any(predictor_nodes<1) .or. any(predictor_nodes>n) .or. &
         any(candidate_nodes<1) .or. any(candidate_nodes>n)) then
         status=-1
         allocate(instruments(0))
         return
      end if
      allocate(ia(n,n))
      ia=-model%a
      do i=1,n
      ia(i,i)=ia(i,i)+1.0_dp
      end do
      call inverse_general(ia,inv,info)
      if(info/=0) then
      status=info
      allocate(instruments(0))
      return
      end if
      xe=matmul(inv,model%s)
      covall=matmul(inv,matmul(model%s,transpose(inv)))
      allocate(tmp(size(candidate_nodes)))
      m=0
      do i=1,size(candidate_nodes)
         node=candidate_nodes(i)
         if(node==outcome_node .or. any(node==predictor_nodes)) cycle
         if(abs(xe(node,outcome_node))>eps) cycle
         if(maxval(abs(covall(node,predictor_nodes)))<=eps) cycle
         m=m+1
         tmp(m)=node
      end do
      allocate(instruments(m))
      if(m>0) instruments=tmp(1:m)
      status=0
   end subroutine ram_miiv_candidates
end module lavaan_miiv
