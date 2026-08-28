module lavaan_efa
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : inverse_general, sym_eigen_jacobi, logdet_spd, trace_matrix
   use lavaan_optimizer, only : bfgs_minimize
   use gpa_rotation, only : rotation_result, rotation_options, gpforth, gpfoblq, rotate_random_starts
   use gpa_criteria, only : criterion_options
   implicit none
   private

   type, public :: efa_result
      real(dp), allocatable :: loadings(:, :), unrotated(:, :), uniqueness(:), communality(:)
      real(dp), allocatable :: phi(:, :), rotation(:, :), reproduced(:, :), residual(:, :)
      real(dp) :: objective = huge(1.0_dp)
      integer :: iterations = 0, status = 0
      logical :: converged = .false., oblique = .false.
      character(len=16) :: extraction = ''
      character(len=32) :: rotation_method = ''
   end type efa_result

   public :: efa_fit_cov, efa_principal_axis, efa_ml

contains

   subroutine efa_fit_cov(s,nfactors,result,extraction,rotation,oblique,nstarts)
      real(dp),intent(in)::s(:,:)
      integer,intent(in)::nfactors
      type(efa_result),intent(out)::result
      character(len=*),intent(in),optional::extraction,rotation
      logical,intent(in),optional::oblique
      integer,intent(in),optional::nstarts
      character(len=16)::ext
      character(len=32)::rot
      logical::obl
      integer::ns
      type(rotation_result)::rr
      type(rotation_options)::ro
      type(criterion_options)::co
      integer::i,info

      ext='PAF'
      if(present(extraction)) ext=upper(adjustl(extraction))
      rot='none'
      if(present(rotation)) rot=lower(adjustl(rotation))
      obl=.false.
      if(present(oblique)) obl=oblique
      ns=0
      if(present(nstarts)) ns=max(0,nstarts)

      select case(trim(ext))
      case('ML')
         call efa_ml(s,nfactors,result)
      case default
         call efa_principal_axis(s,nfactors,result)
      end select
      if(result%status/=0) return
      result%extraction=trim(ext)
      result%rotation_method=trim(rot)
      result%oblique=obl
      if(trim(rot)/='none' .and. trim(rot)/='') then
         ro%random_starts=ns
         if(ns>1) then
            call rotate_random_starts(result%unrotated,trim(rot),.not.obl,ns,rr,co,ro)
         else if(obl) then
            call gpfoblq(result%unrotated,trim(rot),rr,co,ro)
         else
            call gpforth(result%unrotated,trim(rot),rr,co,ro)
         end if
         if(rr%info>1) then
         result%status=20+rr%info
         return
         end if
         result%loadings=rr%loadings
         result%rotation=rr%th
         result%phi=rr%phi
         result%objective=rr%objective
      else
         result%loadings=result%unrotated
         allocate(result%rotation(nfactors,nfactors),result%phi(nfactors,nfactors))
         result%rotation=0.0_dp
         result%phi=0.0_dp
         do i=1,nfactors
         result%rotation(i,i)=1.0_dp
         result%phi(i,i)=1.0_dp
         end do
      end if
      call efa_reproduction(s,result,info)
      if(info/=0) result%status=info
   end subroutine efa_fit_cov

   subroutine efa_principal_axis(s,nfactors,result,maxiter,tol)
      real(dp),intent(in)::s(:,:)
      integer,intent(in)::nfactors
      type(efa_result),intent(out)::result
      integer,intent(in),optional::maxiter
      real(dp),intent(in),optional::tol
      real(dp),allocatable::r(:,:),rinv(:,:),work(:,:),vals(:),vecs(:,:),load(:,:),h2(:),old(:),sd(:)
      integer::p,i,j,it,mx,info
      real(dp)::eps
      p=size(s,1)
      mx=200
      if(present(maxiter)) mx=maxiter
      eps=1e-8_dp
      if(present(tol)) eps=tol
      if(size(s,2)/=p .or. nfactors<1 .or. nfactors>=p) then
      result%status=-1
      return
      end if
      allocate(r(p,p),sd(p))
      r=s
      do i=1,p
         if(s(i,i)<=0.0_dp) then
         result%status=-2
         return
         end if
         sd(i)=sqrt(s(i,i))
      end do
      do j=1,p
      do i=1,p
      r(i,j)=s(i,j)/(sd(i)*sd(j))
      end do
      end do
      call inverse_general(r,rinv,info)
      allocate(h2(p),old(p),work(p,p))
      if(info==0) then
         do i=1,p
            if(rinv(i,i)>1.0_dp) then
            h2(i)=max(0.05_dp,min(0.99_dp,1.0_dp-1.0_dp/rinv(i,i)))
            else
            h2(i)=0.5_dp
            end if
         end do
      else
         h2=0.5_dp
      end if
      do it=1,mx
         old=h2
         work=r
         do i=1,p
         work(i,i)=h2(i)
         end do
         call top_eigen_loadings(work,nfactors,load,vals,vecs,info)
         if(info/=0) then
         result%status=info
         return
         end if
         h2=sum(load*load,dim=2)
         h2=max(0.0_dp,min(0.999999_dp,h2))
         if(maxval(abs(h2-old))<eps) exit
      end do
      allocate(result%unrotated(p,nfactors),result%uniqueness(p),result%communality(p))
      do i=1,p
         result%unrotated(i,:)=sd(i)*load(i,:)
         result%communality(i)=s(i,i)*h2(i)
         result%uniqueness(i)=max(s(i,i)-result%communality(i),1e-10_dp*s(i,i))
      end do
      result%iterations=min(it,mx)
      result%converged=(it<=mx)
      result%objective=sum((r-matmul(load,transpose(load)))**2)
      result%status=0
   end subroutine efa_principal_axis

   subroutine efa_ml(s,nfactors,result,maxiter,tol)
      real(dp),intent(in)::s(:,:)
      integer,intent(in)::nfactors
      type(efa_result),intent(out)::result
      integer,intent(in),optional::maxiter
      real(dp),intent(in),optional::tol
      real(dp),allocatable::x(:),psi(:),load(:,:)
      real(dp)::fval,eps,lds,ldm
      logical::conv
      integer::p,i,info,it,mx
      p=size(s,1)
      mx=1000
      if(present(maxiter)) mx=maxiter
      eps=1e-7_dp
      if(present(tol)) eps=tol
      if(size(s,2)/=p .or. nfactors<1 .or. nfactors>=p) then
      result%status=-1
      return
      end if
      allocate(x(p))
      do i=1,p
      x(i)=log(max(0.3_dp*s(i,i),1.0e-8_dp))
      end do
      call bfgs_minimize(obj,x,fval,conv,it,maxiter=mx,tol=eps)
      call concentrated(x,psi,load,info)
      if(info/=0) then
      result%status=info
      return
      end if
      allocate(result%unrotated(p,nfactors),result%uniqueness(p),result%communality(p))
      result%unrotated=load
      result%uniqueness=psi
      result%communality=sum(load*load,dim=2)
      result%objective=fval
      result%iterations=it
      result%converged=conv
      result%status=0
   contains
      function obj(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         real(dp),allocatable::pp(:),ll(:,:),ss(:,:),sii(:,:)
         integer::istat,jj
         call concentrated(z,pp,ll,istat)
         if(istat/=0) then
         v=huge(1.0_dp)/100.0_dp
         return
         end if
         allocate(ss(p,p))
         ss=matmul(ll,transpose(ll))
         do jj=1,p
         ss(jj,jj)=ss(jj,jj)+pp(jj)
         end do
         call inverse_general(ss,sii,istat)
         if(istat/=0) then
         v=huge(1.0_dp)/100.0_dp
         return
         end if
         ldm=logdet_spd(ss,istat)
         lds=logdet_spd(s,istat)
         if(istat/=0) then
         v=huge(1.0_dp)/100.0_dp
         return
         end if
         v=ldm+trace_matrix(matmul(s,sii))-lds-real(p,dp)
      end function obj
      subroutine concentrated(z,pp,ll,istat)
         real(dp),intent(in)::z(:)
         real(dp),allocatable,intent(out)::pp(:),ll(:,:)
         integer,intent(out)::istat
         real(dp),allocatable::d(:),c(:,:),ev(:),evec(:,:)
         integer::ii,jj,kk
         allocate(pp(p),d(p),c(p,p))
         do ii=1,p
            pp(ii)=max(1e-8_dp*s(ii,ii),min(exp(z(ii)),0.999999_dp*s(ii,ii)))
            d(ii)=sqrt(pp(ii))
         end do
         do jj=1,p
         do ii=1,p
         c(ii,jj)=s(ii,jj)/(d(ii)*d(jj))
         end do
         end do
         call sym_eigen_jacobi(c,ev,evec,istat)
         if(istat/=0) return
         call sort_desc(ev,evec)
         allocate(ll(p,nfactors))
         ll=0.0_dp
         do kk=1,nfactors
            if(ev(kk)>1.0_dp) then
               do ii=1,p
               ll(ii,kk)=d(ii)*evec(ii,kk)*sqrt(ev(kk)-1.0_dp)
               end do
            end if
         end do
      end subroutine concentrated
   end subroutine efa_ml

   subroutine top_eigen_loadings(a,k,load,vals,vecs,info)
      real(dp),intent(in)::a(:,:)
      integer,intent(in)::k
      real(dp),allocatable,intent(out)::load(:,:),vals(:),vecs(:,:)
      integer,intent(out)::info
      real(dp),allocatable::v(:),q(:,:)
      integer::j,p
      p=size(a,1)
      call sym_eigen_jacobi(a,v,q,info)
      if(info/=0) return
      call sort_desc(v,q)
      allocate(load(p,k),vals(k),vecs(p,k))
      vals=v(1:k)
      vecs=q(:,1:k)
      load=0.0_dp
      do j=1,k
         if(vals(j)>0.0_dp) load(:,j)=vecs(:,j)*sqrt(vals(j))
      end do
   end subroutine top_eigen_loadings

   subroutine sort_desc(values,vectors)
      real(dp),intent(inout)::values(:),vectors(:,:)
      integer::i,j
      real(dp)::t
      real(dp),allocatable::col(:)
      allocate(col(size(vectors,1)))
      do i=1,size(values)-1
      do j=i+1,size(values)
         if(values(j)>values(i)) then
         t=values(i)
         values(i)=values(j)
         values(j)=t
            col=vectors(:,i)
            vectors(:,i)=vectors(:,j)
            vectors(:,j)=col
         end if
      end do
      end do
   end subroutine sort_desc

   subroutine efa_reproduction(s,result,info)
      real(dp),intent(in)::s(:,:)
      type(efa_result),intent(inout)::result
      integer,intent(out)::info
      integer::p,i
      p=size(s,1)
      allocate(result%reproduced(p,p),result%residual(p,p))
      if(result%oblique) then
         result%reproduced=matmul(result%loadings,matmul(result%phi,transpose(result%loadings)))
      else
         result%reproduced=matmul(result%loadings,transpose(result%loadings))
      end if
      do i=1,p
      result%reproduced(i,i)=result%reproduced(i,i)+result%uniqueness(i)
      end do
      result%residual=s-result%reproduced
      info=0
   end subroutine efa_reproduction


   pure function lower(s) result(t)
      character(len=*),intent(in)::s
      character(len=len(s))::t
      integer::i,c
      do i=1,len(s)
      c=iachar(s(i:i))
      if(c>=65.and.c<=90) then
      t(i:i)=achar(c+32)
      else
      t(i:i)=s(i:i)
      end if
      end do
   end function lower
   pure function upper(s) result(t)
      character(len=*),intent(in)::s
      character(len=len(s))::t
      integer::i,c
      do i=1,len(s)
      c=iachar(s(i:i))
      if(c>=97.and.c<=122) then
      t(i:i)=achar(c-32)
      else
      t(i:i)=s(i:i)
      end if
      end do
   end function upper
end module lavaan_efa
