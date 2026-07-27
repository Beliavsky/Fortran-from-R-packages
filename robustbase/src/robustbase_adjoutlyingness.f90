! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_adjoutlyingness
   use robustbase_kinds, only: dp
   use robustbase_sort, only: median, quantile_type7
   use robustbase_medcouple, only: medcouple
   use robustbase_linalg, only: solve_linear, matrix_rank
   implicit none
   private
   public :: adjusted_outlyingness_result, adjusted_outlyingness_full

   type :: adjusted_outlyingness_result
      real(dp),allocatable :: outlyingness(:)
      real(dp),allocatable :: center(:)
      real(dp),allocatable :: directions(:,:)
      logical,allocatable :: non_outlier(:)
      real(dp) :: cutoff=0.0_dp
      real(dp) :: medcouple_outlyingness=0.0_dp
      integer :: iterations=0
      integer :: directions_used=0
      logical :: converged=.false.
   end type adjusted_outlyingness_result
contains
   subroutine adjusted_outlyingness_full(x,result,n_directions,p_sample,clower,cupper,alpha_cutoff,coef,max_iterations,central)
      real(dp),intent(in)::x(:,:)
      type(adjusted_outlyingness_result),intent(out)::result
      integer,intent(in),optional::n_directions,p_sample,max_iterations
      real(dp),intent(in),optional::clower,cupper,alpha_cutoff,coef
      logical,intent(in),optional::central
      integer::n,p,ndir,ps,maxit,it,d,i,j,info,used,pair_i,pair_j
      integer,allocatable::ind(:)
      real(dp),allocatable::dirs(:,:),direction(:),pmat(:,:),rhs(:),projection(:),centered(:),q1(:),q3(:),iqr(:),mc(:),upper(:),lower(:),yup(:),ylo(:),all_proj(:,:)
      real(dp)::cl,cu,alpha,k,normd,nx,qlo,qhi,u,denom
      logical::cent
      n=size(x,1);p=size(x,2)
      if(n<1 .or. p<1)error stop 'adjusted_outlyingness_full: empty input'
      ndir=250;if(present(n_directions))ndir=max(1,n_directions)
      ps=p;if(present(p_sample))ps=max(p,p_sample)
      ps=min(ps,n)
      cl=4.0_dp;if(present(clower))cl=clower
      cu=3.0_dp;if(present(cupper))cu=cupper
      alpha=0.75_dp;if(present(alpha_cutoff))alpha=alpha_cutoff
      k=1.5_dp;if(present(coef))k=coef
      maxit=max(100,p)*ndir;if(present(max_iterations))maxit=max(ndir,max_iterations)
      cent=.false.;if(present(central))cent=central
      allocate(dirs(p,ndir),direction(p),projection(n),centered(n),q1(ndir),q3(ndir),iqr(ndir),mc(ndir),upper(ndir),lower(ndir),yup(ndir),ylo(ndir),all_proj(n,ndir))
      dirs=0.0_dp;used=0;it=0
      nx=max(sum(abs(x))/real(n*p,dp),1.0e-14_dp)
      if(p<=n)then
         allocate(ind(ps),pmat(ps,p),rhs(ps))
         rhs=1.0_dp
         do while(used<ndir .and. it<maxit)
            it=it+1
            call random_subset(n,ps,ind)
            pmat=x(ind,:)
            if(matrix_rank(pmat)<p)cycle
            if(ps==p)then
               call solve_linear(pmat,rhs,direction,info)
            else
               call least_squares_direction(pmat,rhs,direction,info)
            end if
            if(info/=0)cycle
            normd=sqrt(sum(direction*direction))
            if(normd*nx<=1.0e-12_dp)cycle
            direction=direction/normd
            if(is_duplicate(direction,dirs,used))cycle
            used=used+1;dirs(:,used)=direction
         end do
      else
         pair_i=1;pair_j=2
         do while(used<ndir .and. pair_i<n)
            direction=x(pair_j,:)-x(pair_i,:)
            normd=sqrt(sum(direction*direction))
            if(normd>1.0e-12_dp)then
               direction=direction/normd
               if(.not.is_duplicate(direction,dirs,used))then
                  used=used+1;dirs(:,used)=direction
               end if
            end if
            pair_j=pair_j+1
            if(pair_j>n)then;pair_i=pair_i+1;pair_j=pair_i+1;end if
         end do
         it=used
      end if
      if(used==0)error stop 'adjusted_outlyingness_full: no valid direction'
      allocate(result%outlyingness(n),result%center(p),result%directions(p,used),result%non_outlier(n))
      result%directions=dirs(:,1:used)
      do j=1,p
         result%center(j)=median(x(:,j))
      end do
      do d=1,used
         projection=matmul(x,dirs(:,d))
         centered=projection-median(projection)
         all_proj(:,d)=centered
         q1(d)=quantile_type7(centered,0.25_dp);q3(d)=quantile_type7(centered,0.75_dp);iqr(d)=max(q3(d)-q1(d),1.0e-14_dp)
         if(cent)then
            mc(d)=0.0_dp
         else
            mc(d)=medcouple(centered)
         end if
         if(mc(d)>=0.0_dp)then
            upper(d)=q3(d)+k*iqr(d)*exp(cu*mc(d))
            lower(d)=q1(d)-k*iqr(d)*exp(-cl*mc(d))
         else
            upper(d)=q3(d)+k*iqr(d)*exp(cl*mc(d))
            lower(d)=q1(d)-k*iqr(d)*exp(-cu*mc(d))
         end if
         yup(d)=-huge(1.0_dp);ylo(d)=huge(1.0_dp)
         do i=1,n
            if(centered(i)<upper(d))yup(d)=max(yup(d),centered(i))
            if(centered(i)>lower(d))ylo(d)=min(ylo(d),centered(i))
         end do
         yup(d)=max(yup(d),1.0e-14_dp)
         ylo(d)=max(-ylo(d),1.0e-14_dp)
      end do
      result%outlyingness=0.0_dp
      do i=1,n
         do d=1,used
            u=all_proj(i,d)
            if(abs(u)<=1.0e-15_dp)cycle
            if(u>0.0_dp)then;denom=yup(d);else;denom=ylo(d);end if
            if(denom>1.0e-14_dp)result%outlyingness(i)=max(result%outlyingness(i),abs(u)/denom)
         end do
      end do
      qlo=quantile_type7(result%outlyingness,1.0_dp-alpha)
      qhi=quantile_type7(result%outlyingness,alpha)
      result%medcouple_outlyingness=medcouple(result%outlyingness)
      result%cutoff=qhi+k*(qhi-qlo)*merge(exp(cu*result%medcouple_outlyingness),1.0_dp,result%medcouple_outlyingness>0.0_dp)
      result%non_outlier=result%outlyingness<=result%cutoff
      result%iterations=it;result%directions_used=used;result%converged=used==ndir .or. p>n
   end subroutine adjusted_outlyingness_full

   subroutine least_squares_direction(a,b,x,info)
      use robustbase_linalg, only: least_squares
      real(dp),intent(in)::a(:,:),b(:)
      real(dp),intent(out)::x(:)
      integer,intent(out)::info
      call least_squares(a,b,x,info)
   end subroutine least_squares_direction

   subroutine random_subset(n,k,subset)
      integer,intent(in)::n,k
      integer,intent(out)::subset(:)
      integer::m,candidate
      real(dp)::u
      m=0
      do while(m<k)
         call random_number(u);candidate=min(n,1+int(u*real(n,dp)))
         if(m==0 .or. .not.any(subset(1:m)==candidate))then;m=m+1;subset(m)=candidate;end if
      end do
      call sort_integer(subset)
   end subroutine random_subset

   subroutine sort_integer(a)
      integer,intent(inout)::a(:)
      integer::i,j,key
      do i=2,size(a)
         key=a(i);j=i-1
         do while(j>=1)
            if(a(j)<=key)exit
            a(j+1)=a(j);j=j-1
         end do
         a(j+1)=key
      end do
   end subroutine sort_integer

   logical function is_duplicate(direction,directions,n_used) result(duplicate)
      real(dp),intent(in)::direction(:),directions(:,:)
      integer,intent(in)::n_used
      integer::j
      duplicate=.false.
      do j=1,n_used
         if(maxval(abs(direction-directions(:,j)))<1.0e-10_dp .or. maxval(abs(direction+directions(:,j)))<1.0e-10_dp)then
            duplicate=.true.;return
         end if
      end do
   end function is_duplicate
end module robustbase_adjoutlyingness
