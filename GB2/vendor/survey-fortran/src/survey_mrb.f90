! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_mrb
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, rep_design_t
  implicit none
  private
  public :: make_mrb
contains
  subroutine make_mrb(design,nrep,rep)
    type(survey_design_t),intent(in) :: design
    integer,intent(in) :: nrep
    type(rep_design_t),intent(out) :: rep
    real(dp),allocatable :: stage_weight(:,:),cumff(:)
    logical,allocatable :: kept(:),selected(:)
    integer,allocatable :: ustrata(:),units(:),chosen(:)
    integer :: rr,s,h,i,j,ns,nstar,nstages
    real(dp) :: fpc,lambda,cff
    if(nrep<2) error stop 'make_mrb: need at least two replicates'
    nstages=design%stages
    if(all(design%pop_size>=huge(1.0_dp)/4.0_dp)) nstages=1
    rep%n=design%n;rep%r=nrep;rep%scale=1.0_dp;rep%mse=.false.
    allocate(rep%weight(rep%n),rep%repweights(rep%n,nrep),rep%rscales(nrep));rep%weight=design%weight;rep%rscales=1.0_dp/real(nrep-1,dp)
    allocate(stage_weight(rep%n,nstages),cumff(rep%n),kept(rep%n),selected(rep%n))
    do rr=1,nrep
      stage_weight=1.0_dp;cumff=1.0_dp;kept=.true.
      do s=1,nstages
        call unique_int(design%strata(:,s),ustrata)
        do h=1,size(ustrata)
          call unique_int(pack(design%cluster(:,s),(design%strata(:,s)==ustrata(h)).and.kept),units)
          ns=size(units);nstar=ns/2;selected=.false.
          if(ns>0) then
            i=find_first((design%strata(:,s)==ustrata(h)).and.kept)
            cff=cumff(i)
            if(design%pop_size(i,s)>=huge(1.0_dp)/4.0_dp) then;fpc=0.0_dp;else;fpc=design%samp_size(i,s)/design%pop_size(i,s);end if
          else
            cff=1.0_dp;fpc=0.0_dp
          end if
          if(nstar>0) then
            call sample_without_replacement(units,nstar,chosen)
            do i=1,rep%n
              if(design%strata(i,s)==ustrata(h)) selected(i)=any(chosen==design%cluster(i,s))
            end do
            lambda=sqrt(max(0.0_dp,cff*real(nstar,dp)*(1.0_dp-fpc)/real(ns-nstar,dp)))
            do i=1,rep%n
              if(design%strata(i,s)==ustrata(h)) then
                stage_weight(i,s)=stage_weight(i,s)*(-lambda+lambda*real(ns,dp)/real(nstar,dp)*merge(1.0_dp,0.0_dp,selected(i)))
              end if
            end do
          else
            do i=1,rep%n
              if(design%strata(i,s)==ustrata(h)) stage_weight(i,s)=0.0_dp
            end do
          end if
          if(s<nstages .and. nstar>0) then
            do i=1,rep%n
              if(design%strata(i,s)==ustrata(h).and.kept(i)) then
                do j=s+1,nstages
                  stage_weight(i,j)=stage_weight(i,j)*sqrt(real(ns,dp)/real(nstar,dp))*merge(1.0_dp,0.0_dp,selected(i))
                end do
              end if
            end do
          end if
          do i=1,rep%n
            if(design%strata(i,s)==ustrata(h)) then
              kept(i)=kept(i).and.selected(i)
              cumff(i)=cumff(i)*fpc
            end if
          end do
        end do
      end do
      do i=1,rep%n;rep%repweights(i,rr)=design%weight(i)*(1.0_dp+sum(stage_weight(i,:)));end do
    end do
  end subroutine make_mrb

  subroutine sample_without_replacement(pool,k,out)
    integer,intent(in) :: pool(:),k
    integer,allocatable,intent(out) :: out(:)
    integer,allocatable :: tmp(:)
    integer :: i,j,t,n
    real(dp) :: u
    n=size(pool);if(k<0.or.k>n) error stop 'sample_without_replacement: invalid k';allocate(tmp(n));tmp=pool
    do i=1,k
      call random_number(u);j=i+int(u*real(n-i+1,dp));if(j>n)j=n;t=tmp(i);tmp(i)=tmp(j);tmp(j)=t
    end do
    allocate(out(k));if(k>0)out=tmp(1:k)
  end subroutine sample_without_replacement

  subroutine unique_int(x,u)
    integer,intent(in)::x(:);integer,allocatable,intent(out)::u(:);integer,allocatable::t(:);integer::i,n
    allocate(t(size(x)));n=0
    do i=1,size(x);if(n==0.or..not.any(t(1:n)==x(i)))then;n=n+1;t(n)=x(i);end if;end do
    allocate(u(n));if(n>0)u=t(1:n)
  end subroutine unique_int

  integer function find_first(mask) result(k)
    logical,intent(in)::mask(:);integer::i;k=0;do i=1,size(mask);if(mask(i))then;k=i;return;end if;end do
  end function find_first
end module survey_mrb
