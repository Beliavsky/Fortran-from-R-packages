! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_design
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, LONELY_REMOVE
  implicit none
  private
  public :: make_design, design_degf, design_weights, validate_design
contains

  subroutine make_design(weight,cluster,design,strata,samp_size,pop_size,lonely_psu,ultimate_cluster,adjust_domain_lonely)
    real(dp), intent(in) :: weight(:)
    integer, intent(in) :: cluster(:,:)
    type(survey_design_t), intent(out) :: design
    integer, intent(in), optional :: strata(:,:)
    real(dp), intent(in), optional :: samp_size(:,:), pop_size(:,:)
    integer, intent(in), optional :: lonely_psu
    logical, intent(in), optional :: ultimate_cluster, adjust_domain_lonely
    integer :: n,s,i,j,k,st,ct,npsu
    n=size(weight); s=size(cluster,2)
    if(size(cluster,1)/=n) error stop 'make_design: cluster rows must match weights'
    design%n=n; design%stages=s
    allocate(design%weight(n),design%cluster(n,s),design%strata(n,s), &
             design%samp_size(n,s),design%pop_size(n,s))
    design%weight=weight; design%cluster=cluster
    design%strata=1
    if(present(strata)) then
      if(size(strata,1)/=n .or. size(strata,2)/=s) error stop 'make_design: strata shape mismatch'
      design%strata=strata
      design%has_strata=any(strata/=1)
    end if
    if(present(samp_size)) then
      if(any(shape(samp_size)/=shape(design%samp_size))) error stop 'make_design: samp_size shape mismatch'
      design%samp_size=samp_size
    else
      do j=1,s
        do i=1,n
          st=design%strata(i,j); npsu=0
          do ct=1,n
            if(design%strata(ct,j)==st) then
              if(.not.any([(design%cluster(k,j)==design%cluster(ct,j) .and. design%strata(k,j)==st,k=1,ct-1)])) npsu=npsu+1
            end if
          end do
          design%samp_size(i,j)=real(npsu,dp)
        end do
      end do
    end if
    if(present(pop_size)) then
      if(any(shape(pop_size)/=shape(design%pop_size))) error stop 'make_design: pop_size shape mismatch'
      design%pop_size=pop_size
    else
      design%pop_size=huge(1.0_dp)
    end if
    design%lonely_psu=LONELY_REMOVE
    if(present(lonely_psu)) design%lonely_psu=lonely_psu
    design%ultimate_cluster=.false.; if(present(ultimate_cluster)) design%ultimate_cluster=ultimate_cluster
    design%adjust_domain_lonely=.false.; if(present(adjust_domain_lonely)) design%adjust_domain_lonely=adjust_domain_lonely
    call validate_design(design)
  end subroutine make_design

  subroutine validate_design(design)
    type(survey_design_t), intent(in) :: design
    if(design%n<=0 .or. design%stages<=0) error stop 'survey design must contain observations and stages'
    if(.not.allocated(design%weight)) error stop 'survey design has no weights'
    if(any(design%weight<0.0_dp)) error stop 'survey weights must be nonnegative'
    if(any(design%samp_size<1.0_dp)) error stop 'sample PSU counts must be >= 1'
    if(any(design%pop_size<design%samp_size .and. design%pop_size<huge(1.0_dp)/2)) &
      error stop 'finite population size smaller than sampled PSU count'
  end subroutine validate_design

  function design_weights(design) result(w)
    type(survey_design_t), intent(in) :: design
    real(dp) :: w(design%n)
    w=design%weight
  end function design_weights

  integer function design_degf(design) result(df)
    type(survey_design_t), intent(in) :: design
    integer, allocatable :: strata_vals(:), cluster_vals(:)
    integer :: h,nstrata,npsu
    call unique_int(design%strata(:,1),strata_vals)
    nstrata=size(strata_vals); npsu=0
    do h=1,nstrata
      call unique_int(pack(design%cluster(:,1),design%strata(:,1)==strata_vals(h)),cluster_vals)
      npsu=npsu+size(cluster_vals)
    end do
    df=max(0,npsu-nstrata)
  end function design_degf

  subroutine unique_int(x,u)
    integer, intent(in) :: x(:)
    integer, allocatable, intent(out) :: u(:)
    integer, allocatable :: tmp(:)
    integer :: i,n
    allocate(tmp(size(x))); n=0
    do i=1,size(x)
      if(n==0 .or. .not.any(tmp(1:n)==x(i))) then
        n=n+1; tmp(n)=x(i)
      end if
    end do
    allocate(u(n)); if(n>0) u=tmp(1:n)
  end subroutine unique_int
end module survey_design
