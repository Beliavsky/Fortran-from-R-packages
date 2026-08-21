! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_multiframe
  use survey_kinds, only : dp
  use survey_types, only : multiframe_design_t, svystat_t
  use survey_pps, only : ht_variance
  implicit none
  private
  public :: make_multiframe_constant, make_multiframe_expected, multiframe_total, multiframe_mean, multiframe_variance
contains
  subroutine make_multiframe_constant(weight1,weight2,overlap1_other,overlap2_other,design,theta)
    real(dp),intent(in) :: weight1(:),weight2(:)
    logical,intent(in) :: overlap1_other(:),overlap2_other(:)
    type(multiframe_design_t),intent(out) :: design
    real(dp),intent(in),optional :: theta
    real(dp) :: scale1,scale2,t
    integer :: i
    if(size(overlap1_other)/=size(weight1).or.size(overlap2_other)/=size(weight2)) error stop 'make_multiframe_constant: shape mismatch'
    design%n1=size(weight1);design%n2=size(weight2)
    allocate(design%frame_weight1(design%n1),design%frame_weight2(design%n2),design%design_weight1(design%n1),design%design_weight2(design%n2))
    design%design_weight1=weight1;design%design_weight2=weight2
    if(present(theta)) then
      t=theta;if(t<0.0_dp.or.t>1.0_dp) error stop 'make_multiframe_constant: theta outside [0,1]';scale1=t;scale2=1.0_dp-t
    else
      scale1=sum(weight1)/real(max(1,size(weight1)),dp);scale2=sum(weight2)/real(max(1,size(weight2)),dp)
      t=scale1+scale2;if(t<=0.0_dp) error stop 'make_multiframe_constant: nonpositive mean weights';scale1=scale1/t;scale2=scale2/t
    end if
    do i=1,design%n1;design%frame_weight1(i)=merge(scale1,1.0_dp,overlap1_other(i));end do
    do i=1,design%n2;design%frame_weight2(i)=merge(scale2,1.0_dp,overlap2_other(i));end do
  end subroutine make_multiframe_constant

  subroutine make_multiframe_expected(weight1,weight2,overlap1,overlap2,design,overlaps_are_weights)
    real(dp),intent(in) :: weight1(:),weight2(:),overlap1(:,:),overlap2(:,:)
    type(multiframe_design_t),intent(out) :: design
    logical,intent(in),optional :: overlaps_are_weights
    logical :: asw
    if(size(overlap1,1)/=size(weight1).or.size(overlap1,2)/=2.or.size(overlap2,1)/=size(weight2).or.size(overlap2,2)/=2) &
      error stop 'make_multiframe_expected: overlap matrices must be n x 2'
    asw=.true.;if(present(overlaps_are_weights))asw=overlaps_are_weights
    design%n1=size(weight1);design%n2=size(weight2)
    allocate(design%frame_weight1(design%n1),design%frame_weight2(design%n2),design%design_weight1(design%n1),design%design_weight2(design%n2))
    design%design_weight1=weight1;design%design_weight2=weight2
    call expected_frame_weights(weight1,overlap1,asw,design%frame_weight1)
    call expected_frame_weights(weight2,overlap2,asw,design%frame_weight2)
  end subroutine make_multiframe_expected

  function multiframe_total(x1,x2,design,dcheck1,dcheck2) result(ans)
    real(dp),intent(in) :: x1(:,:),x2(:,:),dcheck1(:,:),dcheck2(:,:)
    type(multiframe_design_t),intent(in) :: design
    type(svystat_t) :: ans
    real(dp),allocatable :: z1(:,:),z2(:,:)
    integer :: p,j
    call validate_inputs(x1,x2,design,dcheck1,dcheck2);p=size(x1,2);allocate(z1(design%n1,p),z2(design%n2,p),ans%estimate(p),ans%variance(p,p),ans%influence(design%n1+design%n2,p))
    z1=x1*spread(design%design_weight1*design%frame_weight1,2,p);z2=x2*spread(design%design_weight2*design%frame_weight2,2,p)
    do j=1,p;ans%estimate(j)=sum(z1(:,j))+sum(z2(:,j));end do
    ans%variance=multiframe_variance(z1,z2,dcheck1,dcheck2);ans%influence(1:design%n1,:)=z1;ans%influence(design%n1+1:,:)=z2
  end function multiframe_total

  function multiframe_mean(x1,x2,design,dcheck1,dcheck2) result(ans)
    real(dp),intent(in) :: x1(:,:),x2(:,:),dcheck1(:,:),dcheck2(:,:)
    type(multiframe_design_t),intent(in) :: design
    type(svystat_t) :: ans
    real(dp),allocatable :: z1(:,:),z2(:,:)
    real(dp) :: sw
    integer :: p,j
    call validate_inputs(x1,x2,design,dcheck1,dcheck2);p=size(x1,2);sw=sum(design%design_weight1*design%frame_weight1)+sum(design%design_weight2*design%frame_weight2)
    if(sw<=0.0_dp) error stop 'multiframe_mean: nonpositive total weight'
    allocate(z1(design%n1,p),z2(design%n2,p),ans%estimate(p),ans%variance(p,p),ans%influence(design%n1+design%n2,p))
    do j=1,p
      ans%estimate(j)=(dot_product(design%design_weight1*design%frame_weight1,x1(:,j))+dot_product(design%design_weight2*design%frame_weight2,x2(:,j)))/sw
      z1(:,j)=(x1(:,j)-ans%estimate(j))*(design%design_weight1*design%frame_weight1)/sw
      z2(:,j)=(x2(:,j)-ans%estimate(j))*(design%design_weight2*design%frame_weight2)/sw
    end do
    ans%variance=multiframe_variance(z1,z2,dcheck1,dcheck2);ans%influence(1:design%n1,:)=z1;ans%influence(design%n1+1:,:)=z2
  end function multiframe_mean

  function multiframe_variance(z1,z2,dcheck1,dcheck2) result(v)
    real(dp),intent(in) :: z1(:,:),z2(:,:),dcheck1(:,:),dcheck2(:,:)
    real(dp) :: v(size(z1,2),size(z1,2))
    if(size(z2,2)/=size(z1,2)) error stop 'multiframe_variance: column mismatch'
    v=ht_variance(z1,dcheck1)+ht_variance(z2,dcheck2)
  end function multiframe_variance

  subroutine expected_frame_weights(w,overlap,as_weights,fw)
    real(dp),intent(in) :: w(:),overlap(:,:);logical,intent(in)::as_weights;real(dp),intent(out)::fw(:)
    real(dp) :: den,v
    integer :: i,j
    do i=1,size(w)
      den=0.0_dp
      do j=1,2
        if(overlap(i,j)>0.0_dp) then
          if(as_weights) then;v=overlap(i,j);else;v=1.0_dp/overlap(i,j);end if
          den=den+1.0_dp/v
        end if
      end do
      if(den<=0.0_dp.or.w(i)<=0.0_dp) error stop 'make_multiframe_expected: invalid weights/overlap'
      fw(i)=(1.0_dp/den)/w(i)
    end do
  end subroutine expected_frame_weights

  subroutine validate_inputs(x1,x2,d,d1,d2)
    real(dp),intent(in)::x1(:,:),x2(:,:),d1(:,:),d2(:,:);type(multiframe_design_t),intent(in)::d
    if(size(x1,1)/=d%n1.or.size(x2,1)/=d%n2.or.size(x1,2)/=size(x2,2)) error stop 'multiframe estimator: data shape'
    if(any(shape(d1)/=[d%n1,d%n1]).or.any(shape(d2)/=[d%n2,d%n2])) error stop 'multiframe estimator: Dcheck shape'
    if(.not.allocated(d%design_weight1).or..not.allocated(d%design_weight2)) error stop 'multiframe estimator: incomplete design'
  end subroutine validate_inputs
end module survey_multiframe
