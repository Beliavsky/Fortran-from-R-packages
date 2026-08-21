! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_geometry
  use dirichletreg_kinds, only : dp
  implicit none
  private
  public :: to_ternary, to_quaternary
contains
  subroutine to_ternary(comp,xy,stat)
    real(dp),intent(in)::comp(:,:)
    real(dp),intent(out)::xy(:,:)
    integer,intent(out),optional::stat
    real(dp),parameter::sqrt3=1.7320508075688772935274463415058723669_dp
    if(present(stat)) stat=0
    if(size(comp,2)/=3 .or. size(xy,1)/=size(comp,1) .or. size(xy,2)/=2 .or. &
       any(comp<0.0_dp) .or. any(comp>1.0_dp) .or. any(abs(sum(comp,dim=2)-1.0_dp)>1.0e-10_dp)) then
      if(present(stat)) stat=1; xy=0.0_dp; return
    end if
    xy(:,1)=(comp(:,1)+2.0_dp*comp(:,3))/sqrt3
    xy(:,2)=comp(:,1)
  end subroutine to_ternary

  subroutine to_quaternary(comp,xyz,stat)
    real(dp),intent(in)::comp(:,:)
    real(dp),intent(out)::xyz(:,:)
    integer,intent(out),optional::stat
    real(dp),parameter::sqrt3=1.7320508075688772935274463415058723669_dp
    if(present(stat)) stat=0
    if(size(comp,2)/=4 .or. size(xyz,1)/=size(comp,1) .or. size(xyz,2)/=3 .or. &
       any(comp<0.0_dp) .or. any(comp>1.0_dp) .or. any(abs(sum(comp,dim=2)-1.0_dp)>1.0e-10_dp)) then
      if(present(stat)) stat=1; xyz=0.0_dp; return
    end if
    xyz(:,1)=(comp(:,1)+2.0_dp*comp(:,3)+comp(:,4))/sqrt3
    xyz(:,2)=comp(:,1)+comp(:,4)/3.0_dp
    xyz(:,3)=comp(:,4)
  end subroutine to_quaternary
end module dirichletreg_geometry
