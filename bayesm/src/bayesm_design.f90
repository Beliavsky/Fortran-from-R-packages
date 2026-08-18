module bayesm_design
  use bayesm_kinds, only: dp
  implicit none
  private
  public :: create_x
contains
  subroutine create_x(p,xa,xd,xout,info,intc,diff,base)
    integer, intent(in) :: p
    real(dp), intent(in), optional :: xa(:,:),xd(:,:)
    real(dp), allocatable, intent(out) :: xout(:,:)
    integer, intent(out) :: info
    logical, intent(in), optional :: intc,diff
    integer, intent(in), optional :: base
    logical :: use_int,use_diff
    integer :: n,na,nd,b,k1,k2,row,i,j,aidx,col,alt
    real(dp), allocatable :: xd2(:,:),xone(:,:),xtwo(:,:),imod(:,:),vals(:,:)
    use_int=.true.; if (present(intc)) use_int=intc
    use_diff=.false.; if (present(diff)) use_diff=diff
    b=p; if (present(base)) b=base
    info=0
    if (.not.present(xa) .and. .not.present(xd)) then
      info=1; allocate(xout(0,0)); return
    end if
    if (present(xa)) then
      n=size(xa,1)
      if (mod(size(xa,2),p)/=0) then; info=2; allocate(xout(0,0)); return; end if
      na=size(xa,2)/p
    else
      n=size(xd,1); na=0
    end if
    if (present(xd)) then
      if (size(xd,1)/=n) then; info=3; allocate(xout(0,0)); return; end if
      nd=size(xd,2)
    else
      nd=0
    end if
    if (b<1 .or. b>p) then; info=4; allocate(xout(0,0)); return; end if
    if (use_int) then
      allocate(xd2(n,nd+1)); xd2(:,1)=1.0_dp
      if (nd>0) xd2(:,2:)=xd
    else
      allocate(xd2(n,nd)); if (nd>0) xd2=xd
    end if
    k1=size(xd2,2)*(p-1)
    if (use_diff) then
      allocate(imod(p-1,p-1)); imod=0.0_dp
      do i=1,p-1; imod(i,i)=1.0_dp; end do
      allocate(xone(n*(p-1),k1)); xone=0.0_dp
      if (size(xd2,2)>0) then
        do row=1,n
          do i=1,p-1
            do j=1,size(xd2,2)
              xone((row-1)*(p-1)+i,(j-1)*(p-1)+i)=xd2(row,j)
            end do
          end do
        end do
      end if
    else
      allocate(imod(p,p-1)); imod=0.0_dp
      col=0
      do alt=1,p
        if (alt==b) cycle
        col=col+1; imod(alt,col)=1.0_dp
      end do
      allocate(xone(n*p,k1)); xone=0.0_dp
      if (size(xd2,2)>0) then
        do row=1,n
          do alt=1,p
            do j=1,size(xd2,2)
              do i=1,p-1
                xone((row-1)*p+alt,(j-1)*(p-1)+i)=xd2(row,j)*imod(alt,i)
              end do
            end do
          end do
        end do
      end if
    end if
    if (na>0) then
      if (use_diff) then
        allocate(xtwo(n*(p-1),na)); xtwo=0.0_dp
        allocate(vals(n,p))
        do aidx=1,na
          vals=xa(:,(aidx-1)*p+1:aidx*p)
          do row=1,n
            do alt=1,p-1
              j=alt; if (j>=b) j=j+1
              xtwo((row-1)*(p-1)+alt,aidx)=vals(row,j)-vals(row,b)
            end do
          end do
        end do
      else
        allocate(xtwo(n*p,na)); xtwo=0.0_dp
        do aidx=1,na
          do row=1,n
            do alt=1,p
              xtwo((row-1)*p+alt,aidx)=xa(row,(aidx-1)*p+alt)
            end do
          end do
        end do
      end if
    else
      allocate(xtwo(size(xone,1),0))
    end if
    k2=size(xtwo,2)
    allocate(xout(size(xone,1),k1+k2))
    if (k1>0) xout(:,1:k1)=xone
    if (k2>0) xout(:,k1+1:k1+k2)=xtwo
  end subroutine create_x
end module bayesm_design
