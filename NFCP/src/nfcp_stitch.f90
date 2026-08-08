module nfcp_stitch
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use nfcp_types, only : dp, nfcp_ok, nfcp_invalid_input
  implicit none
  private
  public :: stitch_contract_numbers, stitch_by_maturity

contains

  subroutine stitch_contract_numbers(futures, contract_numbers, stitched, selected_columns, status)
    real(dp), intent(in) :: futures(:,:)
    integer, intent(in) :: contract_numbers(:)
    real(dp), allocatable, intent(out) :: stitched(:,:)
    integer, allocatable, intent(out), optional :: selected_columns(:,:)
    integer, intent(out), optional :: status
    integer :: t,j,k,nvalid
    integer,allocatable::valid(:)
    real(dp)::nanv

    nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    if(present(status)) status=nfcp_invalid_input
    if(any(contract_numbers<1)) then
      allocate(stitched(0,0)); if(present(selected_columns)) allocate(selected_columns(0,0)); return
    end if
    allocate(stitched(size(futures,1),size(contract_numbers)),valid(size(futures,2)))
    stitched=nanv
    if(present(selected_columns)) then
      allocate(selected_columns(size(futures,1),size(contract_numbers)))
      selected_columns=0
    end if
    do t=1,size(futures,1)
      nvalid=0
      do j=1,size(futures,2)
        if(ieee_is_finite(futures(t,j))) then
          nvalid=nvalid+1; valid(nvalid)=j
        end if
      end do
      do k=1,size(contract_numbers)
        if(contract_numbers(k)<=nvalid) then
          j=valid(contract_numbers(k))
          stitched(t,k)=futures(t,j)
          if(present(selected_columns)) selected_columns(t,k)=j
        end if
      end do
    end do
    if(present(status)) status=nfcp_ok
  end subroutine stitch_contract_numbers

  subroutine stitch_by_maturity(futures,maturity_matrix,target_maturities,rollover_frequency, &
                                stitched,stitched_maturities,selected_columns,status)
    real(dp),intent(in)::futures(:,:),maturity_matrix(:,:),target_maturities(:),rollover_frequency
    real(dp),allocatable,intent(out)::stitched(:,:),stitched_maturities(:,:)
    integer,allocatable,intent(out),optional::selected_columns(:,:)
    integer,intent(out),optional::status
    integer::nt,nc,t,k,j,best_j,segment_start,next_roll
    real(dp)::best,dist,threshold,nanv

    nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    if(present(status)) status=nfcp_invalid_input
    nt=size(futures,1);nc=size(futures,2)
    if(any(shape(maturity_matrix)/=shape(futures)) .or. rollover_frequency<=0.0_dp .or. &
       any(target_maturities<0.0_dp)) then
      allocate(stitched(0,0),stitched_maturities(0,0));
      if(present(selected_columns)) allocate(selected_columns(0,0)); return
    end if
    allocate(stitched(nt,size(target_maturities)),stitched_maturities(nt,size(target_maturities)))
    stitched=nanv;stitched_maturities=nanv
    if(present(selected_columns)) then
      allocate(selected_columns(nt,size(target_maturities)));selected_columns=0
    end if

    segment_start=1
    do while(segment_start<=nt)
      do k=1,size(target_maturities)
        best=huge(1.0_dp);best_j=0
        do j=1,nc
          if(ieee_is_finite(futures(segment_start,j)) .and. ieee_is_finite(maturity_matrix(segment_start,j))) then
            dist=abs(maturity_matrix(segment_start,j)-target_maturities(k))
            if(dist<best) then;best=dist;best_j=j;end if
          end if
        end do
        if(best_j>0) then
          threshold=maturity_matrix(segment_start,best_j)-rollover_frequency
          next_roll=nt+1
          do t=segment_start+1,nt
            if(ieee_is_finite(maturity_matrix(t,best_j))) then
              if(maturity_matrix(t,best_j)<=threshold) then;next_roll=t;exit;end if
            end if
          end do
          do t=segment_start,min(nt,next_roll-1)
            if(ieee_is_finite(futures(t,best_j))) then
              stitched(t,k)=futures(t,best_j)
              stitched_maturities(t,k)=maturity_matrix(t,best_j)
              if(present(selected_columns)) selected_columns(t,k)=best_j
            end if
          end do
        else
          next_roll=segment_start+1
        end if
      end do
      ! Use the first target contract to define the common rollover boundary.
      best=huge(1.0_dp);best_j=0
      do j=1,nc
        if(ieee_is_finite(futures(segment_start,j)) .and. ieee_is_finite(maturity_matrix(segment_start,j))) then
          dist=abs(maturity_matrix(segment_start,j)-target_maturities(1))
          if(dist<best) then;best=dist;best_j=j;end if
        end if
      end do
      if(best_j==0) then
        segment_start=segment_start+1
      else
        threshold=maturity_matrix(segment_start,best_j)-rollover_frequency
        next_roll=nt+1
        do t=segment_start+1,nt
          if(ieee_is_finite(maturity_matrix(t,best_j)) .and. maturity_matrix(t,best_j)<=threshold) then
            next_roll=t;exit
          end if
        end do
        segment_start=max(segment_start+1,next_roll)
      end if
    end do
    if(present(status)) status=nfcp_ok
  end subroutine stitch_by_maturity

end module nfcp_stitch
