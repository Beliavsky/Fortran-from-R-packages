module boot_envelope
    use boot_kinds, only : dp
    use boot_statistics, only : sort_real
    implicit none
    private
    public :: confidence_envelope
contains
    subroutine confidence_envelope(mat,level_point,level_overall,point_lo,point_hi,overall_lo,overall_hi, &
                                   point_error,overall_error,k_point,k_overall)
        real(dp),intent(in)::mat(:,:),level_point,level_overall
        real(dp),intent(out)::point_lo(size(mat,2)),point_hi(size(mat,2))
        real(dp),intent(out)::overall_lo(size(mat,2)),overall_hi(size(mat,2))
        real(dp),intent(out)::point_error,overall_error
        integer,intent(out)::k_point,k_overall
        integer::r,p,j,kov,kee,kk
        real(dp),allocatable::ranks(:,:),z(:)
        real(dp)::al,eov,eee,etmp
        r=size(mat,1)
        p=size(mat,2)
        allocate(ranks(r,p),z(r))
        do j=1,p
            call column_ranks(mat(:,j),ranks(:,j))
        end do
        k_point=max(1,floor(real(r+1,dp)*(1.0_dp-level_point)/2.0_dp+1.0e-10_dp))
        point_error=empirical_error(ranks,k_point)
        kov=1
        eov=empirical_error(ranks,kov)
        kee=k_point
        eee=point_error
        al=1.0_dp-level_overall
        if(eov<=al)then
            do while(kee>kov+1 .and. eee>=al)
                if(abs(eee-eov)<1.0e-14_dp)then
                kk=(kov+kee)/2
                else
                kk=kov+nint(real(kee-kov,dp)*(al-eov)/(eee-eov))
                end if
                kk=max(kov+1,min(kee-1,kk))
                etmp=empirical_error(ranks,kk)
                if(etmp>al)then
                kee=kk
                eee=etmp
                else
                kov=kk
                eov=etmp
                end if
            end do
        end if
        k_overall=kov
        overall_error=eov
        do j=1,p
            z=mat(:,j)
            call sort_real(z)
            point_lo(j)=z(k_point)
            point_hi(j)=z(r+1-k_point)
            overall_lo(j)=z(k_overall)
            overall_hi(j)=z(r+1-k_overall)
        end do
    end subroutine confidence_envelope

    real(dp) function empirical_error(ranks,k) result(err)
        real(dp),intent(in)::ranks(:,:)
        integer,intent(in)::k
        integer::i,r
        r=size(ranks,1)
        err=0.0_dp
        do i=1,r
            if(minval(ranks(i,:))<=real(k,dp) .or. maxval(ranks(i,:))>=real(r+1-k,dp))err=err+1.0_dp
        end do
        err=err/real(r+1,dp)
    end function empirical_error

    subroutine column_ranks(x,rnk)
        real(dp),intent(in)::x(:)
        real(dp),intent(out)::rnk(size(x))
        integer::i,j,n,cntlt,cnteq
        n=size(x)
        do i=1,n
            cntlt=0
            cnteq=0
            do j=1,n
                if(x(j)<x(i))cntlt=cntlt+1
                if(abs(x(j)-x(i))<=epsilon(1.0_dp)*max(1.0_dp,abs(x(i))))cnteq=cnteq+1
            end do
            rnk(i)=real(cntlt,dp)+(real(cnteq,dp)+1.0_dp)/2.0_dp
        end do
    end subroutine column_ranks
end module boot_envelope
