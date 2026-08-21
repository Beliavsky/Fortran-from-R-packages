module boot_timeseries
    use boot_kinds, only : dp
    implicit none
    private
    public :: fixed_block_indices, geometric_block_indices
contains
    subroutine fixed_block_indices(n,n_sim,r,l,endcorr,indices)
        integer,intent(in)::n,n_sim,r,l
        logical,intent(in)::endcorr
        integer,intent(out)::indices(r,n_sim)
        integer::rep,pos,start,len,endpt,j
        real(dp)::u
        if(l<=0)error stop "fixed_block_indices: l <= 0"
        endpt=merge(n,n-l+1,endcorr)
        do rep=1,r
            pos=1
            do while(pos<=n_sim)
                call random_number(u)
                start=1+int(u*real(endpt,dp))
                if(start>endpt)start=endpt
                len=min(l,n_sim-pos+1)
                do j=0,len-1
                    indices(rep,pos+j)=1+mod(start+j-1,n)
                end do
                pos=pos+len
            end do
        end do
    end subroutine fixed_block_indices

    subroutine geometric_block_indices(n,n_sim,r,mean_l,endcorr,indices)
        integer,intent(in)::n,n_sim,r
        real(dp),intent(in)::mean_l
        logical,intent(in)::endcorr
        integer,intent(out)::indices(r,n_sim)
        integer::rep,pos,start,len,endpt,j
        real(dp)::u,p
        if(mean_l<=1.0_dp)error stop "geometric_block_indices: mean_l must exceed 1"
        p=1.0_dp/mean_l
        endpt=merge(n,max(1,n-int(mean_l)+1),endcorr)
        do rep=1,r
            pos=1
            do while(pos<=n_sim)
                call random_number(u)
                start=1+int(u*real(endpt,dp))
                if(start>endpt)start=endpt
                call random_number(u)
                u=max(u,tiny(1.0_dp))
                len=1+int(log(u)/log(1.0_dp-p))
                len=min(len,n_sim-pos+1)
                do j=0,len-1
                indices(rep,pos+j)=1+mod(start+j-1,n)
                end do
                pos=pos+len
            end do
        end do
    end subroutine geometric_block_indices
end module boot_timeseries
