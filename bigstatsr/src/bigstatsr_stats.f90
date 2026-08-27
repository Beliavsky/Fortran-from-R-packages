! SPDX-License-Identifier: GPL-3.0-only
module bigstatsr_stats
    use bigstatsr_kinds, only: dp
    use bigstatsr_utils, only: sort_pairs, normal_quantile, correlation
    use la_lapack_d, only: dgesv => la_dgesv
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    implicit none
    private

    type, public :: auc_boot_result
        real(dp) :: mean = 0.0_dp
        real(dp) :: lower = 0.0_dp
        real(dp) :: upper = 0.0_dp
        real(dp) :: sd = 0.0_dp
        integer :: n_invalid = 0
    end type auc_boot_result

    type, public :: pcor_result
        real(dp) :: r = 0.0_dp
        real(dp) :: lower = 0.0_dp
        real(dp) :: upper = 0.0_dp
        integer :: df = 0
        integer :: info = 0
    end type pcor_result

    public :: auc, auc_sorted, auc_sorted_weighted, auc_bootstrap, pcor

contains

    function auc(pred,target) result(ans)
        real(dp), intent(in) :: pred(:)
        integer, intent(in) :: target(:)
        real(dp) :: ans
        real(dp), allocatable :: x(:)
        integer, allocatable :: y(:)
        if (size(pred)/=size(target)) error stop 'auc: length mismatch'
        if (any(target<0) .or. any(target>1)) error stop 'auc: target must be 0/1'
        x=pred
        y=target
        call sort_pairs(x,y)
        ans=auc_sorted(x,y)
    end function auc

    function auc_sorted(x,y) result(ans)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: y(:)
        real(dp) :: ans, latest, count_control, count_equal, add, total
        integer :: i,ncase
        logical :: have_latest
        if (size(x)/=size(y)) error stop 'auc_sorted: length mismatch'
        latest=0.0_dp
        have_latest=.false.
        count_control=0.0_dp
        count_equal=0.0_dp
        total=0.0_dp
        ncase=0
        do i=1,size(y)
            if (y(i)==1) then
                ncase=ncase+1
                add=count_control
                if (have_latest) then
                    if (abs(x(i)-latest)<=0.0_dp) add=add-(count_equal+1.0_dp)/2.0_dp
                end if
                total=total+add
            else
                count_control=count_control+1.0_dp
                if (have_latest .and. abs(x(i)-latest)<=0.0_dp) then
                    count_equal=count_equal+1.0_dp
                else
                    latest=x(i)
                    have_latest=.true.
                    count_equal=0.0_dp
                end if
            end if
        end do
        if (count_control<=0.0_dp .or. ncase==0) then
            ans=ieee_value(0.0_dp, ieee_quiet_nan)
        else
            ans=total/(count_control*real(ncase,dp))
        end if
    end function auc_sorted

    function auc_sorted_weighted(x,y,w) result(ans)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: y(:),w(:)
        real(dp) :: ans,latest,count_control,count_equal,add,total,count_case
        integer :: i
        logical :: have_latest
        if (size(x)/=size(y) .or. size(w)/=size(y)) error stop 'auc_sorted_weighted: mismatch'
        latest=0.0_dp
        have_latest=.false.
        count_control=0.0_dp
        count_case=0.0_dp
        count_equal=0.0_dp
        total=0.0_dp
        do i=1,size(y)
            if (w(i)<=0) cycle
            if (y(i)==1) then
                count_case=count_case+real(w(i),dp)
                add=count_control
                if (have_latest) then
                    if (abs(x(i)-latest)<=0.0_dp) add=add-(count_equal+1.0_dp)/2.0_dp
                end if
                total=total+real(w(i),dp)*add
            else
                count_control=count_control+real(w(i),dp)
                if (have_latest .and. abs(x(i)-latest)<=0.0_dp) then
                    count_equal=count_equal+real(w(i),dp)
                else
                    latest=x(i)
                    have_latest=.true.
                    count_equal=0.0_dp
                end if
            end if
        end do
        if (count_control<=0.0_dp .or. count_case<=0.0_dp) then
            ans=ieee_value(0.0_dp, ieee_quiet_nan)
        else
            ans=total/(count_control*count_case)
        end if
    end function auc_sorted_weighted

    function auc_bootstrap(pred,target,nboot) result(res)
        real(dp), intent(in) :: pred(:)
        integer, intent(in) :: target(:),nboot
        type(auc_boot_result) :: res
        real(dp), allocatable :: x(:),rep(:),valid(:)
        integer, allocatable :: y(:),w(:)
        real(dp) :: u,mu
        integer :: b,i,k,nv
        x=pred
        y=target
        call sort_pairs(x,y)
        allocate(rep(nboot),w(size(y)))
        nv=0
        do b=1,nboot
            w=0
            do i=1,size(y)
                call random_number(u)
                k=min(size(y),1+int(u*real(size(y),dp)))
                w(k)=w(k)+1
            end do
            rep(b)=auc_sorted_weighted(x,y,w)
            if (.not. ieee_is_nan(rep(b))) nv=nv+1
        end do
        res%n_invalid=nboot-nv
        if (nv==0) return
        allocate(valid(nv))
        k=0
        do b=1,nboot
            if (.not. ieee_is_nan(rep(b))) then
                k=k+1
                valid(k)=rep(b)
            end if
        end do
        call sort_real(valid)
        mu=sum(valid)/real(nv,dp)
        res%mean=mu
        if (nv>1) res%sd=sqrt(sum((valid-mu)**2)/real(nv-1,dp))
        res%lower=quantile_sorted(valid,0.025_dp)
        res%upper=quantile_sorted(valid,0.975_dp)
    end function auc_bootstrap

    function pcor(x,y,z,alpha) result(res)
        real(dp), intent(in) :: x(:),y(:)
        real(dp), intent(in), optional :: z(:,:)
        real(dp), intent(in), optional :: alpha
        type(pcor_result) :: res
        real(dp), allocatable :: design(:,:),rx(:),ry(:)
        real(dp) :: a,zf,rad
        integer :: n,p,info1,info2
        n=size(x)
        if (size(y)/=n) then
            res%info=-1
            return
        end if
        p=1
        if (present(z)) then
            if (size(z,1)/=n) then
                res%info=-2
                return
            end if
            p=p+size(z,2)
        end if
        allocate(design(n,p))
        design(:,1)=1.0_dp
        if (present(z)) design(:,2:p)=z
        call ls_residual(design,x,rx,info1)
        call ls_residual(design,y,ry,info2)
        if (info1/=0 .or. info2/=0) then
            res%info=1
            return
        end if
        res%df=n-p
        if (res%df<3) then
            res%info=2
            return
        end if
        res%r=correlation(rx,ry)
        res%r=max(-1.0_dp+epsilon(1.0_dp),min(1.0_dp-epsilon(1.0_dp),res%r))
        a=0.05_dp
        if (present(alpha)) a=alpha
        zf=0.5_dp*(log(1.0_dp+res%r)-log(1.0_dp-res%r))
        rad=normal_quantile(1.0_dp-a/2.0_dp)/sqrt(real(res%df-2,dp))
        res%lower=tanh(zf-rad)
        res%upper=tanh(zf+rad)
    end function pcor

    subroutine ls_residual(a,y,r,info)
        real(dp), intent(in) :: a(:,:),y(:)
        real(dp), allocatable, intent(out) :: r(:)
        integer, intent(out) :: info
        real(dp), allocatable :: ata(:,:),rhs(:,:)
        integer, allocatable :: ipiv(:)
        integer :: p
        p=size(a,2)
        allocate(ata(p,p),rhs(p,1),ipiv(p),r(size(y)))
        ata=matmul(transpose(a),a)
        rhs(:,1)=matmul(transpose(a),y)
        call dgesv(p,1,ata,p,ipiv,rhs,p,info)
        if (info==0) then
            r=y-matmul(a,rhs(:,1))
        else
            r=0.0_dp
        end if
    end subroutine ls_residual

    pure real(dp) function quantile_sorted(x,p) result(q)
        real(dp), intent(in) :: x(:),p
        real(dp) :: h,g
        integer :: j,n
        n=size(x)
        if (n==1) then
            q=x(1)
            return
        end if
        h=1.0_dp+(real(n-1,dp))*p
        j=floor(h)
        g=h-real(j,dp)
        if (j>=n) then
            q=x(n)
        else
            q=(1.0_dp-g)*x(j)+g*x(j+1)
        end if
    end function quantile_sorted

    subroutine sort_real(x)
        real(dp), intent(inout) :: x(:)
        integer :: i,j
        real(dp) :: v
        do i=2,size(x)
            v=x(i)
            j=i-1
            do while (j>=1)
                if (x(j)<=v) exit
                x(j+1)=x(j)
                j=j-1
            end do
            x(j+1)=v
        end do
    end subroutine sort_real

end module bigstatsr_stats
