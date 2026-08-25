module greybox_parity
    use greybox_kinds, only: dp, pi
    use greybox_special, only: normal_pdf, normal_cdf, normal_rng, student_t_quantile, beta_inc, gamma_p
    use greybox_regression, only: alm_model, alm_fit
    use greybox_dynamic, only: alm_occurrence_model, alm_fit_occurrence
    implicit none
    private
    public :: rmcb_result, rmcb_test
    public :: dsrboot_result, dsr_bootstrap
    public :: aid_result, aidcat_result, aid_fit, aid_cat

    type :: rmcb_result
        real(dp), allocatable :: mean(:)
        real(dp), allocatable :: interval(:,:)
        logical, allocatable :: groups(:,:)
        real(dp) :: p_value=1.0_dp
        real(dp) :: level=0.95_dp
        integer :: selected=1
    end type rmcb_result

    type :: dsrboot_result
        real(dp), allocatable :: boot(:,:)
        real(dp), allocatable :: smooth(:)
        character(len=16) :: bootstrap_type='multiplicative'
        character(len=16) :: kind='nonparametric'
        real(dp) :: sd=0.0_dp
    end type dsrboot_result

    type :: aid_result
        character(len=40) :: name='regular fractional'
        character(len=16) :: type1='fractional'
        character(len=16) :: type2='regular'
        character(len=16) :: type2a=''
        logical :: new_product=.false.
        logical :: obsolete=.false.
        integer, allocatable :: stockout_start(:)
        integer, allocatable :: stockout_end(:)
        real(dp), allocatable :: candidate_ic(:)
    end type aid_result

    type :: aidcat_result
        integer, allocatable :: category(:)
        integer :: types(2,3)=0
        integer :: anomalies(3)=0
    end type aidcat_result

contains

    subroutine rmcb_test(data,level,result,distribution)
        real(dp), intent(in) :: data(:,:)
        real(dp), intent(in) :: level
        type(rmcb_result), intent(out) :: result
        character(len=*), intent(in), optional :: distribution
        real(dp), allocatable :: ranks(:,:), flat(:), means(:), design(:,:), ones(:,:)
        real(dp) :: qstat, mse, sse, sst, r2, fval, df1, df2, xbeta, diff
        integer :: n,k,i,j,row
        type(alm_model) :: fit_full, fit_null
        character(len=16) :: dist
        n=size(data,1); k=size(data,2)
        if(n<2 .or. k<2) error stop 'rmcb_test: data must be at least 2 by 2'
        dist='tukey'; if(present(distribution))dist=lower_ascii(trim(distribution))
        allocate(ranks(n,k),flat(n*k),means(k))
        if(dist=='tukey'.or.dist=='dnorm')then
            call rank_rows(data,ranks)
            means=sum(ranks,dim=1)/real(n,dp)
        else
            ranks=data
            means=sum(data,dim=1)/real(n,dp)
        end if
        allocate(result%mean(k),result%interval(k,2),result%groups(k,k))
        result%mean=means;result%level=level
        if(dist=='tukey')then
            qstat=studentized_range_quantile(level,k)*sqrt(real(k*(k+1),dp)/(12.0_dp*real(n,dp)))/2.0_dp
            fval=friedman_stat(ranks)
            result%p_value=1.0_dp-gamma_p(0.5_dp*real(k-1,dp),0.5_dp*fval)
        else if(dist=='dnorm')then
            do j=1,k
                flat((j-1)*n+1:j*n)=ranks(:,j)
            end do
            sse=0.0_dp
            do j=1,k
                sse=sse+sum((ranks(:,j)-means(j))**2)
            end do
            mse=sse/real(n*k-k,dp)
            qstat=student_t_quantile((1.0_dp+level)/2.0_dp,real(n*k-k,dp))*sqrt(mse/real(n,dp))
            sst=sum((flat-sum(flat)/real(n*k,dp))**2)
            r2=max(0.0_dp,min(1.0_dp,1.0_dp-sse/max(sst,tiny(1.0_dp))))
            df1=real(k-1,dp);df2=real(n*k-k,dp)
            fval=(r2/df1)/max((1.0_dp-r2)/df2,tiny(1.0_dp))
            xbeta=(df1*fval)/(df1*fval+df2)
            result%p_value=1.0_dp-beta_inc(0.5_dp*df1,0.5_dp*df2,xbeta)
        else
            allocate(design(n*k,k),ones(n*k,1))
            design=0.0_dp;ones=1.0_dp
            do j=1,k
                do i=1,n
                    row=(j-1)*n+i
                    flat(row)=data(i,j)
                    design(row,1)=1.0_dp
                    if(j>1)design(row,j)=1.0_dp
                end do
            end do
            call alm_fit(design,flat,dist,fit_full,max_iter=700,tol=1.0e-6_dp)
            call alm_fit(ones,flat,dist,fit_null,max_iter=700,tol=1.0e-6_dp)
            if(.not.fit_full%converged.or..not.fit_null%converged)then
                result%p_value=1.0_dp;qstat=huge(1.0_dp)/100.0_dp
            else
                means(1)=fit_full%beta(1)
                do j=2,k
                    means(j)=fit_full%beta(1)+fit_full%beta(j)
                end do
                qstat=student_t_quantile((1.0_dp+level)/2.0_dp,real(n*k-k,dp))* &
                    sqrt(max(fit_full%vcov(1,1),0.0_dp))
                result%p_value=1.0_dp-gamma_p(0.5_dp*real(k-1,dp), &
                    0.5_dp*max(0.0_dp,fit_full%loglik-fit_null%loglik))
                if(dist=='dlnorm')then
                    means=exp(means)
                end if
                result%mean=means
            end if
        end if
        do j=1,k
            if(dist=='dlnorm')then
                result%interval(j,1)=exp(log(max(means(j),tiny(1.0_dp)))-qstat)
                result%interval(j,2)=exp(log(max(means(j),tiny(1.0_dp)))+qstat)
            else
                result%interval(j,1)=means(j)-qstat
                result%interval(j,2)=means(j)+qstat
            end if
        end do
        do i=1,k
            do j=1,k
                diff=abs(means(i)-means(j))
                result%groups(i,j)=diff<=2.0_dp*qstat+1.0e-12_dp
            end do
        end do
        result%selected=minloc(means,dim=1)
    end subroutine rmcb_test

    recursive subroutine dsr_bootstrap(y,nsim,result,intermittent,multiplicative,parametric,lag,sd_value,rescale)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: nsim
        type(dsrboot_result), intent(out) :: result
        logical, intent(in), optional :: intermittent,multiplicative,parametric,rescale
        integer, intent(in), optional :: lag
        real(dp), intent(in), optional :: sd_value
        logical :: inter,mult,para,do_scale,is_integer,is_binary
        integer :: n,s,i,l,idx
        real(dp), allocatable :: transformed(:),sorted(:),smooth_sorted(:),diffs(:),work(:),noise(:),colmean(:)
        real(dp) :: sdv,u,target_sd,boot_sd,mu_col

        n=size(y); if(n<1.or.nsim<1)error stop 'dsr_bootstrap: invalid dimensions'
        inter=.false.;if(present(intermittent))inter=intermittent
        mult=.true.;if(present(multiplicative))mult=multiplicative
        para=.false.;if(present(parametric))para=parametric
        do_scale=.true.;if(present(rescale))do_scale=rescale
        l=1;if(present(lag))l=max(1,lag)
        is_integer=all(abs(y-real(nint(y),dp))<1.0e-10_dp)
        is_binary=(maxval(y)>0.0_dp .and. all((abs(y)<1.0e-12_dp).or.abs(y-maxval(y))<1.0e-12_dp))
        allocate(result%boot(n,nsim),result%smooth(n))
        if(is_binary)then
            do s=1,nsim;result%boot(:,s)=y;end do
            result%smooth=y;return
        end if
        if(inter.and.any(abs(y)<=1.0e-14_dp))then
            call intermittent_bootstrap(y,nsim,result,mult,para,l,sd_value,do_scale)
            return
        end if
        allocate(transformed(n),sorted(n),smooth_sorted(n),work(n),noise(n))
        transformed=y
        if(mult)then
            where(transformed>0.0_dp)
                transformed=log(transformed)
            elsewhere
                transformed=0.0_dp
            end where
            result%bootstrap_type='multiplicative'
        else
            result%bootstrap_type='additive'
        end if
        if(para)result%kind='parametric'
        call sort_vector(transformed,sorted)
        call smooth_vector(sorted,smooth_sorted)
        call unsort_smooth(transformed,sorted,smooth_sorted,result%smooth)
        if(.not.para)then
            if(n>1)then
                allocate(diffs(n-1));diffs=abs(smooth_sorted(2:)-smooth_sorted(:n-1));call sort_inplace(diffs)
            else
                allocate(diffs(1));diffs=0.0_dp
            end if
            sdv=0.0_dp
        else
            if(present(sd_value))then
                sdv=sd_value
            else if(n>l)then
                sdv=sum(abs(transformed(l+1:)-transformed(:n-l)))/real(n-l,dp)
            else
                sdv=1.0_dp
            end if
            result%sd=sdv
        end if
        do s=1,nsim
            do i=1,n
                if(para)then
                    noise(i)=normal_rng(0.0_dp,sdv)
                else
                    call random_number(u);idx=1+int(u*real(size(diffs),dp));idx=min(size(diffs),idx)
                    call random_number(u);noise(i)=merge(diffs(idx),-diffs(idx),u>=0.5_dp)
                end if
            end do
            work=smooth_sorted+noise
            call sort_inplace(work)
            call map_sorted_back(transformed,work,result%boot(:,s))
            result%boot(:,s)=transformed+result%boot(:,s)-sum(result%boot(:,s)-transformed)/real(n,dp)
        end do
        if(do_scale.and.n>1.and.nsim>1)then
            allocate(colmean(nsim))
            do s=1,nsim;colmean(s)=sum(result%boot(:,s))/real(n,dp);end do
            target_sd=sample_sd(y)/sqrt(real(n,dp));boot_sd=sample_sd(colmean)
            if(boot_sd>tiny(1.0_dp))then
                do s=1,nsim
                    mu_col=sum(result%boot(:,s))/real(n,dp)
                    result%boot(:,s)=result%smooth+(target_sd/boot_sd)*(result%boot(:,s)-result%smooth)
                    if(mu_col>huge(1.0_dp)) result%boot(:,s)=result%boot(:,s)
                end do
            end if
        end if
        if(mult)then
            result%boot=exp(result%boot);result%smooth=exp(result%smooth)
        end if
        do i=1,n
            if(abs(y(i))<=1.0e-14_dp)result%boot(i,:)=0.0_dp
        end do
        if(is_integer)result%boot=ceiling(result%boot)
    end subroutine dsr_bootstrap

    subroutine aid_fit(y,result,criterion)
        real(dp), intent(in) :: y(:)
        type(aid_result), intent(out) :: result
        character(len=*), intent(in), optional :: criterion
        real(dp), allocatable :: x(:,:), smooth(:), occ(:)
        type(alm_model) :: frac,countmod,occmod
        type(alm_occurrence_model) :: lump_frac,lump_count
        logical :: has_zero,is_integer,is_binary
        real(dp) :: ic(4)
        integer :: n,best
        n=size(y);if(n<4)error stop 'aid_fit: at least four observations are required'
        allocate(smooth(n),x(n,2),occ(n));call smooth_vector(y,smooth)
        x(:,1)=1.0_dp;x(:,2)=smooth
        has_zero=any(abs(y)<=1.0e-14_dp)
        is_integer=all(abs(y-real(nint(y),dp))<1.0e-10_dp)
        is_binary=(maxval(y)>0.0_dp .and. all((abs(y)<1.0e-12_dp).or.abs(y-maxval(y))<1.0e-12_dp))
        call identify_stockouts(y,result)
        ic=huge(1.0_dp)
        if(has_zero)then
            where(y>0.0_dp);occ=1.0_dp;elsewhere;occ=0.0_dp;end where
            call alm_fit(x,occ,'plogis',occmod,max_iter=500,tol=1.0e-6_dp)
            if(is_binary)then
                result%name='smooth intermittent count'
                result%type1='count'
                result%type2='intermittent'
                result%type2a='smooth'
                allocate(result%candidate_ic(1));result%candidate_ic=occmod%aicc;return
            end if
            call alm_fit(x,y,'drectnorm',frac,max_iter=500,tol=1.0e-6_dp);ic(1)=frac%aicc
            call alm_fit_occurrence(x,y,'dnorm','plogis',lump_frac,max_iter=500,tol=1.0e-6_dp);ic(2)=aicc_occ(lump_frac)
            if(is_integer)then
                call alm_fit(x,y,'dnbinom',countmod,max_iter=600,tol=1.0e-6_dp);ic(3)=countmod%aicc
                call alm_fit_occurrence(x,y,'dnbinom','plogis',lump_count, &
                    max_iter=600,tol=1.0e-6_dp)
                ic(4)=aicc_occ(lump_count)
            end if
            best=minloc(ic,dim=1)
            select case(best)
            case(1);result%name='smooth intermittent fractional';result%type1='fractional';result%type2a='smooth'
            case(2);result%name='lumpy intermittent fractional';result%type1='fractional';result%type2a='lumpy'
            case(3);result%name='smooth intermittent count';result%type1='count';result%type2a='smooth'
            case(4);result%name='lumpy intermittent count';result%type1='count';result%type2a='lumpy'
            end select
            result%type2='intermittent';allocate(result%candidate_ic(4));result%candidate_ic=ic
        else
            call alm_fit(x,y,'dnorm',frac,max_iter=500,tol=1.0e-6_dp);ic(1)=frac%aicc
            if(is_integer)then
                call alm_fit(x,y,'dnbinom',countmod,max_iter=600,tol=1.0e-6_dp)
                ic(2)=countmod%aicc
            end if
            if(ic(2)<ic(1))then
                result%name='regular count'
                result%type1='count'
            else
                result%name='regular fractional'
                result%type1='fractional'
            end if
            result%type2='regular';result%type2a='';allocate(result%candidate_ic(2));result%candidate_ic=ic(1:2)
        end if
        if(present(criterion))then
            if(len_trim(criterion)<0)result%name=result%name
        end if
    end subroutine aid_fit

    subroutine aid_cat(data,result)
        real(dp), intent(in) :: data(:,:)
        type(aidcat_result), intent(out) :: result
        type(aid_result) :: one
        integer :: j,c
        allocate(result%category(size(data,2)))
        do j=1,size(data,2)
            call aid_fit(data(:,j),one)
            c=category_code(one%name);result%category(j)=c
            select case(c)
            case(1);result%types(1,1)=result%types(1,1)+1
            case(2);result%types(1,2)=result%types(1,2)+1
            case(3);result%types(1,3)=result%types(1,3)+1
            case(4);result%types(2,1)=result%types(2,1)+1
            case(5);result%types(2,2)=result%types(2,2)+1
            case(6);result%types(2,3)=result%types(2,3)+1
            end select
            if(one%new_product)result%anomalies(1)=result%anomalies(1)+1
            result%anomalies(2)=result%anomalies(2)+size(one%stockout_start)
            if(one%obsolete)result%anomalies(3)=result%anomalies(3)+1
        end do
    end subroutine aid_cat

    subroutine rank_rows(data,ranks)
        real(dp),intent(in)::data(:,:);real(dp),intent(out)::ranks(:,:)
        integer::i,j,k,m;real(dp)::less,equal
        do i=1,size(data,1)
            do j=1,size(data,2)
                less=0.0_dp;equal=0.0_dp
                do k=1,size(data,2)
                    if(data(i,k)<data(i,j))less=less+1.0_dp
                    if(abs(data(i,k)-data(i,j))<=1.0e-12_dp)equal=equal+1.0_dp
                end do
                ranks(i,j)=less+0.5_dp*(equal+1.0_dp)
            end do
        end do
        m=0;if(m<0)ranks=ranks
    end subroutine rank_rows

    pure real(dp) function friedman_stat(ranks) result(q)
        real(dp),intent(in)::ranks(:,:);real(dp)::rj(size(ranks,2));integer::n,k
        n=size(ranks,1);k=size(ranks,2);rj=sum(ranks,dim=1)
        q=12.0_dp/(real(n*k*(k+1),dp))*sum(rj*rj)-3.0_dp*real(n*(k+1),dp)
    end function friedman_stat

    real(dp) function studentized_range_quantile(p,k) result(q)
        real(dp),intent(in)::p;integer,intent(in)::k
        real(dp)::lo,hi,mid;integer::i
        lo=0.0_dp;hi=10.0_dp
        do i=1,70
            mid=0.5_dp*(lo+hi)
            if(studentized_range_cdf(mid,k)<p)then;lo=mid;else;hi=mid;end if
        end do
        q=0.5_dp*(lo+hi)
    end function studentized_range_quantile

    real(dp) function studentized_range_cdf(q,k) result(v)
        real(dp),intent(in)::q;integer,intent(in)::k
        integer,parameter::nstep=2400
        real(dp)::a,b,h,x,fx,s;integer::i
        a=-8.0_dp;b=8.0_dp;h=(b-a)/real(nstep,dp);s=0.0_dp
        do i=0,nstep
            x=a+h*real(i,dp)
            fx=real(k,dp)*normal_pdf(x,0.0_dp,1.0_dp)* &
                max(0.0_dp,normal_cdf(x+q,0.0_dp,1.0_dp)-normal_cdf(x,0.0_dp,1.0_dp))**(k-1)
            if(i==0.or.i==nstep)then;s=s+fx;else if(mod(i,2)==0)then;s=s+2.0_dp*fx;else;s=s+4.0_dp*fx;end if
        end do
        v=min(1.0_dp,max(0.0_dp,s*h/3.0_dp))
    end function studentized_range_cdf

    subroutine smooth_vector(y,s)
        real(dp),intent(in)::y(:);real(dp),intent(out)::s(:);integer::i,j1,j2,w
        w=max(1,min(5,size(y)/5))
        do i=1,size(y);j1=max(1,i-w);j2=min(size(y),i+w);s(i)=sum(y(j1:j2))/real(j2-j1+1,dp);end do
    end subroutine smooth_vector

    subroutine sort_vector(x,y)
        real(dp),intent(in)::x(:);real(dp),intent(out)::y(:);y=x;call sort_inplace(y)
    end subroutine sort_vector
    subroutine sort_inplace(a)
        real(dp),intent(inout)::a(:);real(dp)::t;integer::i,j
        do i=2,size(a);t=a(i);j=i-1;do while(j>=1);if(a(j)<=t)exit;a(j+1)=a(j);j=j-1;end do;a(j+1)=t;end do
    end subroutine sort_inplace
    subroutine map_sorted_back(original,sorted_values,out)
        real(dp),intent(in)::original(:),sorted_values(:);real(dp),intent(out)::out(:)
        integer::i,j,rank;logical::used(size(original))
        used=.false.
        do i=1,size(original)
            rank=1
            do j=1,size(original)
                if(original(j)<original(i))rank=rank+1
                if(abs(original(j)-original(i))<=1.0e-14_dp.and.j<i)rank=rank+1
            end do
            out(i)=sorted_values(rank)
        end do
        if(any(used))out=out
    end subroutine map_sorted_back
    subroutine unsort_smooth(original,sorted,ss,out)
        real(dp),intent(in)::original(:),sorted(:),ss(:);real(dp),intent(out)::out(:)
        call map_sorted_back(original,ss,out)
        if(sum(sorted)>huge(1.0_dp))out=out
    end subroutine unsort_smooth

    pure real(dp) function sample_sd(x) result(s)
        real(dp),intent(in)::x(:);real(dp)::m
        if(size(x)<2)then;s=0.0_dp;return;end if
        m=sum(x)/real(size(x),dp);s=sqrt(sum((x-m)**2)/real(size(x)-1,dp))
    end function sample_sd

    subroutine intermittent_bootstrap(y,nsim,result,mult,para,lag,sd_value,do_scale)
        real(dp),intent(in)::y(:);integer,intent(in)::nsim;type(dsrboot_result),intent(out)::result
        logical,intent(in)::mult,para,do_scale;integer,intent(in)::lag;real(dp),intent(in),optional::sd_value
        real(dp),allocatable::sizes(:),intervals(:);type(dsrboot_result)::sb,ib
        integer,allocatable::nz(:);integer::i,j,k,nnz,pos
        if (lag < 0) error stop 'intermittent_bootstrap: invalid lag'
        nnz=count(y>0.0_dp);allocate(sizes(nnz),intervals(nnz),nz(nnz));j=0;pos=0
        do i=1,size(y);if(y(i)>0.0_dp)then;j=j+1;nz(j)=i;sizes(j)=y(i);intervals(j)=real(i-pos,dp);pos=i;end if;end do
        call dsr_bootstrap(sizes,nsim,sb,intermittent=.false.,multiplicative=mult, &
            parametric=para,lag=1,sd_value=sd_value,rescale=do_scale)
        call dsr_bootstrap(intervals,nsim,ib,intermittent=.false.,multiplicative=.false., &
            parametric=para,lag=1,sd_value=sd_value,rescale=do_scale)
        allocate(result%boot(size(y),nsim),result%smooth(size(y)));result%boot=0.0_dp;result%smooth=y
        do k=1,nsim
            pos=0
            do j=1,nnz
                pos=pos+max(1,nint(abs(ib%boot(j,k))))
                if(pos>size(y))exit
                result%boot(pos,k)=max(0.0_dp,sb%boot(j,k))
            end do
        end do
        result%bootstrap_type=merge('multiplicative  ','additive        ',mult)
        result%kind=merge('parametric      ','nonparametric   ',para)
    end subroutine intermittent_bootstrap

    subroutine identify_stockouts(y,result)
        real(dp),intent(in)::y(:);type(aid_result),intent(inout)::result
        integer,allocatable::starts(:),ends(:);integer::i,n,nrun,c,thr,firstnz,lastnz
        real(dp)::mean_gap,p
        n=size(y);firstnz=0;lastnz=0
        do i=1,n;if(y(i)>0.0_dp)then;if(firstnz==0)firstnz=i;lastnz=i;end if;end do
        if(firstnz==0)then
            result%new_product=.true.
            result%obsolete=.true.
            allocate(result%stockout_start(0),result%stockout_end(0))
            return
        end if
        result%new_product=firstnz>1;result%obsolete=lastnz<n
        mean_gap=real(n,dp)/real(max(1,count(y>0.0_dp)),dp);p=1.0_dp/max(mean_gap,1.0_dp)
        if(p>=1.0_dp)then;thr=2;else;thr=max(2,ceiling(log(0.01_dp)/log(1.0_dp-p)-1.0_dp));end if
        nrun=0;i=firstnz+1
        do while(i<lastnz)
            if(abs(y(i))<=1.0e-14_dp)then
                c=i
                do while(i<=lastnz.and.abs(y(i))<=1.0e-14_dp)
                    i=i+1
                end do
                if(i-c>=thr)nrun=nrun+1
            else
                i=i+1
            end if
        end do
        allocate(starts(nrun),ends(nrun));nrun=0;i=firstnz+1
        do while(i<lastnz)
            if(abs(y(i))<=1.0e-14_dp)then
                c=i
                do while(i<=lastnz.and.abs(y(i))<=1.0e-14_dp)
                    i=i+1
                end do
                if(i-c>=thr)then
                    nrun=nrun+1
                    starts(nrun)=c
                    ends(nrun)=i-1
                end if
            else
                i=i+1
            end if
        end do
        call move_alloc(starts,result%stockout_start);call move_alloc(ends,result%stockout_end)
    end subroutine identify_stockouts

    pure real(dp) function aicc_occ(m) result(v)
        type(alm_occurrence_model),intent(in)::m;real(dp)::n,k
        n=real(m%nobs,dp);k=real(m%nparam,dp);v=m%aic
        if(n>k+1.0_dp)v=v+2.0_dp*k*(k+1.0_dp)/(n-k-1.0_dp)
    end function aicc_occ
    pure integer function category_code(name) result(c)
        character(len=*),intent(in)::name
        select case(trim(name))
        case('regular count');c=1
        case('smooth intermittent count');c=2
        case('lumpy intermittent count');c=3
        case('regular fractional');c=4
        case('smooth intermittent fractional');c=5
        case default;c=6
        end select
    end function category_code
    pure function lower_ascii(s) result(out)
        character(len=*),intent(in)::s;character(len=len(s))::out;integer::i,k
        out=s;do i=1,len(s);k=iachar(out(i:i));if(k>=65.and.k<=90)out(i:i)=achar(k+32);end do
    end function lower_ascii

end module greybox_parity
