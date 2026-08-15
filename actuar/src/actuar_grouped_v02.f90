module actuar_grouped_v02
    use actuar_kinds, only: dp
    implicit none
    private
    public :: empirical_moments, grouped_moments, grouped_mean, grouped_variance
    public :: grouped_quantile, ogive_eval, elev_individual, elev_grouped
contains
    pure function empirical_moments(x,orders) result(mom)
        real(dp),intent(in)::x(:),orders(:)
        real(dp),allocatable::mom(:)
        integer::i
        allocate(mom(size(orders)));mom=0.0_dp
        if(size(x)==0) return
        do i=1,size(orders)
            if(orders(i)<0.0_dp) then
                mom(i)=huge(1.0_dp)
            else
                mom(i)=sum(x**orders(i))/real(size(x),dp)
            end if
        end do
    end function empirical_moments

    pure function grouped_moments(boundaries,counts,orders) result(mom)
        real(dp),intent(in)::boundaries(:),counts(:),orders(:)
        real(dp),allocatable::mom(:)
        real(dp)::f,den
        integer::i,j
        allocate(mom(size(orders)));mom=0.0_dp
        if(size(boundaries)/=size(counts)+1) return
        den=sum(counts);if(den<=0.0_dp) return
        do j=1,size(orders)
            if(orders(j)<0.0_dp) then;mom(j)=huge(1.0_dp);cycle;end if
            do i=1,size(counts)
                if(abs(boundaries(i+1)-boundaries(i))<=tiny(1.0_dp)) cycle
                f=(boundaries(i+1)**(orders(j)+1.0_dp)-boundaries(i)**(orders(j)+1.0_dp)) &
                  /((orders(j)+1.0_dp)*(boundaries(i+1)-boundaries(i)))
                mom(j)=mom(j)+counts(i)*f
            end do
            mom(j)=mom(j)/den
        end do
    end function grouped_moments

    pure real(dp) function grouped_mean(boundaries,counts) result(mu)
        real(dp),intent(in)::boundaries(:),counts(:)
        real(dp)::den
        integer::i
        mu=0.0_dp;den=sum(counts)
        if(size(boundaries)/=size(counts)+1 .or. den<=0.0_dp) return
        do i=1,size(counts)
            mu=mu+counts(i)*0.5_dp*(boundaries(i)+boundaries(i+1))
        end do
        mu=mu/den
    end function grouped_mean

    pure real(dp) function grouped_variance(boundaries,counts) result(v)
        real(dp),intent(in)::boundaries(:),counts(:)
        real(dp)::mu,den,mid
        integer::i
        v=0.0_dp;den=sum(counts)
        if(size(boundaries)/=size(counts)+1 .or. den<=1.0_dp) return
        mu=grouped_mean(boundaries,counts)
        do i=1,size(counts)
            mid=0.5_dp*(boundaries(i)+boundaries(i+1))
            v=v+counts(i)*(mid-mu)**2
        end do
        v=v/(den-1.0_dp)
    end function grouped_variance

    pure real(dp) function ogive_eval(x,boundaries,counts) result(p)
        real(dp),intent(in)::x,boundaries(:),counts(:)
        real(dp)::total,cum,frac
        integer::i
        total=sum(counts);p=0.0_dp
        if(size(boundaries)/=size(counts)+1 .or. total<=0.0_dp) return
        if(x<=boundaries(1)) return
        if(x>=boundaries(size(boundaries))) then;p=1.0_dp;return;end if
        cum=0.0_dp
        do i=1,size(counts)
            if(x>=boundaries(i+1)) then
                cum=cum+counts(i)
            else if(x>boundaries(i)) then
                frac=(x-boundaries(i))/(boundaries(i+1)-boundaries(i))
                cum=cum+counts(i)*frac
                exit
            else
                exit
            end if
        end do
        p=cum/total
    end function ogive_eval

    pure real(dp) function grouped_quantile(prob,boundaries,counts) result(x)
        real(dp),intent(in)::prob,boundaries(:),counts(:)
        real(dp)::target,cum,next,total,frac
        integer::i
        if(size(boundaries)/=size(counts)+1 .or. sum(counts)<=0.0_dp) then;x=0.0_dp;return;end if
        if(prob<=0.0_dp) then;x=boundaries(1);return;end if
        if(prob>=1.0_dp) then;x=boundaries(size(boundaries));return;end if
        total=sum(counts);target=prob*total;cum=0.0_dp
        do i=1,size(counts)
            next=cum+counts(i)
            if(target<=next .and. counts(i)>0.0_dp) then
                frac=(target-cum)/counts(i)
                x=boundaries(i)+frac*(boundaries(i+1)-boundaries(i));return
            end if
            cum=next
        end do
        x=boundaries(size(boundaries))
    end function grouped_quantile

    pure real(dp) function elev_individual(limit,x) result(v)
        real(dp),intent(in)::limit,x(:)
        if(size(x)==0) then;v=0.0_dp;else;v=sum(min(limit,x))/real(size(x),dp);end if
    end function elev_individual

    pure real(dp) function elev_grouped(limit,boundaries,counts) result(v)
        real(dp),intent(in)::limit,boundaries(:),counts(:)
        real(dp)::l,total,lower,upper,p
        integer::i
        v=0.0_dp;total=sum(counts)
        if(size(boundaries)/=size(counts)+1 .or. total<=0.0_dp) return
        l=min(limit,boundaries(size(boundaries)))
        do i=1,size(counts)
            lower=boundaries(i);upper=boundaries(i+1)
            if(l<=lower) then
                v=v+counts(i)*l
            else if(l>=upper) then
                v=v+counts(i)*0.5_dp*(lower+upper)
            else
                p=(l-lower)/(upper-lower)
                v=v+counts(i)*(p*0.5_dp*(lower+l)+(1.0_dp-p)*l)
            end if
        end do
        v=v/total
    end function elev_grouped
end module actuar_grouped_v02
