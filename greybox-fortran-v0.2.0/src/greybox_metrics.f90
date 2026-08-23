module greybox_metrics
    use greybox_kinds, only: dp, pi
    implicit none
    private
    public :: me, mae, mse, mre, mis, mpe, mape, mase, rmsse, same
    public :: rmae, rrmse, rame, rmis, smse, spis, sce, smis, gmrae
    public :: pinball, hm, ham, asymmetry, extremity, cextremity, measures

contains

    pure real(dp) function me(actual,forecast) result(v)
        real(dp),intent(in)::actual(:),forecast(:)
        if(size(actual)/=size(forecast))then;v=huge(1.0_dp);return;end if
        v=sum(actual-forecast)/real(size(actual),dp)
    end function me
    pure real(dp) function mae(actual,forecast) result(v)
        real(dp),intent(in)::actual(:),forecast(:)
        if(size(actual)/=size(forecast))then;v=huge(1.0_dp);return;end if
        v=sum(abs(actual-forecast))/real(size(actual),dp)
    end function mae
    pure real(dp) function mse(actual,forecast) result(v)
        real(dp),intent(in)::actual(:),forecast(:)
        if(size(actual)/=size(forecast))then;v=huge(1.0_dp);return;end if
        v=sum((actual-forecast)**2)/real(size(actual),dp)
    end function mse
    pure complex(dp) function mre(actual,forecast) result(v)
        real(dp),intent(in)::actual(:),forecast(:)
        integer::i
        v=(0.0_dp,0.0_dp)
        if(size(actual)/=size(forecast).or.size(actual)==0)return
        do i=1,size(actual)
            v=v+sqrt(cmplx(actual(i)-forecast(i),0.0_dp,kind=dp))
        end do
        v=v/real(size(actual),dp)
    end function mre
    pure real(dp) function mis(actual,lower,upper,level) result(v)
        real(dp),intent(in)::actual(:),lower(:),upper(:),level
        real(dp)::lev,alpha
        integer::i
        if(size(actual)/=size(lower).or.size(actual)/=size(upper))then;v=huge(1.0_dp);return;end if
        lev=level;if(lev>1.0_dp)lev=lev/100.0_dp;alpha=1.0_dp-lev
        v=0.0_dp
        do i=1,size(actual)
            v=v+upper(i)-lower(i)
            if(actual(i)<lower(i))v=v+2.0_dp/alpha*(lower(i)-actual(i))
            if(actual(i)>upper(i))v=v+2.0_dp/alpha*(actual(i)-upper(i))
        end do
        v=v/real(size(actual),dp)
    end function mis
    pure real(dp) function mpe(actual,forecast) result(v)
        real(dp),intent(in)::actual(:),forecast(:)
        if(size(actual)/=size(forecast))then;v=huge(1.0_dp);return;end if
        v=sum((actual-forecast)/actual)/real(size(actual),dp)
    end function mpe
    pure real(dp) function mape(actual,forecast) result(v)
        real(dp),intent(in)::actual(:),forecast(:)
        if(size(actual)/=size(forecast))then;v=huge(1.0_dp);return;end if
        v=sum(abs((actual-forecast)/actual))/real(size(actual),dp)
    end function mape
    pure real(dp) function mase(actual,forecast,scale) result(v)
        real(dp),intent(in)::actual(:),forecast(:),scale
        v=mae(actual,forecast)/scale
    end function mase
    pure real(dp) function rmsse(actual,forecast,scale) result(v)
        real(dp),intent(in)::actual(:),forecast(:),scale
        v=sqrt(mse(actual,forecast)/scale)
    end function rmsse
    pure real(dp) function same(actual,forecast,scale) result(v)
        real(dp),intent(in)::actual(:),forecast(:),scale
        v=abs(me(actual,forecast))/scale
    end function same
    pure real(dp) function rmae(actual,forecast,benchmark) result(v)
        real(dp),intent(in)::actual(:),forecast(:),benchmark(:)
        if (same_forecast(forecast,benchmark)) then
            v = 1.0_dp
        else
            v = mae(actual,forecast)/mae(actual,benchmark)
        end if
    end function rmae
    pure real(dp) function rrmse(actual,forecast,benchmark) result(v)
        real(dp),intent(in)::actual(:),forecast(:),benchmark(:)
        if (same_forecast(forecast,benchmark)) then
            v = 1.0_dp
        else
            v = sqrt(mse(actual,forecast)/mse(actual,benchmark))
        end if
    end function rrmse
    pure real(dp) function rame(actual,forecast,benchmark) result(v)
        real(dp),intent(in)::actual(:),forecast(:),benchmark(:)
        if (same_forecast(forecast,benchmark)) then
            v = 1.0_dp
        else
            v = abs(me(actual,forecast))/abs(me(actual,benchmark))
        end if
    end function rame
    pure real(dp) function rmis(actual,lower,upper,benchmark_lower,benchmark_upper,level) result(v)
        real(dp),intent(in)::actual(:),lower(:),upper(:),benchmark_lower(:),benchmark_upper(:),level
        v=mis(actual,lower,upper,level)/mis(actual,benchmark_lower,benchmark_upper,level)
    end function rmis
    pure real(dp) function smse(actual,forecast,scale) result(v)
        real(dp),intent(in)::actual(:),forecast(:),scale
        v=mse(actual,forecast)/scale
    end function smse
    pure real(dp) function spis(actual,forecast,scale) result(v)
        real(dp),intent(in)::actual(:),forecast(:),scale
        real(dp)::stock
        integer::i
        stock=0.0_dp;v=0.0_dp
        do i=1,size(actual);stock=stock+forecast(i)-actual(i);v=v+stock;end do
        v=v/scale
    end function spis
    pure real(dp) function sce(actual,forecast,scale) result(v)
        real(dp),intent(in)::actual(:),forecast(:),scale
        v=sum(actual-forecast)/scale
    end function sce
    pure real(dp) function smis(actual,lower,upper,scale,level) result(v)
        real(dp),intent(in)::actual(:),lower(:),upper(:),scale,level
        v=mis(actual,lower,upper,level)/scale
    end function smis
    pure real(dp) function gmrae(actual,forecast,benchmark) result(v)
        real(dp),intent(in)::actual(:),forecast(:),benchmark(:)
        integer::i,n
        real(dp)::s,r
        s=0.0_dp;n=0
        do i=1,size(actual)
            if(abs(actual(i)-benchmark(i)) > epsilon(1.0_dp)*max(1.0_dp,abs(benchmark(i))))then
                r=abs((actual(i)-forecast(i))/(actual(i)-benchmark(i)))
                if(r>0.0_dp)then;s=s+log(r);n=n+1;end if
            end if
        end do
        if(n==0)then;v=0.0_dp;else;v=exp(s/real(n,dp));end if
    end function gmrae

    pure complex(dp) function hm(x,c) result(v)
        real(dp),intent(in)::x(:),c
        integer::i
        v=(0.0_dp,0.0_dp)
        do i=1,size(x);v=v+sqrt(cmplx(x(i)-c,0.0_dp,kind=dp));end do
        v=v/real(size(x),dp)
    end function hm
    pure real(dp) function ham(x,c) result(v)
        real(dp),intent(in)::x(:),c
        v=sum(sqrt(abs(x-c)))/real(size(x),dp)
    end function ham
    pure real(dp) function asymmetry(x,c) result(v)
        real(dp),intent(in)::x(:),c
        complex(dp)::z
        z=hm(x,c);v=1.0_dp-atan2(aimag(z),real(z,dp))/(pi/4.0_dp)
    end function asymmetry
    pure real(dp) function extremity(x,c) result(v)
        real(dp),intent(in)::x(:),c
        real(dp)::den,power
        den=(sum((x-c)**2)/real(size(x),dp))**0.25_dp
        power=log(0.5_dp)/log(2.0_dp*3.0_dp**(-0.75_dp))
        v=2.0_dp*(ham(x,c)/den)**power-1.0_dp
    end function extremity
    pure complex(dp) function cextremity(x,c) result(v)
        real(dp),intent(in)::x(:),c
        complex(dp)::ch
        real(dp)::den,power,vr,vi
        den=(sum((x-c)**2)/real(size(x),dp))**0.25_dp
        ch=hm(x,c)/den
        power=log(0.5_dp)/log(2.0_dp*3.0_dp**(-0.75_dp))
        vr=2.0_dp*(max(0.0_dp,2.0_dp*real(ch,dp)))**power-1.0_dp
        vi=2.0_dp*(max(0.0_dp,2.0_dp*aimag(ch)))**power-1.0_dp
        v=cmplx(vr,vi,kind=dp)
    end function cextremity

    pure real(dp) function pinball(actual,forecast,level,loss) result(v)
        real(dp),intent(in)::actual(:),forecast(:),level,loss
        real(dp)::e
        integer::i
        v=0.0_dp
        do i=1,size(actual)
            e=abs(actual(i)-forecast(i))**loss
            if(actual(i)<=forecast(i))then;v=v+(1.0_dp-level)*e;else;v=v+level*e;end if
        end do
    end function pinball

    pure function measures(holdout,forecast,actual,use_mean_benchmark) result(v)
        real(dp),intent(in)::holdout(:),forecast(:),actual(:)
        logical,intent(in),optional::use_mean_benchmark
        real(dp)::v(16),benchmark(size(holdout)),a_mean,ma,md,md2
        logical::bm
        integer::i
        bm=.false.;if(present(use_mean_benchmark))bm=use_mean_benchmark
        a_mean=sum(actual)/real(size(actual),dp)
        if(bm)then;benchmark=a_mean;else;benchmark=actual(size(actual));end if
        ma=sum(abs(actual))/real(size(actual),dp)
        md=0.0_dp;md2=0.0_dp
        if(size(actual)>1)then
            do i=2,size(actual);md=md+abs(actual(i)-actual(i-1));md2=md2+(actual(i)-actual(i-1))**2;end do
            md=md/real(size(actual)-1,dp);md2=md2/real(size(actual)-1,dp)
        end if
        v(1)=me(holdout,forecast);v(2)=mae(holdout,forecast);v(3)=mse(holdout,forecast)
        v(4)=mpe(holdout,forecast);v(5)=mape(holdout,forecast)
        v(6)=sce(holdout,forecast,ma);v(7)=mase(holdout,forecast,ma);v(8)=smse(holdout,forecast,ma*ma)
        v(9)=mase(holdout,forecast,md);v(10)=rmsse(holdout,forecast,md2);v(11)=same(holdout,forecast,md)
        v(12)=rmae(holdout,forecast,benchmark);v(13)=rrmse(holdout,forecast,benchmark);v(14)=rame(holdout,forecast,benchmark)
        v(15)=asymmetry(holdout-forecast,0.0_dp);v(16)=spis(holdout,forecast,ma)
    end function measures

    pure logical function same_forecast(x,y) result(same)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: scale
        scale = max(1.0_dp,maxval(abs(y)))
        same = maxval(abs(x-y)) <= epsilon(1.0_dp)*scale
    end function same_forecast

end module greybox_metrics
