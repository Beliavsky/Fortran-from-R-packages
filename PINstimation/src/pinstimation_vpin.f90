! SPDX-License-Identifier: GPL-3.0-or-later
module pinstimation_vpin
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use pinstimation_kinds, only : dp
   use pinstimation_types, only : vpin_result
   use pinstimation_math, only : normal_cdf, sample_sd, logistic, logit, softplus, inv_softplus, log_sum_exp
   use pinstimation_optimization, only : optimizer_result, minimize_nelder_mead
   implicit none
   private
   public :: compute_vpin_from_buckets, build_volume_buckets, compute_vpin
   public :: ivpin_loglik, compute_ivpin_from_buckets

contains

   subroutine compute_vpin_from_buckets(buy_volume, sell_volume, sample_length, vpin, bucket_size, status)
      real(dp), intent(in) :: buy_volume(:), sell_volume(:)
      integer, intent(in) :: sample_length
      real(dp), allocatable, intent(out) :: vpin(:)
      real(dp), intent(in), optional :: bucket_size
      integer, intent(out), optional :: status
      real(dp), allocatable :: imbalance(:), cumulative(:)
      real(dp) :: vbs, window_sum
      integer :: i, n
      if (present(status)) status = 0
      n = size(buy_volume)
      allocate(vpin(n)); vpin = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(sell_volume) /= n .or. sample_length < 1 .or. n < sample_length) then
         if (present(status)) status = 1
         return
      end if
      if (present(bucket_size)) then
         vbs = bucket_size
      else
         vbs = sum(buy_volume + sell_volume)/real(n,dp)
      end if
      if (vbs <= 0.0_dp) then
         if (present(status)) status = 2
         return
      end if
      allocate(imbalance(n), cumulative(n))
      imbalance = abs(buy_volume - sell_volume)
      cumulative(1) = imbalance(1)
      do i = 2, n
         cumulative(i) = cumulative(i-1) + imbalance(i)
      end do
      do i = sample_length, n
         if (i == sample_length) then
            window_sum = cumulative(i)
         else
            window_sum = cumulative(i) - cumulative(i-sample_length)
         end if
         vpin(i) = window_sum/(real(sample_length,dp)*vbs)
      end do
   end subroutine compute_vpin_from_buckets

   subroutine build_volume_buckets(price_change, volume, duration, trading_days, buckets_per_day, &
         buy_volume, sell_volume, bucket_duration, bucket_size, status)
      real(dp), intent(in) :: price_change(:), volume(:), duration(:)
      integer, intent(in) :: trading_days, buckets_per_day
      real(dp), allocatable, intent(out) :: buy_volume(:), sell_volume(:), bucket_duration(:)
      real(dp), intent(out) :: bucket_size
      integer, intent(out), optional :: status
      real(dp), allocatable :: btmp(:), stmp(:), dtmp(:)
      real(dp) :: sdp, fill, remaining, take, fraction, buy_probability
      integer :: i, k, max_buckets, ncomplete
      if (present(status)) status = 0
      if (size(volume) /= size(price_change) .or. size(duration) /= size(volume) .or. &
          trading_days < 1 .or. buckets_per_day < 1 .or. any(volume < 0.0_dp) .or. any(duration < 0.0_dp)) then
         allocate(buy_volume(0),sell_volume(0),bucket_duration(0)); bucket_size=0.0_dp
         if (present(status)) status=1
         return
      end if
      bucket_size = sum(volume)/(real(trading_days*buckets_per_day,dp))
      if (bucket_size <= 0.0_dp) then
         allocate(buy_volume(0),sell_volume(0),bucket_duration(0)); if(present(status)) status=2; return
      end if
      max_buckets = int(ceiling(sum(volume)/bucket_size)) + 1
      allocate(btmp(max_buckets),stmp(max_buckets),dtmp(max_buckets))
      btmp=0.0_dp;stmp=0.0_dp;dtmp=0.0_dp
      sdp=sample_sd(price_change); if(sdp<=tiny(1.0_dp)) sdp=1.0_dp
      k=1;fill=0.0_dp
      do i=1,size(volume)
         remaining=volume(i)
         buy_probability=normal_cdf(price_change(i)/sdp)
         do while(remaining>1.0e-12_dp)
            take=min(remaining,bucket_size-fill)
            if(volume(i)>0.0_dp) then;fraction=take/volume(i);else;fraction=0.0_dp;end if
            btmp(k)=btmp(k)+take*buy_probability
            stmp(k)=stmp(k)+take*(1.0_dp-buy_probability)
            dtmp(k)=dtmp(k)+duration(i)*fraction
            fill=fill+take;remaining=remaining-take
            if(fill>=bucket_size*(1.0_dp-1.0e-10_dp)) then;k=k+1;fill=0.0_dp;end if
         end do
      end do
      ncomplete=k-1
      allocate(buy_volume(ncomplete),sell_volume(ncomplete),bucket_duration(ncomplete))
      buy_volume=btmp(:ncomplete);sell_volume=stmp(:ncomplete);bucket_duration=dtmp(:ncomplete)
   end subroutine build_volume_buckets

   subroutine compute_vpin(price_change, volume, duration, trading_days, buckets_per_day, sample_length, result, &
         improved, status)
      real(dp), intent(in) :: price_change(:), volume(:), duration(:)
      integer, intent(in) :: trading_days, buckets_per_day, sample_length
      type(vpin_result), intent(out) :: result
      logical, intent(in), optional :: improved
      integer, intent(out), optional :: status
      logical :: do_ivpin
      integer :: st
      do_ivpin=.false.;if(present(improved))do_ivpin=improved
      call build_volume_buckets(price_change,volume,duration,trading_days,buckets_per_day,result%buy_volume,&
         result%sell_volume,result%duration,result%volume_bucket_size,st)
      if(st/=0) then;if(present(status))status=st;return;end if
      call compute_vpin_from_buckets(result%buy_volume,result%sell_volume,sample_length,result%vpin,&
         result%volume_bucket_size,st)
      if(do_ivpin) call compute_ivpin_from_buckets(result%buy_volume,result%sell_volume,result%duration,sample_length,result%ivpin,st)
      if(present(status))status=st
   end subroutine compute_vpin

   real(dp) function ivpin_loglik(buy_volume,sell_volume,duration,parameters) result(value)
      real(dp),intent(in)::buy_volume(:),sell_volume(:),duration(:),parameters(5)
      real(dp)::a,d,mu,eb,es,e(3)
      integer::i
      a=parameters(1);d=parameters(2);mu=parameters(3);eb=parameters(4);es=parameters(5)
      if(size(sell_volume)/=size(buy_volume).or.size(duration)/=size(buy_volume).or.a<=0.0_dp.or.a>=1.0_dp.or.&
         d<=0.0_dp.or.d>=1.0_dp.or.min(mu,eb,es)<=0.0_dp) then;value=-huge(1.0_dp);return;end if
      value=0.0_dp
      do i=1,size(buy_volume)
         e(1)=log(a*d)+buy_volume(i)*log(eb)+sell_volume(i)*log(es+mu)-(eb+es+mu)*duration(i)
         e(2)=log(a*(1.0_dp-d))+buy_volume(i)*log(eb+mu)+sell_volume(i)*log(es)-(eb+es+mu)*duration(i)
         e(3)=log(1.0_dp-a)+buy_volume(i)*log(eb)+sell_volume(i)*log(es)-(eb+es)*duration(i)
         value=value+log_sum_exp(e)
      end do
   end function ivpin_loglik

   subroutine compute_ivpin_from_buckets(buy_volume,sell_volume,duration,sample_length,ivpin,status)
      real(dp),intent(in)::buy_volume(:),sell_volume(:),duration(:)
      integer,intent(in)::sample_length
      real(dp),allocatable,intent(out)::ivpin(:)
      integer,intent(out),optional::status
      type(optimizer_result)::opt
      real(dp)::u0(5),p(5),iar,uar,total_rate,informed_rate
      integer::n,k,first
      if(present(status))status=0
      n=size(buy_volume);allocate(ivpin(n));ivpin=ieee_value(0.0_dp,ieee_quiet_nan)
      if(size(sell_volume)/=n.or.size(duration)/=n.or.sample_length<1.or.n<sample_length.or.any(duration<=0.0_dp)) then
         if(present(status))status=1;return
      end if
      do k=sample_length,n
         first=k-sample_length+1
         total_rate=sum((buy_volume(first:k)+sell_volume(first:k))/duration(first:k))/real(sample_length,dp)
         informed_rate=sum(abs(buy_volume(first:k)-sell_volume(first:k))/duration(first:k))/real(sample_length,dp)
         iar=max(informed_rate,1.0e-8_dp);uar=max(total_rate-informed_rate,1.0e-8_dp)
         if(k==sample_length) then
            p=[0.5_dp,0.5_dp,2.0_dp*iar,max(uar/2.0_dp,1.0e-8_dp),max(uar/2.0_dp,1.0e-8_dp)]
            call ivpin_pack(p,u0)
         end if
         call minimize_nelder_mead(objective,u0,opt,1.0e-7_dp,1200)
         call ivpin_unpack(opt%parameters,p)
         ivpin(k)=p(1)*p(3)/(p(4)+p(5)+p(3))
         u0=opt%parameters
      end do
      if(present(status))status=0
   contains
      real(dp) function objective(u) result(v)
         real(dp),intent(in)::u(:)
         real(dp)::pp(5)
         call ivpin_unpack(u,pp)
         v=-ivpin_loglik(buy_volume(first:k),sell_volume(first:k),duration(first:k),pp)
         if(.not.ieee_is_finite(v))v=huge(1.0_dp)/100.0_dp
      end function objective
   end subroutine compute_ivpin_from_buckets

   subroutine ivpin_pack(p,u)
      real(dp),intent(in)::p(5)
      real(dp),intent(out)::u(5)
      u=[logit(p(1)),logit(p(2)),inv_softplus(p(3)),inv_softplus(p(4)),inv_softplus(p(5))]
   end subroutine ivpin_pack

   subroutine ivpin_unpack(u,p)
      real(dp),intent(in)::u(:)
      real(dp),intent(out)::p(5)
      p=[logistic(u(1)),logistic(u(2)),softplus(u(3))+1.0e-12_dp,softplus(u(4))+1.0e-12_dp,softplus(u(5))+1.0e-12_dp]
   end subroutine ivpin_unpack

end module pinstimation_vpin
