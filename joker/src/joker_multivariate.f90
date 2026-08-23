
module joker_multivariate
  use joker_special
  implicit none
  private
  public :: dcat, rcat, llcat, ddir, rdir, lldir
  public :: dmultinom, rmultinom, llmultinom
  public :: dmultigam, rmultigam, llmultigam
contains
  pure real(dp) function dcat(x,prob,logd) result(v)
    integer,intent(in)::x
    real(dp),intent(in)::prob(:)
    logical,intent(in),optional::logd
    real(dp)::lv
    if(x<1 .or. x>size(prob))then
      lv=-huge(1.0_dp)
    else
      lv=log(max(prob(x),tiny(1.0_dp)))
    end if
    if(present(logd))then
      if(logd)then;v=lv;return;end if
    end if
    if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function

  integer function rcat(prob) result(x)
    real(dp),intent(in)::prob(:)
    real(dp)::u,c
    integer::i
    call random_number(u);c=0
    do i=1,size(prob)
      c=c+prob(i)
      if(u<=c)then;x=i;return;end if
    end do
    x=size(prob)
  end function

  pure real(dp) function llcat(x,prob) result(v)
    integer,intent(in)::x(:)
    real(dp),intent(in)::prob(:)
    integer::i
    v=0
    do i=1,size(x);v=v+dcat(x(i),prob,.true.);end do
  end function

  pure real(dp) function ddir(x,alpha,logd) result(v)
    real(dp),intent(in)::x(:),alpha(:)
    logical,intent(in),optional::logd
    real(dp)::lv
    integer::i
    if(size(x)/=size(alpha) .or. any(x<=0) .or. any(alpha<=0) .or. abs(sum(x)-1)>1e-10_dp)then
      lv=-huge(1.0_dp)
    else
      lv=log_gamma(sum(alpha))
      do i=1,size(alpha)
        lv=lv-log_gamma(alpha(i))+(alpha(i)-1)*log(x(i))
      end do
    end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if
    if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function

  subroutine rdir(alpha,x)
    real(dp),intent(in)::alpha(:)
    real(dp),intent(out)::x(:)
    integer::i
    do i=1,size(alpha);x(i)=rng_gamma(alpha(i),1.0_dp);end do
    x=x/sum(x)
  end subroutine

  pure real(dp) function lldir(x,alpha) result(v)
    real(dp),intent(in)::x(:,:),alpha(:)
    integer::i
    v=0
    do i=1,size(x,1);v=v+ddir(x(i,:),alpha,.true.);end do
  end function

  pure real(dp) function dmultinom(x,sizep,prob,logd) result(v)
    integer,intent(in)::x(:),sizep
    real(dp),intent(in)::prob(:)
    logical,intent(in),optional::logd
    real(dp)::lv
    integer::i
    if(size(x)/=size(prob) .or. sum(x)/=sizep .or. any(x<0) .or. any(prob<0))then
      lv=-huge(1.0_dp)
    else
      lv=log_gamma(real(sizep+1,dp))
      do i=1,size(x)
        lv=lv-log_gamma(real(x(i)+1,dp))
        if(x(i)>0)lv=lv+x(i)*log(max(prob(i),tiny(1.0_dp)))
      end do
    end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if
    if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function

  subroutine rmultinom(sizep,prob,x)
    integer,intent(in)::sizep
    real(dp),intent(in)::prob(:)
    integer,intent(out)::x(:)
    integer::i,k
    x=0
    do i=1,sizep
      k=rcat(prob);x(k)=x(k)+1
    end do
  end subroutine

  pure real(dp) function llmultinom(x,sizep,prob) result(v)
    integer,intent(in)::x(:,:),sizep
    real(dp),intent(in)::prob(:)
    integer::i
    v=0
    do i=1,size(x,1);v=v+dmultinom(x(i,:),sizep,prob,.true.);end do
  end function

  pure real(dp) function dmultigam(x,shape,scale,logd) result(v)
    real(dp),intent(in)::x(:),shape(:),scale
    logical,intent(in),optional::logd
    real(dp)::lv,inc
    integer::i
    if(size(x)/=size(shape) .or. scale<=0 .or. any(shape<=0))then
      lv=-huge(1.0_dp)
    else
      lv=0
      inc=x(1)
      if(inc<=0)then
        lv=-huge(1.0_dp)
      else
        lv=lv+(shape(1)-1)*log(inc)-inc/scale-log_gamma(shape(1))-shape(1)*log(scale)
        do i=2,size(shape)
          inc=x(i)-x(i-1)
          if(inc<=0)then;lv=-huge(1.0_dp);exit;end if
          lv=lv+(shape(i)-1)*log(inc)-inc/scale-log_gamma(shape(i))-shape(i)*log(scale)
        end do
      end if
    end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if
    if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function

  subroutine rmultigam(shape,scale,x)
    real(dp),intent(in)::shape(:),scale
    real(dp),intent(out)::x(:)
    integer::i
    real(dp)::s
    s=0
    do i=1,size(shape)
      s=s+rng_gamma(shape(i),scale);x(i)=s
    end do
  end subroutine

  pure real(dp) function llmultigam(x,shape,scale) result(v)
    real(dp),intent(in)::x(:,:),shape(:),scale
    integer::i
    v=0
    do i=1,size(x,1);v=v+dmultigam(x(i,:),shape,scale,.true.);end do
  end function
end module joker_multivariate
