module skewlap_distribution
  use gh_math, only: dp
  implicit none
  private
  public :: dskewlap, pskewlap, qskewlap, rskewlap
contains
  pure function dskewlap(x,mu,alpha,beta) result(d)
    real(dp),intent(in)::x,mu,alpha,beta;real(dp)::d
    if(x<=mu)then;d=exp((x-mu)/alpha)/(alpha+beta);else;d=exp((mu-x)/beta)/(alpha+beta);endif
  end function
  pure function pskewlap(x,mu,alpha,beta) result(p)
    real(dp),intent(in)::x,mu,alpha,beta;real(dp)::p
    if(x<mu)then;p=alpha*exp((x-mu)/alpha)/(alpha+beta);else;p=1-beta*exp((mu-x)/beta)/(alpha+beta);endif
  end function
  pure function qskewlap(p,mu,alpha,beta) result(q)
    real(dp),intent(in)::p,mu,alpha,beta;real(dp)::q
    if(p<=0)then;q=-huge(1.0_dp);elseif(p>=1)then;q=huge(1.0_dp)
    elseif(p<alpha/(alpha+beta))then;q=mu+alpha*log(p*(alpha+beta)/alpha)
    else;q=mu-beta*log((alpha+beta)*(1-p)/beta);endif
  end function
  subroutine rskewlap(x,mu,alpha,beta)
    real(dp),intent(out)::x(:);real(dp),intent(in)::mu,alpha,beta;real(dp)::u;integer::i
    do i=1,size(x);call random_number(u);x(i)=qskewlap(u,mu,alpha,beta);enddo
  end subroutine
end module skewlap_distribution
