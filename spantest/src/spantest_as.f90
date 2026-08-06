! SPDX-License-Identifier: GPL-3.0-only
module spantest_as
  use spantest_kinds, only : dp, pi_dp
  use spantest_types, only : as_result, span_ok, span_invalid_input, span_singular, span_numerical_failure
  use spantest_linalg, only : residualize_matrix
  use spantest_probability, only : student_t_cdf
  use spantest_random, only : rng_state, rng_seed, rng_normal
  implicit none
  private
  public :: span_as, cauchy_pvalue

contains

  pure real(dp) function cauchy_pvalue(p) result(out)
    real(dp), intent(in) :: p(:)
    if (size(p)==0) then
      out=0.0_dp
    else
      out=0.5_dp-atan(sum(tan((0.5_dp-p)*pi_dp))/real(size(p),dp))/pi_dp
      out=max(0.0_dp,min(1.0_dp,out))
    end if
  end function cauchy_pvalue

  subroutine student_subseries_pvalues(score,kexp,p,info)
    real(dp), intent(in) :: score(:,:)
    real(dp), intent(in) :: kexp
    real(dp), intent(out) :: p(:)
    integer, intent(out) :: info
    real(dp), allocatable :: fm(:,:), means(:),sd(:)
    integer, allocatable :: count(:)
    integer :: n,m,nf,i,j,g
    real(dp) :: tstat
    n=size(score,1); m=size(score,2)
    nf=int(floor(real(n,dp)**kexp))
    info=0; p=1.0_dp
    if (nf<2 .or. nf>n .or. size(p)/=m) then
      info=1; return
    end if
    allocate(fm(nf,m),means(m),sd(m),count(nf))
    fm=0.0_dp; count=0
    do i=1,n
      g=1+((i-1)*nf)/n
      fm(g,:)=fm(g,:)+score(i,:)
      count(g)=count(g)+1
    end do
    do g=1,nf
      if (count(g)<=0) then
        info=1; return
      end if
      fm(g,:)=fm(g,:)/real(count(g),dp)
    end do
    means=sum(fm,dim=1)/real(nf,dp)
    do j=1,m
      sd(j)=sqrt(sum((fm(:,j)-means(j))**2)/real(nf-1,dp))
      if (sd(j)<=sqrt(tiny(1.0_dp))) then
        if (abs(means(j))<=100.0_dp*epsilon(1.0_dp)) then
          p(j)=1.0_dp
        else
          p(j)=0.0_dp
        end if
      else
        tstat=means(j)/(sd(j)/sqrt(real(nf,dp)))
        p(j)=2.0_dp*student_t_cdf(-abs(tstat),real(nf-1,dp))
        p(j)=max(0.0_dp,min(1.0_dp,p(j)))
      end if
    end do
  end subroutine student_subseries_pvalues

  subroutine random_product_weights(n,l,seed,w)
    integer, intent(in) :: n,l,seed
    real(dp), intent(out) :: w(n)
    type(rng_state) :: rng
    integer :: i,j
    if (l<=0) then
      w=1.0_dp
      return
    end if
    call rng_seed(rng,seed)
    w=1.0_dp
    do j=1,l
      do i=1,n
        w(i)=w(i)*(1.0_dp+rng_normal(rng))
      end do
    end do
  end subroutine random_product_weights

  function span_as(bench,test,ks,l_values,b_draws,seed) result(res)
    real(dp), intent(in) :: bench(:,:),test(:,:)
    real(dp), intent(in), optional :: ks(:)
    integer, intent(in), optional :: l_values(:),b_draws,seed
    type(as_result) :: res
    real(dp), allocatable :: kvals(:),x1(:),y(:,:),xc(:,:),xaug(:,:),main(:,:),bd(:,:),ba(:,:)
    real(dp), allocatable :: resid(:,:),x1m(:,:),x1p(:,:),yp1(:,:),onem(:,:),onep(:,:),yp2(:,:)
    real(dp), allocatable :: ew(:,:),ew1(:,:),dscore(:,:),ascore(:,:),w(:),dw(:,:),aw(:,:)
    real(dp), allocatable :: pd(:),pa(:),pad(:),draw_d(:,:),draw_a(:,:),draw_ad(:,:)
    real(dp), allocatable :: final_d(:,:),final_a(:,:),final_ad(:,:),merged(:)
    integer, allocatable :: lvals(:)
    integer :: n,k,m,nk,nl,nb,base_b,seed0,info,i,j,li,ki,b,idx
    real(dp) :: denom
    character(len=16) :: sl,sk

    if (size(bench,1)/=size(test,1) .or. size(bench,1)<1 .or. size(bench,2)<1 .or. size(test,2)<1) then
      res%status=span_invalid_input; res%message='incompatible or empty return matrices'; return
    end if
    n=size(bench,1); k=size(bench,2); m=size(test,2)
    if (present(ks)) then
      allocate(kvals(size(ks))); kvals=ks
    else
      allocate(kvals(1)); kvals=1.0_dp/3.0_dp
    end if
    if (present(l_values)) then
      allocate(lvals(size(l_values))); lvals=l_values
    else
      allocate(lvals(2)); lvals=[0,2]
    end if
    nb=1; if (present(b_draws)) nb=b_draws
    seed0=123; if (present(seed)) seed0=seed
    if (nb<1 .or. any(kvals<=0.0_dp) .or. any(kvals>1.0_dp) .or. any(lvals<0)) then
      res%status=span_invalid_input; res%message='invalid ks, L, or B controls'; return
    end if
    nk=size(kvals); nl=size(lvals)

    allocate(x1(n),y(n,m),xc(n,max(0,k-1)),xaug(n,k),main(n,k+1),bd(n,k),ba(n,k))
    x1=bench(:,1); y=test-spread(x1,2,m)
    xaug(:,1)=x1
    if (k>1) then
      do j=2,k
        xc(:,j-1)=bench(:,j)-x1
        xaug(:,j)=xc(:,j-1)
      end do
    end if
    main(:,1)=1.0_dp; main(:,2:k+1)=xaug
    bd(:,1)=1.0_dp
    if (k>1) bd(:,2:k)=xc
    ba=xaug

    allocate(resid(n,m))
    call residualize_matrix(main,y,resid,info)
    if (info/=0) then
      res%status=span_singular; res%message='singular main subseries design'; return
    end if

    allocate(x1m(n,1),x1p(n,1),yp1(n,m),onem(n,1),onep(n,1),yp2(n,m))
    x1m(:,1)=x1; onem(:,1)=1.0_dp
    call residualize_matrix(bd,x1m,x1p,info)
    if (info/=0) then
      res%status=span_singular; res%message='singular delta base design'; return
    end if
    call residualize_matrix(bd,y,yp1,info)
    if (info/=0) then
      res%status=span_singular; res%message='singular delta partialling design'; return
    end if
    allocate(ew(n,m))
    do j=1,m
      denom=dot_product(yp1(:,j),yp1(:,j))
      if (denom<=tiny(1.0_dp)) then
        res%status=span_singular; res%message='degenerate delta partial residual'; return
      end if
      ew(:,j)=x1p(:,1)-yp1(:,j)*(dot_product(yp1(:,j),x1p(:,1))/denom)
    end do

    call residualize_matrix(ba,onem,onep,info)
    if (info/=0) then
      res%status=span_singular; res%message='singular alpha base design'; return
    end if
    call residualize_matrix(ba,y,yp2,info)
    if (info/=0) then
      res%status=span_singular; res%message='singular alpha partialling design'; return
    end if
    allocate(ew1(n,m))
    do j=1,m
      denom=dot_product(yp2(:,j),yp2(:,j))
      if (denom<=tiny(1.0_dp)) then
        res%status=span_singular; res%message='degenerate alpha partial residual'; return
      end if
      ew1(:,j)=onep(:,1)-yp2(:,j)*(dot_product(yp2(:,j),onep(:,1))/denom)
    end do

    allocate(dscore(n,m),ascore(n,m),w(n),dw(n,m),aw(n,m),pd(m),pa(m),pad(m))
    dscore=resid*ew; ascore=resid*ew1
    allocate(final_d(nl,nk),final_a(nl,nk),final_ad(nl,nk),merged(m))

    do li=1,nl
      base_b=1
      if (lvals(li)>0) base_b=nb
      do ki=1,nk
        allocate(draw_d(m,base_b),draw_a(m,base_b),draw_ad(m,base_b))
        do b=1,base_b
          call random_product_weights(n,lvals(li),seed0+b-1,w)
          dw=dscore*spread(w,2,m); aw=ascore*spread(w,2,m)
          call student_subseries_pvalues(dw,kvals(ki),pd,info)
          if (info/=0) then
            res%status=span_numerical_failure; res%message='sample too short for requested subseries exponent'; return
          end if
          call student_subseries_pvalues(aw,kvals(ki),pa,info)
          if (info/=0) then
            res%status=span_numerical_failure; res%message='sample too short for requested subseries exponent'; return
          end if
          do j=1,m
            pad(j)=cauchy_pvalue([pd(j),pa(j)])
          end do
          draw_d(:,b)=pd; draw_a(:,b)=pa; draw_ad(:,b)=pad
        end do
        do j=1,m
          merged(j)=cauchy_pvalue(draw_d(j,:))
        end do
        final_d(li,ki)=cauchy_pvalue(merged)
        do j=1,m
          merged(j)=cauchy_pvalue(draw_a(j,:))
        end do
        final_a(li,ki)=cauchy_pvalue(merged)
        do j=1,m
          merged(j)=cauchy_pvalue(draw_ad(j,:))
        end do
        final_ad(li,ki)=cauchy_pvalue(merged)
        deallocate(draw_d,draw_a,draw_ad)
      end do
    end do

    allocate(res%pvalues(3*nl*nk),res%names(3*nl*nk))
    idx=0
    do i=1,3
      do li=1,nl
        write(sl,'(i0)') lvals(li)
        do ki=1,nk
          write(sk,'(i0)') ki
          idx=idx+1
          select case(i)
          case(1)
            res%names(idx)='CCTd_L'//trim(sl)//'_k'//trim(sk); res%pvalues(idx)=final_d(li,ki)
          case(2)
            res%names(idx)='CCTad_L'//trim(sl)//'_k'//trim(sk); res%pvalues(idx)=final_ad(li,ki)
          case(3)
            res%names(idx)='CCTa_L'//trim(sl)//'_k'//trim(sk); res%pvalues(idx)=final_a(li,ki)
          end select
        end do
      end do
    end do
    res%status=span_ok
  end function span_as

end module spantest_as
