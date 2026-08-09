! SPDX-License-Identifier: BSD-2-Clause
module smoof_single
  use smoof_kinds, only : dp, pi
  implicit none
  private

  public :: ackley, adjiman, alpine01, alpine02, aluffi_pentini
  public :: bartels_conn, beale, bent_cigar, bird, bohachevsky_n1, booth
  public :: branin, brent, brown, bukin_n2, bukin_n4, bukin_n6
  public :: carrom_table, chichinadze, chung_reynolds, complex_fn
  public :: cosine_mixture, cross_in_tray, cube, deckkers_aarts
  public :: deflected_corrugated_spring, dixon_price, double_sum, drop_wave
  public :: easom, eggcrate, eggholder, el_attar_vidyasagar_dutta, engvall
  public :: exponential_fn, freudenstein_roth, giunta, goldstein_price
  public :: griewank, hansen, hartmann, himmelblau, holder_table_n1
  public :: holder_table_n2, hosaki, hyper_ellipsoid, jennrich_sampson
  public :: judge, keane, kearfott, leon, matyas, mccormick, michalewicz
  public :: modified_rastrigin, periodic_fn, powell_sum, price_n1, price_n2
  public :: price_n4, rastrigin, rosenbrock, schaffer_n2, schaffer_n4
  public :: schwefel, shekel, shubert, six_hump_camel, sphere
  public :: styblinski_tang, sum_different_powers, three_hump_camel
  public :: trecanni, vincent, swiler2014, zettl

contains

  pure real(dp) function ackley(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: d, e
    d = sum(x*x) / real(size(x), dp)
    e = sum(cos(2.0_dp*pi*x)) / real(size(x), dp)
    f = -20.0_dp*exp(-0.2_dp*sqrt(d)) - exp(e) + 20.0_dp + exp(1.0_dp)
  end function ackley

  pure real(dp) function adjiman(x) result(f)
    real(dp), intent(in) :: x(:)
    f = cos(x(1))*sin(x(2)) - x(1)/(x(2)**2 + 1.0_dp)
  end function adjiman

  pure real(dp) function alpine01(x) result(f)
    real(dp), intent(in) :: x(:)
    f = sum(abs(x*sin(x) + 0.1_dp*x))
  end function alpine01

  pure real(dp) function alpine02(x) result(f)
    real(dp), intent(in) :: x(:)
    f = product(sqrt(x)*sin(x))
  end function alpine02

  pure real(dp) function aluffi_pentini(x) result(f)
    real(dp), intent(in) :: x(:)
    f = 0.25_dp*x(1)**4 - 0.5_dp*x(1)**2 + 0.1_dp*x(1) + 0.5_dp*x(2)**2
  end function aluffi_pentini

  pure real(dp) function bartels_conn(x) result(f)
    real(dp), intent(in) :: x(:)
    f = abs(x(1)**2 + x(2)**2 + x(1)*x(2)) + abs(sin(x(1))) + abs(cos(x(2)))
  end function bartels_conn

  pure real(dp) function beale(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: a, b, c
    a = x(1)*x(2)
    b = a*x(2)
    c = b*x(2)
    f = (1.5_dp-x(1)+a)**2 + (2.25_dp-x(1)+b)**2 + (2.625_dp-x(1)+c)**2
  end function beale

  pure real(dp) function bent_cigar(x) result(f)
    real(dp), intent(in) :: x(:)
    if (size(x) == 1) then
      f = x(1)**2
    else
      f = x(1)**2 + 1.0e6_dp*sum(x(2:)**2)
    end if
  end function bent_cigar

  pure real(dp) function bird(x) result(f)
    real(dp), intent(in) :: x(:)
    f = (x(1)-x(2))**2 + exp((1.0_dp-sin(x(1)))**2)*cos(x(2)) &
      + exp((1.0_dp-cos(x(2)))**2)*sin(x(1))
  end function bird

  pure real(dp) function bohachevsky_n1(x) result(f)
    real(dp), intent(in) :: x(:)
    integer :: i
    f = 0.0_dp
    do i = 1, size(x)-1
      f = f + x(i)**2 + 2.0_dp*x(i+1)**2 - 0.3_dp*cos(3.0_dp*pi*x(i)) &
        - 0.4_dp*cos(4.0_dp*pi*x(i+1)) + 0.7_dp
    end do
  end function bohachevsky_n1

  pure real(dp) function booth(x) result(f)
    real(dp), intent(in) :: x(:)
    f = (x(1)+2.0_dp*x(2)-7.0_dp)**2 + (2.0_dp*x(1)+x(2)-5.0_dp)**2
  end function booth

  pure real(dp) function branin(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), parameter :: b = 5.1_dp/(4.0_dp*pi*pi), c = 5.0_dp/pi
    f = (x(2)-b*x(1)**2+c*x(1)-6.0_dp)**2 &
      + 10.0_dp*(1.0_dp-1.0_dp/(8.0_dp*pi))*cos(x(1)) + 10.0_dp
  end function branin

  pure real(dp) function brent(x) result(f)
    real(dp), intent(in) :: x(:)
    f = (x(1)+10.0_dp)**2 + (x(2)+10.0_dp)**2 + exp(-x(1)**2-x(2)**2)
  end function brent

  pure real(dp) function brown(x) result(f)
    real(dp), intent(in) :: x(:)
    integer :: i
    real(dp) :: a, b
    f = 0.0_dp
    do i = 1, size(x)-1
      a = x(i)**2
      b = x(i+1)**2
      f = f + a**(b+1.0_dp) + b**(a+1.0_dp)
    end do
  end function brown

  pure real(dp) function bukin_n2(x) result(f)
    real(dp), intent(in) :: x(:)
    f = 100.0_dp*(x(2)**2 - 0.01_dp*x(1)**2 + 1.0_dp) + 0.01_dp*(x(1)+10.0_dp)**2
  end function bukin_n2

  pure real(dp) function bukin_n4(x) result(f)
    real(dp), intent(in) :: x(:)
    f = 100.0_dp*x(2)**2 + 0.01_dp*abs(x(1)+10.0_dp)
  end function bukin_n4

  pure real(dp) function bukin_n6(x) result(f)
    real(dp), intent(in) :: x(:)
    f = 100.0_dp*sqrt(abs(x(2)-0.01_dp*x(1)**2)) + 0.01_dp*abs(x(1)+10.0_dp)
  end function bukin_n6

  pure real(dp) function carrom_table(x) result(f)
    real(dp), intent(in) :: x(:)
    f = -(1.0_dp/30.0_dp)*exp(2.0_dp*abs(1.0_dp-sqrt(sum(x(1:2)**2))/pi)) &
      * cos(x(1))**2*cos(x(2))**2
  end function carrom_table

  pure real(dp) function chichinadze(x) result(f)
    real(dp), intent(in) :: x(:)
    f = x(1)**2 - 12.0_dp*x(1) + 11.0_dp + 10.0_dp*cos(0.5_dp*pi*x(1)) &
      + 8.0_dp*sin(5.0_dp*pi*x(1)) - exp(-0.5_dp*(x(2)-0.5_dp)**2)/sqrt(5.0_dp)
  end function chichinadze

  pure real(dp) function chung_reynolds(x) result(f)
    real(dp), intent(in) :: x(:)
    f = sum(x*x)**2
  end function chung_reynolds

  pure real(dp) function complex_fn(x) result(f)
    real(dp), intent(in) :: x(:)
    f = (x(1)**3-3.0_dp*x(1)*x(2)**2-1.0_dp)**2 &
      + (3.0_dp*x(2)*x(1)**2-x(2)**3)**2
  end function complex_fn

  pure real(dp) function cosine_mixture(x) result(f)
    real(dp), intent(in) :: x(:)
    f = 0.1_dp*sum(cos(5.0_dp*pi*x)) - sum(x*x)
  end function cosine_mixture

  pure real(dp) function cross_in_tray(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: a
    a = exp(abs(100.0_dp-sqrt(x(1)**2+x(2)**2)/pi))
    f = -0.0001_dp*(abs(a*sin(x(1))*sin(x(2)))+1.0_dp)**0.1_dp
  end function cross_in_tray

  pure real(dp) function cube(x) result(f)
    real(dp), intent(in) :: x(:)
    f = 100.0_dp*(x(2)-x(1)**3)**2 + (1.0_dp-x(1))**2
  end function cube

  pure real(dp) function deckkers_aarts(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: a, b, s
    a = x(1)**2; b = x(2)**2; s = a+b
    f = 1.0e5_dp*a + b - s**2 + 1.0e-5_dp*s**4
  end function deckkers_aarts

  pure real(dp) function deflected_corrugated_spring(x, k, alpha) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: k, alpha
    real(dp) :: kk, aa, s
    kk = 5.0_dp; aa = 5.0_dp
    if (present(k)) kk = k
    if (present(alpha)) aa = alpha
    s = sum((x-aa)**2)
    f = 0.1_dp*s - cos(kk*sqrt(s))
  end function deflected_corrugated_spring

  pure real(dp) function dixon_price(x) result(f)
    real(dp), intent(in) :: x(:)
    integer :: i
    f = (x(1)-1.0_dp)**2
    do i = 2, size(x)
      f = f + real(i,dp)*(2.0_dp*x(i)**2-x(i-1))**2
    end do
  end function dixon_price

  pure real(dp) function double_sum(x) result(f)
    real(dp), intent(in) :: x(:)
    integer :: i
    real(dp) :: s
    f=0.0_dp; s=0.0_dp
    do i=1,size(x)
      s=s+x(i); f=f+s*s
    end do
  end function double_sum

  pure real(dp) function drop_wave(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: a
    a = sum(x*x)
    f = -(1.0_dp+cos(12.0_dp*sqrt(a)))/(0.5_dp*a+2.0_dp)
  end function drop_wave

  pure real(dp) function easom(x) result(f)
    real(dp), intent(in) :: x(:)
    f = -cos(x(1))*cos(x(2))*exp(-((x(1)-pi)**2+(x(2)-pi)**2))
  end function easom

  pure real(dp) function eggcrate(x) result(f)
    real(dp), intent(in) :: x(:)
    f = x(1)**2+x(2)**2+25.0_dp*(sin(x(1))**2+sin(x(2))**2)
  end function eggcrate

  pure real(dp) function eggholder(x) result(f)
    real(dp), intent(in) :: x(:)
    f = -(x(2)+47.0_dp)*sin(sqrt(abs(x(2)+0.5_dp*x(1)+47.0_dp))) &
      - x(1)*sin(sqrt(abs(x(1)-(x(2)+47.0_dp))))
  end function eggholder

  pure real(dp) function el_attar_vidyasagar_dutta(x) result(f)
    real(dp), intent(in) :: x(:)
    f = (x(1)**2+x(2)-10.0_dp)**2 + (x(1)+x(2)**2-7.0_dp)**2 &
      + (x(1)**2+x(2)**3-1.0_dp)**2
  end function el_attar_vidyasagar_dutta

  pure real(dp) function engvall(x) result(f)
    real(dp), intent(in) :: x(:)
    f = x(1)**4+x(2)**4+2.0_dp*x(1)**2*x(2)**2-4.0_dp*x(1)+3.0_dp
  end function engvall

  pure real(dp) function exponential_fn(x) result(f)
    real(dp), intent(in) :: x(:)
    f = -exp(-0.5_dp*sum(x*x))
  end function exponential_fn

  pure real(dp) function freudenstein_roth(x) result(f)
    real(dp), intent(in) :: x(:)
    f = (x(1)-13.0_dp+((5.0_dp-x(2))*x(2)-2.0_dp)*x(2))**2 &
      + (x(1)-29.0_dp+((x(2)+1.0_dp)*x(2)-14.0_dp)*x(2))**2
  end function freudenstein_roth

  pure real(dp) function giunta(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: a(size(x)), b(size(x))
    a=1.067_dp*x-1.0_dp; b=sin(a)
    f=0.6_dp+sum(b+b*b+0.02_dp*sin(4.0_dp*a))
  end function giunta

  pure real(dp) function goldstein_price(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: a,b,xx1,xx2,xx12
    xx1=x(1)**2; xx2=x(2)**2; xx12=x(1)*x(2)
    a=1.0_dp+(x(1)+x(2)+1.0_dp)**2*(19.0_dp-14.0_dp*x(1)+3.0_dp*xx1 &
      -14.0_dp*x(2)+6.0_dp*xx12+3.0_dp*xx2)
    b=30.0_dp+(2.0_dp*x(1)-3.0_dp*x(2))**2*(18.0_dp-32.0_dp*x(1)+12.0_dp*xx1 &
      +48.0_dp*x(2)-36.0_dp*xx12+27.0_dp*xx2)
    f=a*b
  end function goldstein_price

  pure real(dp) function griewank(x) result(f)
    real(dp), intent(in) :: x(:)
    integer :: i
    real(dp) :: p
    p=1.0_dp
    do i=1,size(x)
      p=p*cos(x(i)/sqrt(real(i,dp)))
    end do
    f=sum(x*x)/4000.0_dp-p+1.0_dp
  end function griewank

  pure real(dp) function hansen(x) result(f)
    real(dp), intent(in) :: x(:)
    integer :: i
    real(dp) :: a,b
    a=0.0_dp; b=0.0_dp
    do i=0,4
      a=a+real(i+1,dp)*cos(real(i,dp)*x(1)+real(i+1,dp))
      b=b+real(i+1,dp)*cos(real(i+2,dp)*x(2)+real(i+1,dp))
    end do
    f=a*b
  end function hansen

  pure real(dp) function hartmann(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), parameter :: alpha(4)=[1.0_dp,1.2_dp,3.0_dp,3.2_dp]
    real(dp), parameter :: a3(4,3)=reshape([ &
      3.0_dp,0.1_dp,3.0_dp,0.1_dp,10.0_dp,10.0_dp,10.0_dp,10.0_dp, &
      30.0_dp,35.0_dp,30.0_dp,35.0_dp],[4,3])
    real(dp), parameter :: p3(4,3)=reshape([ &
      0.36890_dp,0.46990_dp,0.10910_dp,0.03815_dp, &
      0.11700_dp,0.43870_dp,0.87320_dp,0.57430_dp, &
      0.26730_dp,0.74700_dp,0.55470_dp,0.88280_dp],[4,3])
    real(dp), parameter :: a6(4,6)=reshape([ &
      10.0_dp,0.05_dp,3.0_dp,17.0_dp,3.0_dp,10.0_dp,3.5_dp,8.0_dp, &
      17.0_dp,17.0_dp,1.7_dp,0.05_dp,3.5_dp,0.1_dp,10.0_dp,10.0_dp, &
      1.7_dp,8.0_dp,17.0_dp,0.1_dp,8.0_dp,14.0_dp,8.0_dp,14.0_dp],[4,6])
    real(dp), parameter :: p6(4,6)=reshape([ &
      0.1312_dp,0.2329_dp,0.2348_dp,0.4047_dp,0.1696_dp,0.4135_dp,0.1451_dp,0.8828_dp, &
      0.5569_dp,0.8307_dp,0.3522_dp,0.8732_dp,0.0124_dp,0.3736_dp,0.2883_dp,0.5743_dp, &
      0.8283_dp,0.1004_dp,0.3047_dp,0.1091_dp,0.5886_dp,0.9991_dp,0.6650_dp,0.0381_dp],[4,6])
    integer :: i,j,n
    real(dp) :: s
    n=size(x); f=0.0_dp
    do i=1,4
      s=0.0_dp
      do j=1,n
        if (n==3) then
          s=s+a3(i,j)*(x(j)-p3(i,j))**2
        else
          s=s+a6(i,j)*(x(j)-p6(i,j))**2
        end if
      end do
      f=f-alpha(i)*exp(-s)
    end do
  end function hartmann

  pure real(dp) function himmelblau(x) result(f)
    real(dp), intent(in) :: x(:)
    f=(x(1)**2+x(2)-11.0_dp)**2+(x(1)+x(2)**2-7.0_dp)**2
  end function himmelblau

  pure real(dp) function holder_table_n1(x) result(f)
    real(dp), intent(in) :: x(:)
    f=-abs(cos(x(1))*cos(x(2))*exp(abs(1.0_dp-sqrt(x(1)**2+x(2)**2)/3.1415_dp)))
  end function holder_table_n1

  pure real(dp) function holder_table_n2(x) result(f)
    real(dp), intent(in) :: x(:)
    f=-abs(sin(x(1))*cos(x(2))*exp(abs(1.0_dp-sqrt(x(1)**2+x(2)**2)/3.1415_dp)))
  end function holder_table_n2

  pure real(dp) function hosaki(x) result(f)
    real(dp), intent(in) :: x(:)
    f=(1.0_dp-8.0_dp*x(1)+7.0_dp*x(1)**2-7.0_dp*x(1)**3/3.0_dp+0.25_dp*x(1)**4) &
      *x(2)**2*exp(-x(2))
  end function hosaki

  pure real(dp) function hyper_ellipsoid(x) result(f)
    real(dp), intent(in) :: x(:)
    integer :: i
    f=0.0_dp
    do i=1,size(x); f=f+real(i,dp)*x(i)**2; end do
  end function hyper_ellipsoid

  pure real(dp) function jennrich_sampson(x) result(f)
    real(dp), intent(in) :: x(:)
    integer :: i
    f=0.0_dp
    do i=1,10
      f=f+(2.0_dp+2.0_dp*i-(exp(real(i,dp)*x(1))+exp(real(i,dp)*x(2))))**2
    end do
  end function jennrich_sampson

  pure real(dp) function judge(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), parameter :: a(20)=[4.284_dp,4.149_dp,3.877_dp,0.533_dp,2.211_dp,2.389_dp, &
      2.145_dp,3.231_dp,1.998_dp,1.379_dp,2.106_dp,1.428_dp,1.011_dp,2.179_dp,2.858_dp, &
      1.388_dp,1.651_dp,1.593_dp,1.046_dp,2.152_dp]
    real(dp), parameter :: b(20)=[0.286_dp,0.973_dp,0.384_dp,0.276_dp,0.973_dp,0.543_dp, &
      0.957_dp,0.948_dp,0.543_dp,0.797_dp,0.936_dp,0.889_dp,0.006_dp,0.828_dp,0.399_dp, &
      0.617_dp,0.939_dp,0.784_dp,0.072_dp,0.889_dp]
    real(dp), parameter :: c(20)=[0.645_dp,0.585_dp,0.310_dp,0.058_dp,0.455_dp,0.779_dp, &
      0.259_dp,0.202_dp,0.028_dp,0.099_dp,0.142_dp,0.296_dp,0.175_dp,0.180_dp,0.842_dp, &
      0.039_dp,0.103_dp,0.620_dp,0.158_dp,0.704_dp]
    f=sum((x(1)+b*x(2)+c*x(2)**2-a)**2)
  end function judge

  pure real(dp) function keane(x) result(f)
    real(dp), intent(in) :: x(:)
    f=sin(x(1)-x(2))**2*sin(x(1)+x(2))**2/sqrt(x(1)**2+x(2)**2)
  end function keane

  pure real(dp) function kearfott(x) result(f)
    real(dp), intent(in) :: x(:)
    f=(x(1)**2+x(2)**2-2.0_dp)**2+(x(1)**2-x(2)**2-1.0_dp)**2
  end function kearfott

  pure real(dp) function leon(x) result(f)
    real(dp), intent(in) :: x(:)
    f=100.0_dp*(x(2)-x(1)**2)**2+(1.0_dp-x(1))**2
  end function leon

  pure real(dp) function matyas(x) result(f)
    real(dp), intent(in) :: x(:)
    f=0.26_dp*(x(1)**2+x(2)**2)-0.48_dp*x(1)*x(2)
  end function matyas

  pure real(dp) function mccormick(x) result(f)
    real(dp), intent(in) :: x(:)
    f=sin(x(1)+x(2))+(x(1)-x(2))**2-1.5_dp*x(1)+2.5_dp*x(2)+1.0_dp
  end function mccormick

  pure real(dp) function michalewicz(x, m) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: m
    real(dp) :: mm
    integer :: i
    mm=10.0_dp; if(present(m)) mm=m
    f=0.0_dp
    do i=1,size(x)
      f=f-sin(x(i))*sin(real(i,dp)*x(i)**2/pi)**(2.0_dp*mm)
    end do
  end function michalewicz

  pure real(dp) function modified_rastrigin(x, k) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: k(:)
    if (present(k)) then
      f=sum(10.0_dp*(1.0_dp+cos(2.0_dp*pi*k*x))+2.0_dp*k*x*x)
    else
      f=sum(10.0_dp*(1.0_dp+cos(2.0_dp*pi*x))+2.0_dp*x*x)
    end if
  end function modified_rastrigin

  pure real(dp) function periodic_fn(x) result(f)
    real(dp), intent(in) :: x(:)
    f=1.0_dp+sin(x(1))**2+sin(x(2))**2-0.1_dp*exp(-x(1)**2-x(2)**2)
  end function periodic_fn

  pure real(dp) function powell_sum(x) result(f)
    real(dp), intent(in) :: x(:)
    integer :: i
    f=0.0_dp
    do i=1,size(x); f=f+abs(x(i))**real(i+1,dp); end do
  end function powell_sum

  pure real(dp) function price_n1(x) result(f)
    real(dp), intent(in) :: x(:)
    f=sum((abs(x)-5.0_dp)**2)
  end function price_n1

  pure real(dp) function price_n2(x) result(f)
    real(dp), intent(in) :: x(:)
    f=periodic_fn(x)
  end function price_n2

  pure real(dp) function price_n4(x) result(f)
    real(dp), intent(in) :: x(:)
    f=(2.0_dp*x(1)**3*x(2)-x(2)**3)**2+(6.0_dp*x(1)-x(2)**2+x(2))**2
  end function price_n4

  pure real(dp) function rastrigin(x) result(f)
    real(dp), intent(in) :: x(:)
    f=10.0_dp*real(size(x),dp)+sum(x*x-10.0_dp*cos(2.0_dp*pi*x))
  end function rastrigin

  pure real(dp) function rosenbrock(x) result(f)
    real(dp), intent(in) :: x(:)
    integer :: i
    f=0.0_dp
    do i=1,size(x)-1
      f=f+100.0_dp*(x(i)**2-x(i+1))**2+(x(i)-1.0_dp)**2
    end do
  end function rosenbrock

  pure real(dp) function schaffer_n2(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: a,b
    a=x(1)**2; b=x(2)**2
    f=0.5_dp+(sin(a-b)**2-0.5_dp)/(1.0_dp+0.001_dp*(a+b))**2
  end function schaffer_n2

  pure real(dp) function schaffer_n4(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: a,b
    a=x(1)**2; b=x(2)**2
    f=0.5_dp+(cos(sin(abs(a-b)))**2-0.5_dp)/(1.0_dp+0.001_dp*(a+b))**2
  end function schaffer_n4

  pure real(dp) function schwefel(x) result(f)
    real(dp), intent(in) :: x(:)
    f=sum(-x*sin(sqrt(abs(x))))
  end function schwefel

  pure real(dp) function shekel(x, m) result(f)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: m
    real(dp), parameter :: c(4,10)=reshape([ &
      4.0_dp,4.0_dp,4.0_dp,4.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
      8.0_dp,8.0_dp,8.0_dp,8.0_dp,6.0_dp,6.0_dp,6.0_dp,6.0_dp, &
      3.0_dp,7.0_dp,3.0_dp,7.0_dp,2.0_dp,9.0_dp,2.0_dp,9.0_dp, &
      5.0_dp,3.0_dp,5.0_dp,3.0_dp,8.0_dp,1.0_dp,8.0_dp,1.0_dp, &
      6.0_dp,2.0_dp,6.0_dp,2.0_dp,7.0_dp,3.6_dp,7.0_dp,3.6_dp],[4,10])
    real(dp), parameter :: beta(10)=[0.1_dp,0.2_dp,0.2_dp,0.4_dp,0.4_dp, &
      0.6_dp,0.3_dp,0.7_dp,0.5_dp,0.5_dp]
    integer :: i
    f=0.0_dp
    do i=1,m
      f=f-1.0_dp/(sum((x(1:4)-c(:,i))**2)+beta(i))
    end do
  end function shekel

  pure real(dp) function shubert(x) result(f)
    real(dp), intent(in) :: x(:)
    integer :: j
    real(dp) :: a,b
    a=0.0_dp; b=0.0_dp
    do j=1,5
      a=a+real(j,dp)*cos(real(j+1,dp)*x(1)+real(j,dp))
      b=b+real(j,dp)*cos(real(j+1,dp)*x(2)+real(j,dp))
    end do
    f=a*b
  end function shubert

  pure real(dp) function six_hump_camel(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: a,b
    a=x(1)**2; b=x(2)**2
    f=(4.0_dp-2.1_dp*a+a*a/3.0_dp)*a+x(1)*x(2)+(4.0_dp*b-4.0_dp)*b
  end function six_hump_camel

  pure real(dp) function sphere(x) result(f)
    real(dp), intent(in) :: x(:)
    f=sum(x*x)
  end function sphere

  pure real(dp) function styblinski_tang(x) result(f)
    real(dp), intent(in) :: x(:)
    f=0.5_dp*sum(x**4-16.0_dp*x*x+5.0_dp*x)
  end function styblinski_tang

  pure real(dp) function sum_different_powers(x) result(f)
    real(dp), intent(in) :: x(:)
    integer :: i
    f=0.0_dp
    do i=1,size(x); f=f+abs(x(i))**real(i+1,dp); end do
  end function sum_different_powers

  pure real(dp) function three_hump_camel(x) result(f)
    real(dp), intent(in) :: x(:)
    f=2.0_dp*x(1)**2-1.05_dp*x(1)**4+x(1)**6/6.0_dp+x(1)*x(2)+x(2)**2
  end function three_hump_camel

  pure real(dp) function trecanni(x) result(f)
    real(dp), intent(in) :: x(:)
    f=x(1)**4+4.0_dp*(x(1)**3+x(1)**2)+x(2)**2
  end function trecanni

  pure real(dp) function vincent(x) result(f)
    real(dp), intent(in) :: x(:)
    f=sum(sin(10.0_dp*log(x)))/real(size(x),dp)
  end function vincent

  pure real(dp) function swiler2014(x1,x2,x3) result(f)
    integer, intent(in) :: x1
    real(dp), intent(in) :: x2,x3
    real(dp), parameter :: fac(5)=[0.0_dp,12.0_dp,0.5_dp,8.0_dp,3.5_dp]
    real(dp) :: a,b
    a=sin(2.0_dp*pi*x3-pi)
    b=7.0_dp*sin(2.0_dp*pi*x2-pi)**2
    f=a+b+fac(x1)*a
  end function swiler2014

  pure real(dp) function zettl(x) result(f)
    real(dp), intent(in) :: x(:)
    f=(x(1)**2+x(2)**2-2.0_dp*x(1))**2+0.25_dp*x(1)
  end function zettl

end module smoof_single
