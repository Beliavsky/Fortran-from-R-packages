module dice_design_uniformity
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use dice_design_kinds, only : dp
  use dice_design_utils, only : sort_real, normal_quantile
  implicit none
  private

  type, public :: rss2d_result
    real(dp), allocatable :: global_stat(:, :)
    integer :: worst_case(2) = 0
    real(dp) :: worst_dir(2) = 0.0_dp
    real(dp), allocatable :: stat(:)
    real(dp), allocatable :: angle(:)
    real(dp), allocatable :: curve(:, :)
    real(dp) :: gof_test_stat = 0.0_dp
  end type rss2d_result

  type, public :: rss3d_result
    real(dp), allocatable :: global_stat(:, :, :)
    integer :: worst_case(3) = 0
    real(dp) :: worst_dir(3) = 0.0_dp
    real(dp), allocatable :: stat(:, :)
    real(dp), allocatable :: theta(:), phi(:)
    real(dp) :: gof_test_stat = 0.0_dp
  end type rss3d_result

  public :: unif_test_statistic, unif_test_quantile, rss2d, rss3d
  public :: uniform_cdf, sumof2uniforms_cdf, sumof3uniforms_cdf

contains

  pure function uniform_cdf(p) result(f)
    real(dp), intent(in) :: p
    real(dp) :: f
    f = 0.5_dp*(max(p+1.0_dp,0.0_dp)-max(p-1.0_dp,0.0_dp))
  end function uniform_cdf

  pure function sumof2uniforms_cdf(p, ax, ay) result(f)
    real(dp), intent(in) :: p, ax, ay
    real(dp) :: f, s
    real(dp), parameter :: precision = 1.0e-12_dp

    if (abs(ax)>precision .and. abs(ay)>precision) then
      s = max(p-ax-ay,0.0_dp)**2 + max(p+ax+ay,0.0_dp)**2 - &
          max(p-ax+ay,0.0_dp)**2 - max(p+ax-ay,0.0_dp)**2
      f = s/(8.0_dp*ax*ay)
    else
      f = uniform_cdf(p)
    end if
  end function sumof2uniforms_cdf

  pure function sumof3uniforms_cdf(p, ax, ay, az) result(f)
    real(dp), intent(in) :: p, ax, ay, az
    real(dp) :: f, s
    real(dp), parameter :: precision = 1.0e-12_dp

    if (abs(ax)>precision .and. abs(ay)>precision .and. abs(az)>precision) then
      s = max(p+ax+ay+az,0.0_dp)**3 + max(p+ax-ay-az,0.0_dp)**3 + &
          max(p-ax+ay-az,0.0_dp)**3 + max(p-ax-ay+az,0.0_dp)**3 - &
          max(p-ax-ay-az,0.0_dp)**3 - max(p-ax+ay+az,0.0_dp)**3 - &
          max(p+ax-ay+az,0.0_dp)**3 - max(p+ax+ay-az,0.0_dp)**3
      f = s/(48.0_dp*ax*ay*az)
    else if (abs(ax)<=precision) then
      f = sumof2uniforms_cdf(p,ay,az)
    else if (abs(ay)<=precision) then
      f = sumof2uniforms_cdf(p,ax,az)
    else
      f = sumof2uniforms_cdf(p,ax,ay)
    end if
  end function sumof3uniforms_cdf

  function unif_test_statistic(x, kind, transform_spacings) result(stat)
    real(dp), intent(in) :: x(:)
    character(len=*), intent(in) :: kind
    logical, intent(in), optional :: transform_spacings
    real(dp) :: stat
    real(dp), allocatable :: z(:), s(:), ordered(:), spacings(:), v(:)
    real(dp) :: dplus, dminus
    integer :: n, i
    logical :: trans
    character(len=:), allocatable :: k

    allocate(z(size(x)))
    z = x
    trans = .false.
    if (present(transform_spacings)) trans = transform_spacings
    if (trans) then
      n = size(z)
      allocate(ordered(n),spacings(n+1),s(n+1),v(n+1))
      ordered = z
      call sort_real(ordered)
      spacings(1) = ordered(1)
      do i = 1, n-1
        spacings(i+1) = ordered(i+1)-ordered(i)
      end do
      spacings(n+1) = 1.0_dp-ordered(n)
      s = spacings
      call sort_real(s)
      v(1) = real(n+1,dp)*s(1)
      do i = 1, n
        v(i+1) = v(i) + real(n+1-i,dp)*(s(i+1)-s(i))
      end do
      deallocate(z)
      allocate(z(n+1))
      z = v
      deallocate(ordered,spacings,s,v)
    end if

    k = trim(adjustl(kind))
    n = size(z)
    select case (k)
    case ('greenwood','GREENWOOD')
      allocate(ordered(n),spacings(n+1))
      ordered = z
      call sort_real(ordered)
      spacings(1) = ordered(1)
      do i = 1, n-1
        spacings(i+1) = ordered(i+1)-ordered(i)
      end do
      spacings(n+1) = 1.0_dp-ordered(n)
      stat = sum(spacings**2)
    case ('spacings.max','SPACINGS.MAX')
      allocate(ordered(n),spacings(n+1))
      ordered = z
      call sort_real(ordered)
      spacings(1) = ordered(1)
      do i = 1, n-1
        spacings(i+1) = ordered(i+1)-ordered(i)
      end do
      spacings(n+1) = 1.0_dp-ordered(n)
      stat = maxval(abs(spacings-1.0_dp/real(n,dp)))
    case ('qm','QM')
      allocate(ordered(n),spacings(n+1))
      ordered = z
      call sort_real(ordered)
      spacings(1) = ordered(1)
      do i = 1, n-1
        spacings(i+1) = ordered(i+1)-ordered(i)
      end do
      spacings(n+1) = 1.0_dp-ordered(n)
      stat = sum(spacings**2)+sum(spacings(1:n)*spacings(2:n+1))
    case ('ks','KS')
      allocate(ordered(n))
      ordered = z
      call sort_real(ordered)
      dplus = 0.0_dp
      dminus = 0.0_dp
      do i = 1, n
        dplus = max(dplus,real(i,dp)/real(n,dp)-ordered(i))
        dminus = max(dminus,ordered(i)-real(i-1,dp)/real(n,dp))
      end do
      stat = max(dplus,dminus)
    case ('V','v')
      allocate(ordered(n))
      ordered = z
      call sort_real(ordered)
      dplus = 0.0_dp
      dminus = 0.0_dp
      do i = 1, n
        dplus = max(dplus,real(i,dp)/real(n,dp)-ordered(i))
        dminus = max(dminus,ordered(i)-real(i-1,dp)/real(n,dp))
      end do
      stat = dplus+dminus
    case ('cvm','CVM')
      allocate(ordered(n))
      ordered = z
      call sort_real(ordered)
      stat = 1.0_dp/(12.0_dp*real(n,dp))
      do i = 1, n
        stat = stat+(ordered(i)-real(2*i-1,dp)/(2.0_dp*real(n,dp)))**2
      end do
    case default
      error stop 'unif_test_statistic: unknown test type'
    end select
  end function unif_test_statistic

  function unif_test_quantile(kind, n, alpha) result(q)
    character(len=*), intent(in) :: kind
    integer, intent(in) :: n
    real(dp), intent(in), optional :: alpha
    real(dp) :: q, a, qtab
    integer :: col, i
    character(len=:), allocatable :: k
    integer, parameter :: ntab(23) = [ &
      2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,40,50,60,80,100,200,500 ]
    real(dp), parameter :: q10(23) = [ &
      1.381_dp,1.635_dp,1.800_dp,1.915_dp,1.995_dp,2.053_dp,2.097_dp, &
      2.131_dp,2.157_dp,2.204_dp,2.227_dp,2.242_dp,2.251_dp,2.258_dp, &
      2.265_dp,2.265_dp,2.258_dp,2.248_dp,2.238_dp,2.220_dp,2.205_dp, &
      2.159_dp,2.107_dp ]
    real(dp), parameter :: q05(23) = [ &
      1.539_dp,1.852_dp,2.037_dp,2.160_dp,2.246_dp,2.306_dp,2.349_dp, &
      2.381_dp,2.404_dp,2.441_dp,2.457_dp,2.464_dp,2.466_dp,2.465_dp, &
      2.456_dp,2.443_dp,2.415_dp,2.389_dp,2.367_dp,2.331_dp,2.304_dp, &
      2.226_dp,2.147_dp ]
    real(dp), parameter :: q025(23) = [ &
      1.673_dp,2.075_dp,2.311_dp,2.461_dp,2.559_dp,2.615_dp,2.670_dp, &
      2.700_dp,2.717_dp,2.683_dp,2.691_dp,2.691_dp,2.685_dp,2.677_dp, &
      2.651_dp,2.624_dp,2.573_dp,2.531_dp,2.495_dp,2.441_dp,2.400_dp, &
      2.289_dp,2.183_dp ]
    real(dp), parameter :: q01(23) = [ &
      1.780_dp,2.269_dp,2.560_dp,2.737_dp,2.849_dp,2.921_dp,2.967_dp, &
      2.997_dp,3.008_dp,3.015_dp,3.014_dp,3.003_dp,2.988_dp,2.970_dp, &
      2.920_dp,2.873_dp,2.790_dp,2.723_dp,2.669_dp,2.587_dp,2.528_dp, &
      2.371_dp,2.228_dp ]
    real(dp), parameter :: ks_q(4) = [1.224_dp,1.358_dp,1.480_dp,1.628_dp]
    real(dp), parameter :: v_q(4) = [1.620_dp,1.747_dp,1.862_dp,2.001_dp]
    real(dp), parameter :: cvm_q(4) = [0.347_dp,0.461_dp,0.581_dp,0.743_dp]
    real(dp) :: qcol(23)

    if (n <= 0) error stop 'unif_test_quantile: n must be positive'
    a = 0.05_dp
    if (present(alpha)) a = alpha
    if (abs(a-0.1_dp)<1e-14_dp) then
      col = 1
    else if (abs(a-0.05_dp)<1e-14_dp) then
      col = 2
    else if (abs(a-0.025_dp)<1e-14_dp) then
      col = 3
    else if (abs(a-0.01_dp)<1e-14_dp) then
      col = 4
    else
      error stop 'unif_test_quantile: alpha must be 0.1, 0.05, 0.025, or 0.01'
    end if
    select case (col)
    case (1)
      qcol = q10
    case (2)
      qcol = q05
    case (3)
      qcol = q025
    case (4)
      qcol = q01
    end select
    k = trim(adjustl(kind))
    select case (k)
    case ('greenwood','GREENWOOD')
      if (n > 500) then
        q = (2.0_dp*real(n,dp)/real(n+2,dp) + sqrt(4.0_dp/real(n,dp))*normal_quantile(1.0_dp-a/2.0_dp))/real(n,dp)
      else
        if (n < ntab(1)) error stop 'unif_test_quantile: Greenwood table starts at n=2'
        if (n == ntab(23)) then
          qtab = qcol(23)
        else
          qtab = qcol(23)
          do i = 1, 22
            if (n >= ntab(i) .and. n <= ntab(i+1)) then
              qtab = qcol(i)+(qcol(i+1)-qcol(i))* &
                real(n-ntab(i),dp)/real(ntab(i+1)-ntab(i),dp)
              exit
            end if
          end do
        end if
        q = qtab/real(n,dp)
      end if
    case ('qm','QM')
      error stop 'unif_test_quantile: no default value for Quesenberry-Miller statistic'
    case ('ks','KS')
      q = ks_q(col)/(sqrt(real(n,dp))+0.12_dp+0.11_dp/sqrt(real(n,dp)))
    case ('V','v')
      q = v_q(col)/(sqrt(real(n,dp))+0.155_dp+0.24_dp/sqrt(real(n,dp)))
    case ('cvm','CVM')
      q = cvm_q(col)/(1.0_dp+1.0_dp/real(n,dp)) + &
          0.4_dp/real(n,dp)-0.6_dp/real(n*n,dp)
    case default
      error stop 'unif_test_quantile: unknown test type'
    end select
  end function unif_test_quantile

  subroutine rss2d(design_in, lower, upper, result, gof_test_type, gof_test_stat, transform_spacings, n_angle)
    real(dp), intent(in) :: design_in(:, :), lower(:), upper(:)
    type(rss2d_result), intent(out) :: result
    character(len=*), intent(in), optional :: gof_test_type
    real(dp), intent(in), optional :: gof_test_stat
    logical, intent(in), optional :: transform_spacings
    integer, intent(in), optional :: n_angle
    real(dp), allocatable :: design(:, :), theta(:), ct(:), st(:), fproj(:), angle_stat(:)
    real(dp) :: gstat, gmax, p, threshold, nanv
    integer :: n, d, na, i, j, a, kmax
    logical :: trans
    character(len=16) :: kind

    n = size(design_in,1)
    d = size(design_in,2)
    if (d < 2) error stop 'rss2d: design must have at least two columns'
    if (size(lower)/=d .or. size(upper)/=d) error stop 'rss2d: bound size mismatch'
    na = 360
    if (present(n_angle)) na = n_angle
    if (na <= 0) error stop 'rss2d: n_angle must be positive'
    kind = 'greenwood'
    if (present(gof_test_type)) kind = gof_test_type
    trans = .false.
    if (present(transform_spacings)) trans = transform_spacings
    allocate(design(n,d),theta(na),ct(na),st(na),fproj(n),angle_stat(na))
    design = design_in
    do j = 1, d
      if (minval(design(:,j))<lower(j) .or. maxval(design(:,j))>upper(j)) error stop 'rss2d: design outside bounds'
      if (upper(j)<=lower(j)) error stop 'rss2d: invalid bounds'
      design(:,j) = 2.0_dp*((design(:,j)-lower(j))/(upper(j)-lower(j))-0.5_dp)
    end do
    do a = 1, na
      theta(a) = real(a-1,dp)*acos(-1.0_dp)/real(na,dp)
      ct(a) = cos(theta(a))
      st(a) = sin(theta(a))
    end do
    nanv = ieee_value(0.0_dp,ieee_quiet_nan)
    allocate(result%global_stat(d,d))
    result%global_stat = nanv
    gmax = -huge(1.0_dp)
    do i = 1, d-1
      do j = i+1, d
        do a = 1, na
          do kmax = 1, n
            p = design(kmax,i)*ct(a)+design(kmax,j)*st(a)
            fproj(kmax) = sumof2uniforms_cdf(p,ct(a),st(a))
          end do
          angle_stat(a) = unif_test_statistic(fproj,trim(kind),trans)
        end do
        gstat = maxval(angle_stat)
        result%global_stat(i,j) = gstat
        result%global_stat(j,i) = gstat
        if (gstat > gmax) then
          gmax = gstat
          result%worst_case = [i,j]
          if (allocated(result%stat)) deallocate(result%stat)
          allocate(result%stat(na))
          result%stat = angle_stat
        end if
      end do
    end do
    kmax = maxloc(result%stat,dim=1)
    result%worst_dir = [ct(kmax),st(kmax)]
    allocate(result%angle(na),result%curve(2*na,2))
    result%angle = theta
    do a = 1, na
      result%curve(a,1) = result%stat(a)*ct(a)
      result%curve(a,2) = result%stat(a)*st(a)
      result%curve(na+a,:) = -result%curve(a,:)
    end do
    if (present(gof_test_stat)) then
      threshold = gof_test_stat
    else
      threshold = unif_test_quantile(trim(kind),n,0.05_dp)
    end if
    result%gof_test_stat = threshold
  end subroutine rss2d

  subroutine rss3d(design_in, lower, upper, result, gof_test_type, gof_test_stat, transform_spacings, n_angle)
    real(dp), intent(in) :: design_in(:, :), lower(:), upper(:)
    type(rss3d_result), intent(out) :: result
    character(len=*), intent(in), optional :: gof_test_type
    real(dp), intent(in), optional :: gof_test_stat
    logical, intent(in), optional :: transform_spacings
    integer, intent(in), optional :: n_angle
    real(dp), allocatable :: design(:, :), theta(:), phi(:), ct(:), st(:), cp(:), sp(:)
    real(dp), allocatable :: angle_stat(:, :), best_stat(:, :), fproj(:)
    real(dp) :: ax, ay, az, p, gstat, gmax, nanv
    integer :: n, d, na, i1, i2, i3, it, ip, obs, imax(2)
    logical :: trans
    character(len=16) :: kind

    n = size(design_in,1)
    d = size(design_in,2)
    if (d < 3) error stop 'rss3d: design must have at least three columns'
    if (size(lower)/=d .or. size(upper)/=d) error stop 'rss3d: bound size mismatch'
    na = 60
    if (present(n_angle)) na = n_angle
    if (na <= 0) error stop 'rss3d: n_angle must be positive'
    kind = 'greenwood'
    if (present(gof_test_type)) kind = gof_test_type
    trans = .false.
    if (present(transform_spacings)) trans = transform_spacings
    allocate(design(n,d),theta(na),phi(na),ct(na),st(na),cp(na),sp(na),angle_stat(na,na),best_stat(na,na),fproj(n))
    design = design_in
    do i1 = 1, d
      if (minval(design(:,i1))<lower(i1) .or. maxval(design(:,i1))>upper(i1)) error stop 'rss3d: design outside bounds'
      if (upper(i1)<=lower(i1)) error stop 'rss3d: invalid bounds'
      design(:,i1) = 2.0_dp*((design(:,i1)-lower(i1))/(upper(i1)-lower(i1))-0.5_dp)
    end do
    do it = 1, na
      theta(it) = real(it-1,dp)*acos(-1.0_dp)/real(na,dp)
      ct(it) = cos(theta(it))
      st(it) = sin(theta(it))
      phi(it) = -0.5_dp*acos(-1.0_dp)+real(it-1,dp)*acos(-1.0_dp)/real(na+9,dp)
      cp(it) = cos(phi(it))
      sp(it) = sin(phi(it))
    end do
    nanv = ieee_value(0.0_dp,ieee_quiet_nan)
    allocate(result%global_stat(d,d,d))
    result%global_stat = nanv
    gmax = -huge(1.0_dp)
    do i1 = 1, d-2
      do i2 = i1+1, d-1
        do i3 = i2+1, d
          do ip = 1, na
            do it = 1, na
              ax = ct(it)*cp(ip)
              ay = st(it)*cp(ip)
              az = sp(ip)
              do obs = 1, n
                p = design(obs,i1)*ax+design(obs,i2)*ay+design(obs,i3)*az
                fproj(obs) = sumof3uniforms_cdf(p,ax,ay,az)
              end do
              angle_stat(it,ip) = unif_test_statistic(fproj,trim(kind),trans)
            end do
          end do
          gstat = maxval(angle_stat)
          result%global_stat(i1,i2,i3) = gstat
          if (gstat > gmax) then
            gmax = gstat
            result%worst_case = [i1,i2,i3]
            best_stat = angle_stat
          end if
        end do
      end do
    end do
    imax = maxloc(best_stat)
    ax = ct(imax(1))*cp(imax(2))
    ay = st(imax(1))*cp(imax(2))
    az = sp(imax(2))
    result%worst_dir = [ax,ay,az]
    allocate(result%stat(na,na),result%theta(na),result%phi(na))
    result%stat = best_stat
    result%theta = theta
    result%phi = phi
    if (present(gof_test_stat)) then
      result%gof_test_stat = gof_test_stat
    else
      result%gof_test_stat = unif_test_quantile(trim(kind),n,0.05_dp)
    end if
  end subroutine rss3d

end module dice_design_uniformity
