module rcppnumerical_cuhre
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rcppnumerical_kinds, only : dp
  use rcppnumerical_callbacks, only : multivariate_function_interface
  implicit none
  private

  real(dp), parameter :: pi = acos(-1.0_dp)
  integer, parameter :: nrules = 5

  type :: cubature_set_t
    integer :: n = 0
    real(dp) :: weight(nrules) = 0.0_dp
    real(dp) :: scale(nrules) = 0.0_dp
    real(dp) :: norm(nrules) = 0.0_dp
    real(dp), allocatable :: gen(:)
  end type cubature_set_t

  type :: cubature_rule_t
    type(cubature_set_t), allocatable :: sets(:)
    real(dp) :: errcoeff(3) = 0.0_dp
    integer :: n = 0
  end type cubature_rule_t

  type :: region_t
    integer :: div = 0
    real(dp), allocatable :: lower(:), upper(:)
    real(dp) :: avg = 0.0_dp
    real(dp) :: err = 0.0_dp
    integer :: bisectdim = 1
  end type region_t

  type, public :: cubature_result_t
    real(dp) :: value = 0.0_dp
    real(dp) :: error_estimate = huge(1.0_dp)
    integer :: error_code = -1
    integer :: evaluations = 0
    integer :: regions = 0
  contains
    procedure :: successful => cubature_successful
  end type cubature_result_t

  public :: integrate_nd

contains

  pure logical function cubature_successful(self)
    class(cubature_result_t), intent(in) :: self
    cubature_successful = self%error_code == 0
  end function cubature_successful

  subroutine integrate_nd(f, lower, upper, result, maxeval, eps_abs, eps_rel, user_data)
    procedure(multivariate_function_interface) :: f
    real(dp), intent(in) :: lower(:), upper(:)
    type(cubature_result_t), intent(out) :: result
    integer, intent(in), optional :: maxeval
    real(dp), intent(in), optional :: eps_abs, eps_rel
    class(*), intent(inout), optional :: user_data

    type(cubature_rule_t) :: rule
    type(region_t), allocatable :: regions(:)
    type(region_t) :: left_child, right_child, old_region
    integer :: d, max_evaluations, max_regions, nregions, idx, split_dim
    real(dp) :: abs_tol, rel_tol, target, total_avg, total_err
    real(dp) :: diff, combined_err, correction

    result = cubature_result_t()
    d = size(lower)
    max_evaluations = 1000
    if (present(maxeval)) max_evaluations = maxeval
    abs_tol = 1.0e-6_dp
    if (present(eps_abs)) abs_tol = eps_abs
    rel_tol = 1.0e-6_dp
    if (present(eps_rel)) rel_tol = eps_rel

    if (size(upper) /= d .or. d < 2 .or. d > 20 .or. max_evaluations < 1 .or. &
        abs_tol < 0.0_dp .or. rel_tol < 0.0_dp) then
      return
    end if
    if (any(upper <= lower)) then
      result%error_code = -2
      return
    end if
    call build_rule(d, rule)
    if (rule%n <= 0) return
    max_regions = max(2, max_evaluations/max(1, rule%n) + 3)
    allocate(regions(max_regions))

    nregions = 1
    allocate(regions(1)%lower(d), regions(1)%upper(d))
    regions(1)%lower = 0.0_dp
    regions(1)%upper = 1.0_dp
    regions(1)%div = 0
    call sample_region(unit_eval, rule, regions(1), result%evaluations)
    total_avg = regions(1)%avg
    total_err = regions(1)%err
    result%regions = 1

    do
      target = max(abs_tol, rel_tol*abs(total_avg))
      if (total_err <= target) then
        result%error_code = 0
        exit
      end if
      if (result%evaluations >= max_evaluations) then
        result%error_code = 1
        exit
      end if
      idx = maxloc(regions(1:nregions)%err, dim=1)
      old_region = regions(idx)
      split_dim = old_region%bisectdim
      call split_region(old_region, split_dim, left_child, right_child)
      call sample_region(unit_eval, rule, left_child, result%evaluations)
      call sample_region(unit_eval, rule, right_child, result%evaluations)

      diff = abs(0.25_dp*(left_child%avg + right_child%avg - old_region%avg))
      combined_err = left_child%err + right_child%err
      if (combined_err > 0.0_dp) then
        correction = 1.0_dp + 2.0_dp*diff/combined_err
        left_child%err = left_child%err*correction
        right_child%err = right_child%err*correction
      end if
      left_child%err = left_child%err + diff
      right_child%err = right_child%err + diff

      total_avg = total_avg + left_child%avg + right_child%avg - old_region%avg
      total_err = total_err + left_child%err + right_child%err - old_region%err
      regions(idx) = left_child
      nregions = nregions + 1
      if (nregions > size(regions)) then
        result%error_code = 1
        exit
      end if
      regions(nregions) = right_child
      result%regions = nregions
    end do

    result%value = total_avg
    result%error_estimate = total_err

  contains

    function unit_eval(t, ignored) result(value)
      real(dp), intent(in) :: t(:)
      class(*), intent(inout), optional :: ignored
      real(dp) :: value
      real(dp) :: x(d), jacobian, transform, angle
      integer :: i

      jacobian = 1.0_dp
      do i = 1, d
        if (ieee_is_finite(lower(i)) .and. ieee_is_finite(upper(i))) then
          x(i) = lower(i) + (upper(i) - lower(i))*t(i)
          jacobian = jacobian*(upper(i) - lower(i))
        else if (ieee_is_finite(lower(i))) then
          transform = (1.0_dp - t(i))/t(i)
          x(i) = lower(i) + transform
          jacobian = jacobian/(t(i)*t(i))
        else if (ieee_is_finite(upper(i))) then
          transform = (1.0_dp - t(i))/t(i)
          x(i) = upper(i) - transform
          jacobian = jacobian/(t(i)*t(i))
        else
          angle = pi*(t(i) - 0.5_dp)
          x(i) = tan(angle)
          jacobian = jacobian*pi/(cos(angle)*cos(angle))
        end if
      end do
      if (present(user_data)) then
        value = f(x, user_data)*jacobian
      else
        value = f(x)*jacobian
      end if
    end function unit_eval

  end subroutine integrate_nd

  subroutine split_region(parent, dim, left_child, right_child)
    type(region_t), intent(in) :: parent
    integer, intent(in) :: dim
    type(region_t), intent(out) :: left_child, right_child
    real(dp) :: midpoint
    integer :: d

    d = size(parent%lower)
    allocate(left_child%lower(d), left_child%upper(d))
    allocate(right_child%lower(d), right_child%upper(d))
    left_child%lower = parent%lower
    left_child%upper = parent%upper
    right_child%lower = parent%lower
    right_child%upper = parent%upper
    midpoint = 0.5_dp*(parent%lower(dim) + parent%upper(dim))
    left_child%upper(dim) = midpoint
    right_child%lower(dim) = midpoint
    left_child%div = parent%div + 1
    right_child%div = parent%div + 1
  end subroutine split_region

  subroutine sample_region(f, rule, region, evaluations)
    procedure(multivariate_function_interface) :: f
    type(cubature_rule_t), intent(in) :: rule
    type(region_t), intent(inout) :: region
    integer, intent(inout) :: evaluations

    real(dp) :: sums(nrules), set_sum, center_value
    real(dp), allocatable :: axis1_minus(:), axis1_plus(:)
    real(dp), allocatable :: axis2_minus(:), axis2_plus(:)
    real(dp) :: range, maxrange, ratio, base, fourthdiff, maxdiff
    real(dp) :: maxerr, volume
    integer :: s, r, xset, d, dim, maxdim, count_eval

    d = size(region%lower)
    allocate(axis1_minus(d), axis1_plus(d), axis2_minus(d), axis2_plus(d))
    axis1_minus = 0.0_dp
    axis1_plus = 0.0_dp
    axis2_minus = 0.0_dp
    axis2_plus = 0.0_dp
    sums = 0.0_dp
    count_eval = 0
    center_value = 0.0_dp

    maxrange = -1.0_dp
    maxdim = 1
    do dim = 1, d
      range = region%upper(dim) - region%lower(dim)
      if (range > maxrange) then
        maxrange = range
        maxdim = dim
      end if
    end do

    do s = 1, size(rule%sets)
      if (s == 1) then
        call evaluate_center(f, region, set_sum, count_eval)
        center_value = set_sum
      else if (count(rule%sets(s)%gen /= 0.0_dp) == 1) then
        call evaluate_axis_set(f, region, maxval(rule%sets(s)%gen), set_sum, &
                               count_eval, axis1_minus, axis1_plus, &
                               axis2_minus, axis2_plus, s)
      else
        call evaluate_symmetric_set(f, region, rule%sets(s)%gen, set_sum, count_eval)
      end if
      sums = sums + set_sum*rule%sets(s)%weight
    end do
    evaluations = evaluations + count_eval

    ratio = (maxval(rule%sets(3)%gen)/maxval(rule%sets(2)%gen))**2
    base = center_value*2.0_dp*(1.0_dp - ratio)
    maxdiff = 0.0_dp
    region%bisectdim = maxdim
    do dim = 1, d
      fourthdiff = abs(base + ratio*(axis1_minus(dim) + axis1_plus(dim)) - &
                       (axis2_minus(dim) + axis2_plus(dim)))
      if (fourthdiff > maxdiff) then
        maxdiff = fourthdiff
        region%bisectdim = dim
      end if
    end do

    do r = 2, nrules - 1
      maxerr = 0.0_dp
      do xset = 1, size(rule%sets)
        maxerr = max(maxerr, abs(sums(r + 1) + rule%sets(xset)%scale(r)*sums(r))* &
                     rule%sets(xset)%norm(r))
      end do
      sums(r) = maxerr
    end do

    volume = 2.0_dp**(-region%div)
    region%avg = volume*sums(1)
    if (rule%errcoeff(1)*sums(2) <= sums(3) .and. &
        rule%errcoeff(1)*sums(3) <= sums(4)) then
      region%err = volume*rule%errcoeff(2)*sums(2)
    else
      region%err = volume*rule%errcoeff(3)*max(sums(2), sums(3), sums(4))
    end if
  end subroutine sample_region

  subroutine evaluate_center(f, region, value, count_eval)
    procedure(multivariate_function_interface) :: f
    type(region_t), intent(in) :: region
    real(dp), intent(out) :: value
    integer, intent(inout) :: count_eval
    real(dp) :: x(size(region%lower))
    x = 0.5_dp*(region%lower + region%upper)
    value = f(x)
    count_eval = count_eval + 1
  end subroutine evaluate_center

  subroutine evaluate_axis_set(f, region, generator, value, count_eval, &
                               axis1_minus, axis1_plus, axis2_minus, axis2_plus, set_index)
    procedure(multivariate_function_interface) :: f
    type(region_t), intent(in) :: region
    real(dp), intent(in) :: generator
    real(dp), intent(out) :: value
    integer, intent(inout) :: count_eval
    real(dp), intent(inout) :: axis1_minus(:), axis1_plus(:)
    real(dp), intent(inout) :: axis2_minus(:), axis2_plus(:)
    integer, intent(in) :: set_index
    real(dp) :: x(size(region%lower)), center(size(region%lower)), width(size(region%lower))
    real(dp) :: fm, fp
    integer :: dim

    center = 0.5_dp*(region%lower + region%upper)
    width = region%upper - region%lower
    value = 0.0_dp
    do dim = 1, size(center)
      x = center
      x(dim) = center(dim) - generator*width(dim)
      fm = f(x)
      x(dim) = center(dim) + generator*width(dim)
      fp = f(x)
      value = value + fm + fp
      count_eval = count_eval + 2
      if (set_index == 2) then
        axis1_minus(dim) = fm
        axis1_plus(dim) = fp
      else if (set_index == 3) then
        axis2_minus(dim) = fm
        axis2_plus(dim) = fp
      end if
    end do
  end subroutine evaluate_axis_set

  subroutine evaluate_symmetric_set(f, region, generator, value, count_eval)
    procedure(multivariate_function_interface) :: f
    type(region_t), intent(in) :: region
    real(dp), intent(in) :: generator(:)
    real(dp), intent(out) :: value
    integer, intent(inout) :: count_eval
    integer :: d
    real(dp), allocatable :: perm(:), signed_perm(:), center(:), width(:), point(:)
    logical, allocatable :: used(:)

    d = size(generator)
    allocate(perm(d), signed_perm(d), center(d), width(d), point(d), used(d))
    center = 0.5_dp*(region%lower + region%upper)
    width = region%upper - region%lower
    used = .false.
    value = 0.0_dp
    call permute_values(1)

  contains

    recursive subroutine permute_values(pos)
      integer, intent(in) :: pos
      integer :: i, j
      logical :: duplicate
      if (pos > d) then
        signed_perm = perm
        call assign_signs(1)
        return
      end if
      do i = 1, d
        if (used(i)) cycle
        duplicate = .false.
        do j = 1, i - 1
          if (.not. used(j) .and. generator(j) == generator(i)) then
            duplicate = .true.
            exit
          end if
        end do
        if (duplicate) cycle
        used(i) = .true.
        perm(pos) = generator(i)
        call permute_values(pos + 1)
        used(i) = .false.
      end do
    end subroutine permute_values

    recursive subroutine assign_signs(pos)
      integer, intent(in) :: pos
      if (pos > d) then
        point = center - signed_perm*width
        value = value + f(point)
        count_eval = count_eval + 1
        return
      end if
      if (perm(pos) == 0.0_dp) then
        signed_perm(pos) = 0.0_dp
        call assign_signs(pos + 1)
      else
        signed_perm(pos) = perm(pos)
        call assign_signs(pos + 1)
        signed_perm(pos) = -perm(pos)
        call assign_signs(pos + 1)
      end if
    end subroutine assign_signs

  end subroutine evaluate_symmetric_set

  subroutine build_rule(d, rule)
    integer, intent(in) :: d
    type(cubature_rule_t), intent(out) :: rule
    if (d == 2) then
      call build_rule13(rule)
    else if (d == 3) then
      call build_rule11(rule)
    else
      call build_rule9(d, rule)
    end if
    call finalize_rule(rule)
  end subroutine build_rule

  subroutine allocate_sets(rule, nsets, d)
    type(cubature_rule_t), intent(out) :: rule
    integer, intent(in) :: nsets, d
    integer :: i
    allocate(rule%sets(nsets))
    do i = 1, nsets
      allocate(rule%sets(i)%gen(d))
      rule%sets(i)%gen = 0.0_dp
    end do
  end subroutine allocate_sets

  subroutine finalize_rule(rule)
    type(cubature_rule_t), intent(inout) :: rule
    real(dp) :: denom, scale
    integer :: s, r, x
    rule%n = sum(rule%sets%n)
    do s = 1, size(rule%sets)
      do r = 2, nrules - 1
        if (rule%sets(s)%weight(r) == 0.0_dp) then
          scale = 100.0_dp
        else
          scale = -rule%sets(s)%weight(r + 1)/rule%sets(s)%weight(r)
        end if
        denom = 0.0_dp
        do x = 1, size(rule%sets)
          denom = denom + rule%sets(x)%n* &
            abs(rule%sets(x)%weight(r + 1) + scale*rule%sets(x)%weight(r))
        end do
        rule%sets(s)%scale(r) = scale
        if (denom > 0.0_dp) rule%sets(s)%norm(r) = 1.0_dp/denom
      end do
    end do
  end subroutine finalize_rule

  subroutine set_entry(set, n, weight, generators)
    type(cubature_set_t), intent(inout) :: set
    integer, intent(in) :: n
    real(dp), intent(in) :: weight(nrules)
    real(dp), intent(in), optional :: generators(:)
    set%n = n
    set%weight = weight
    if (present(generators)) set%gen(1:size(generators)) = generators
  end subroutine set_entry

  subroutine build_rule13(rule)
    type(cubature_rule_t), intent(out) :: rule
    real(dp), parameter :: w(5,14) = reshape([ &
      .00844923090033615_dp,.3213775489050763_dp,.3372900883288987_dp,-.8264123822525677_dp,.6539094339575232_dp, &
      .023771474018994404_dp,-.1767341636743844_dp,-.1644903060344491_dp,.306583861409436_dp,-.2041614154424632_dp, &
      .02940016170142405_dp,.07347600537466073_dp,.07707849911634623_dp,.002389292538329435_dp,-.174698151579499_dp, &
      .006644436465817374_dp,-.03638022004364754_dp,-.03804478358506311_dp,-.1343024157997222_dp,.03937939671417803_dp, &
      .0042536044255016_dp,.021252979220987123_dp,.02223559940380806_dp,.08833366840533902_dp,.006974520545933992_dp, &
      0.0_dp,.1460984204026913_dp,.1480693879765931_dp,0.0_dp,0.0_dp, &
      .0040664827465935255_dp,.017476132861520992_dp,4.467143702185815e-6_dp,.0009786283074168292_dp,.0066677021717782585_dp, &
      .03362231646315497_dp,.1444954045641582_dp,.150894476707413_dp,-.1319227889147519_dp,.05512960621544304_dp, &
      .033200804136503725_dp,.0001307687976001325_dp,3.6472001075162155e-5_dp,.00799001220015063_dp,.05443846381278608_dp, &
      .014093686924979677_dp,.0005380992313941161_dp,.000577719899901388_dp,.0033917470797606257_dp,.02310903863953934_dp, &
      .000977069770327625_dp,.0001042259576889814_dp,.0001041757313688177_dp,.0022949157182832643_dp,.01506937747477189_dp, &
      .007531996943580376_dp,-.001401152865045733_dp,-.001452822267047819_dp,-.01358584986119197_dp,-.060570216489018905_dp, &
      .02577183086722915_dp,.008041788181514763_dp,.008338339968783704_dp,.04025866859057809_dp,.04225737654686337_dp, &
      .015625_dp,-.1420416552759383_dp,-.147279632923196_dp,.003760268580063992_dp,.02561989142123099_dp ], [5,14])
    real(dp), parameter :: g(16) = [ .12585646717265545_dp,.3506966822267133_dp,.4795480315809981_dp, &
      .4978005239276064_dp,.25_dp,.07972723291487795_dp,.1904495567970094_dp,.3291384627633596_dp, &
      .43807365825146577_dp,.499121592026599_dp,.4895111329084231_dp,.32461421628226944_dp, &
      .43637106005656195_dp,.1791307322940614_dp,.2833333333333333_dp,.1038888888888889_dp ]
    integer :: i
    call allocate_sets(rule, 14, 2)
    call set_entry(rule%sets(1), 1, w(:,1))
    do i = 2, 6
      call set_entry(rule%sets(i), 4, w(:,i), [g(i-1)])
    end do
    do i = 7, 11
      call set_entry(rule%sets(i), 4, w(:,i), [g(i-1), g(i-1)])
    end do
    call set_entry(rule%sets(12), 8, w(:,12), [g(11), g(12)])
    call set_entry(rule%sets(13), 8, w(:,13), [g(13), g(14)])
    call set_entry(rule%sets(14), 8, w(:,14), [g(15), g(16)])
    rule%errcoeff = [10.0_dp, 1.0_dp, 5.0_dp]
  end subroutine build_rule13

  subroutine build_rule11(rule)
    type(cubature_rule_t), intent(out) :: rule
    real(dp), parameter :: w(5,13) = reshape([ &
      .0009903847688882167_dp,1.715006248224684_dp,1.936014978949526_dp,.517082819560576_dp,2.05440450381852_dp, &
      .0084964717409851_dp,-.3755893815889209_dp,-.3673449403754268_dp,.01445269144914044_dp,.013777599884901202_dp, &
      .00013587331735072814_dp,.1488632145140549_dp,.02929778657898176_dp,-.3601489663995932_dp,-.576806291790441_dp, &
      .022982920777660364_dp,-.2497046640620823_dp,-.1151883520260315_dp,.3628307003418485_dp,.03726835047700328_dp, &
      .004202649722286289_dp,.1792501419135204_dp,.05086658220872218_dp,.007148802650872729_dp,.0068148789397772195_dp, &
      .0012671889041675774_dp,.0034461267589738897_dp,.04453911087786469_dp,-.09222852896022966_dp,.057231697338518496_dp, &
      .0002109560854981544_dp,-.005140483185555825_dp,-.022878282571259_dp,.01719339732471725_dp,-.044930187438112855_dp, &
      .016830857056410086_dp,.006536017839876424_dp,.02908926216345833_dp,-.102141653746035_dp,.027292365738663484_dp, &
      .00021876823557504823_dp,-.00065134549392297_dp,-.002898884350669207_dp,-.007504397861080493_dp,.000354747395055699_dp, &
      .009690420479796819_dp,-.006304672433547204_dp,-.028059634133074954_dp,.01648362537726711_dp,.01571366799739551_dp, &
      .030773311284628138_dp,.01266959399788263_dp,.05638741361145884_dp,.05234610158469334_dp,.049900992192785674_dp, &
      .0084974310856038_dp,-.005454241018647931_dp,-.02427469611942451_dp,.014454323316130661_dp,.0137791555266677_dp, &
      .0017749535291258914_dp,.004826995274768427_dp,.021483070341828822_dp, &
      .003019236275367777_dp,.0028782064230998723_dp ], [5,13])
    real(dp), parameter :: g(14) = [ .095_dp,.25_dp,.375_dp,.4_dp,.4975_dp,.49936724991757_dp, &
      .38968518428362114_dp,.49998494965443835_dp,.3951318612385894_dp,.22016983438253684_dp, &
      .4774686911397297_dp,.2189239229503431_dp,.4830546566815374_dp,.2288552938881567_dp ]
    integer :: i
    call allocate_sets(rule, 13, 3)
    call set_entry(rule%sets(1), 1, w(:,1))
    do i = 2, 6
      call set_entry(rule%sets(i), 6, w(:,i), [g(i-1)])
    end do
    call set_entry(rule%sets(7), 12, w(:,7), [g(6), g(6)])
    call set_entry(rule%sets(8), 12, w(:,8), [g(7), g(7)])
    call set_entry(rule%sets(9), 8, w(:,9), [g(8), g(8), g(8)])
    call set_entry(rule%sets(10), 8, w(:,10), [g(9), g(9), g(9)])
    call set_entry(rule%sets(11), 8, w(:,11), [g(10), g(10), g(10)])
    call set_entry(rule%sets(12), 24, w(:,12), [g(11), g(12), g(12)])
    call set_entry(rule%sets(13), 24, w(:,13), [g(13), g(13), g(14)])
    rule%errcoeff = [4.0_dp, 0.5_dp, 3.0_dp]
  end subroutine build_rule11

  subroutine build_rule9(d, rule)
    integer, intent(in) :: d
    type(cubature_rule_t), intent(out) :: rule
    real(dp), parameter :: w(42) = [ &
      -.0023611709677855118_dp,.11415390023857325_dp,-.6383392007670239_dp,.7484998850468521_dp, &
      -.0014324017033399125_dp,.057471507864489726_dp,-.14225104571434243_dp,-.06287502873828698_dp, &
      .2545911332489591_dp,-1.2073285666782363_dp,.8956736576416068_dp,-.36479356986049147_dp, &
      .0035417564516782677_dp,-.07260936739589368_dp,.10557491625218991_dp,.0021486025550098688_dp, &
      -.03226856389295395_dp,.010636783990231217_dp,.01468910249614349_dp,.5113470834646759_dp, &
      .45976448120806345_dp,.18239678493024573_dp,-.04508628929435784_dp,.21415883524352793_dp, &
      -.027351546526545645_dp,.054941067048711234_dp,.11937596202570775_dp,.6508951939192025_dp, &
      .1474493982943446_dp,.05769338449097348_dp,.034999626602143584_dp,-1.3868627719278281_dp, &
      -.2386668732575009_dp,.015532417276607053_dp,.003532809960709087_dp,.09231719987444222_dp, &
      .02254314464717892_dp,.013675773263272822_dp,-.32544759695960125_dp,.0017708782258391338_dp, &
      .0010743012775049344_dp,.2515001149531479_dp ]
    real(dp), parameter :: g(5) = [ .4779536579022695_dp,.20302858736911987_dp, &
      .44762735462617813_dp,.125_dp,.34303789878087815_dp ]
    integer :: twod, i
    real(dp) :: wt(5)

    twod = 2**d
    call allocate_sets(rule, 9, d)
    wt(1) = d*(d*(d*w(1) + w(2)) + w(3)) + w(4)
    wt(2) = d*(d*(d*w(5) + w(6)) + w(7)) - w(8)
    wt(3) = d*w(9) - wt(2)
    wt(4) = d*(d*w(10) + w(11)) - 1.0_dp + wt(1)
    wt(5) = d*w(12) + 1.0_dp - wt(1)
    call set_entry(rule%sets(1), 1, wt)

    wt(1) = d*(d*w(13) + w(14)) + w(15)
    wt(2) = d*(d*w(16) + w(17)) + w(18)
    wt(3) = w(19) - wt(2)
    wt(4) = d*w(20) + w(21) + wt(1)
    wt(5) = w(22) - wt(1)
    call set_entry(rule%sets(2), 2*d, wt, [g(1)])

    wt(1) = d*w(23) + w(24)
    wt(2) = d*w(25) + w(26)
    wt(3) = w(27) - wt(2)
    wt(4) = d*w(28) + w(29)
    wt(5) = -wt(1)
    call set_entry(rule%sets(3), 2*d, wt, [g(2)])

    wt = [w(30), w(31), -w(30), w(32), -w(30)]
    call set_entry(rule%sets(4), 2*d, wt, [g(3)])
    wt = [0.0_dp, 0.0_dp, w(33), 0.0_dp, 0.0_dp]
    call set_entry(rule%sets(5), 2*d, wt, [g(4)])

    wt = [w(34) - d*w(13), w(35) - d*w(16), 0.0_dp, 0.0_dp, 0.0_dp]
    wt(3) = -wt(2)
    wt(4) = w(36) + wt(1)
    wt(5) = -wt(1)
    call set_entry(rule%sets(6), 2*d*(d - 1), wt, [g(1), g(1)])

    wt = [w(37), w(38), -w(38), w(39), -w(37)]
    call set_entry(rule%sets(7), 4*d*(d - 1), wt, [g(1), g(2)])
    wt = [w(40), w(41), -w(41), w(40), -w(40)]
    call set_entry(rule%sets(8), 4*d*(d - 1)*(d - 2)/3, wt, [g(1), g(1), g(1)])
    wt = [w(42)/twod, w(8)/twod, -w(8)/twod, w(42)/twod, -w(42)/twod]
    call set_entry(rule%sets(9), twod, wt, [(g(5), i=1,d)])
    rule%errcoeff = [5.0_dp, 1.0_dp, 5.0_dp]
  end subroutine build_rule9

end module rcppnumerical_cuhre
