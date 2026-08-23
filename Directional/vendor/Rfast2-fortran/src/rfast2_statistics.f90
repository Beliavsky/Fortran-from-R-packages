module rfast2_statistics
   use rfast_special, only : dp, pi, normal_cdf, normal_pdf, normal_quantile, student_t_cdf, chisq_cdf
   use rfast_arrays, only : mean_r, variance_r, skewness_r, kurtosis_r, sort_real, colmeans, colvars
   use rfast_directional, only : vm_mle, circular_fit
   use rfast2_arrays, only : quantile_rfast2
   use rfast2_random, only : sample_int, runif
   use rfast2_types, only : scalar_test_result, km_result, meta_result, silhouette_result
   implicit none
   private

   public :: jarque_bera, jarque_bera_cols, empirical_entropy
   public :: cor_test_pearson, covar_cols, pooled_colvars
   public :: pinar1, col_pinar1, circular_cor, circular_cors
   public :: kaplan_meier, moran_i, wald_pois_ratio, walter_ci
   public :: wls_meta, ref_meta, silhouette_euclidean
   public :: permutation_ttest1, permutation_ttest2, bootstrap_ttest1
   public :: energy_equal_univariate

contains

   function jarque_bera(x) result(res)
      real(dp), intent(in) :: x(:)
      type(scalar_test_result) :: res
      real(dp) :: s, k

      s = skewness_r(x)
      k = kurtosis_r(x,.false.)
      res%statistic = real(size(x),dp)/6.0_dp * (s*s + 0.25_dp*(k-3.0_dp)**2)
      res%df = 2.0_dp
      res%pvalue = 1.0_dp - chisq_cdf(res%statistic,2.0_dp)
   end function jarque_bera

   function jarque_bera_cols(x) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: out(size(x,2),2)
      type(scalar_test_result) :: r
      integer :: j

      do j = 1, size(x,2)
         r = jarque_bera(x(:,j))
         out(j,1) = r%statistic
         out(j,2) = r%pvalue
      end do
   end function jarque_bera_cols

   real(dp) function empirical_entropy(x, k) result(h)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: k
      integer :: nb, i, bin, n
      integer, allocatable :: freq(:)
      real(dp) :: lo, hi, width, iqr, p

      n = size(x)
      lo = minval(x)
      hi = maxval(x)
      if (present(k)) then
         nb = max(1,k)
      else
         iqr = quantile_rfast2(x,0.75_dp) - quantile_rfast2(x,0.25_dp)
         width = 2.0_dp*iqr/real(max(1,n),dp)**(1.0_dp/3.0_dp)
         if (width <= tiny(1.0_dp)) then
            nb = max(1,int(sqrt(real(n,dp))))
         else
            nb = max(1,ceiling((hi-lo)/width))
         end if
      end if
      allocate(freq(nb))
      freq = 0
      if (hi <= lo) then
         h = 0.0_dp
         return
      end if
      do i = 1, n
         bin = 1 + int((x(i)-lo)/(hi-lo)*real(nb,dp))
         bin = max(1,min(nb,bin))
         freq(bin) = freq(bin) + 1
      end do
      h = 0.0_dp
      do i = 1, nb
         if (freq(i) > 0) then
            p = real(freq(i),dp)/real(n,dp)
            h = h - p*log(p)
         end if
      end do
   end function empirical_entropy

   function cor_test_pearson(y, x, rho, alpha) result(res)
      real(dp), intent(in) :: y(:), x(:)
      real(dp), intent(in), optional :: rho, alpha
      type(scalar_test_result) :: res
      real(dp) :: r0, a, mx, my, sx, sy, r, zh0, zh1, se, zcrit, e1, e2
      integer :: n

      n = size(y)
      r0 = 0.0_dp
      a = 0.05_dp
      if (present(rho)) r0 = rho
      if (present(alpha)) a = alpha
      mx = mean_r(x)
      my = mean_r(y)
      sx = sqrt(sum((x-mx)**2))
      sy = sqrt(sum((y-my)**2))
      r = sum((x-mx)*(y-my))/max(tiny(1.0_dp),sx*sy)
      r = max(-1.0_dp+1.0e-14_dp,min(1.0_dp-1.0e-14_dp,r))
      zh0 = 0.5_dp*log((1.0_dp+r0)/(1.0_dp-r0))
      zh1 = 0.5_dp*log((1.0_dp+r)/(1.0_dp-r))
      se = 1.0_dp/sqrt(real(n-3,dp))
      res%statistic = (zh1-zh0)/se
      res%df = real(n-3,dp)
      res%estimate = r
      res%pvalue = 2.0_dp*(1.0_dp-student_t_cdf(abs(res%statistic),res%df))
      zcrit = normal_quantile(1.0_dp-a/2.0_dp)
      e1 = exp(2.0_dp*(zh1-zcrit*se))
      e2 = exp(2.0_dp*(zh1+zcrit*se))
      res%lower = (e1-1.0_dp)/(e1+1.0_dp)
      res%upper = (e2-1.0_dp)/(e2+1.0_dp)
   end function cor_test_pearson

   function covar_cols(y, x) result(v)
      real(dp), intent(in) :: y(:), x(:,:)
      real(dp) :: v(size(x,2)), my
      integer :: j, n

      n = size(y)
      my = sum(y)
      do j = 1, size(x,2)
         v(j) = (dot_product(x(:,j),y)-sum(x(:,j))*my/real(n,dp))/real(n-1,dp)
      end do
   end function covar_cols

   function pooled_colvars(x, group, std) result(v)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: group(:)
      logical, intent(in), optional :: std
      real(dp) :: v(size(x,2)), m
      integer :: g, j, ng, ni, denom
      logical :: sd

      sd = .false.
      if (present(std)) sd = std
      ng = maxval(group)
      v = 0.0_dp
      denom = 0
      do g = 1, ng
         ni = count(group == g)
         if (ni <= 1) cycle
         denom = denom + ni - 1
         do j = 1, size(x,2)
            m = sum(pack(x(:,j),group==g))/real(ni,dp)
            v(j) = v(j) + sum((pack(x(:,j),group==g)-m)**2)
         end do
      end do
      if (denom > 0) v = v/real(denom,dp)
      if (sd) v = sqrt(max(0.0_dp,v))
   end function pooled_colvars

   function pinar1(x, unbiased) result(par)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: unbiased
      real(dp) :: par(2), sx1, sxn, num, den, alpha
      integer :: n
      logical :: ub

      n = size(x)
      sx1 = sum(x(2:n))
      sxn = sum(x(1:n-1))
      num = dot_product(x(2:n),x(1:n-1))-sx1*sxn/real(n-1,dp)
      den = sum(x(1:n-1)**2)-sxn*sxn/real(n-1,dp)
      alpha = num/den
      ub = .false.
      if (present(unbiased)) ub = unbiased
      if (ub .and. n > 3) alpha = (real(n,dp)*alpha+1.0_dp)/real(n-3,dp)
      par(1) = (sx1-alpha*sxn)/real(n-1,dp)
      par(2) = alpha
   end function pinar1

   function col_pinar1(x, unbiased) result(out)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: unbiased
      real(dp) :: out(size(x,2),2)
      integer :: j

      do j = 1, size(x,2)
         if (present(unbiased)) then
            out(j,:) = pinar1(x(:,j),unbiased)
         else
            out(j,:) = pinar1(x(:,j))
         end if
      end do
   end function col_pinar1

   function circular_cor(theta, phi, with_pvalue) result(out)
      real(dp), intent(in) :: theta(:), phi(:)
      logical, intent(in), optional :: with_pvalue
      real(dp) :: out(2), s1(size(theta)), s2(size(theta)), rho, lam22, lam02, lam20, z
      type(circular_fit) :: f1, f2
      logical :: pv
      integer :: n

      n = size(theta)
      f1 = vm_mle(theta)
      f2 = vm_mle(phi)
      s1 = sin(theta-f1%mu)
      s2 = sin(phi-f2%mu)
      rho = dot_product(s1,s2)/sqrt(max(tiny(1.0_dp),sum(s1*s1)*sum(s2*s2)))
      out = [rho,1.0_dp]
      pv = .false.
      if (present(with_pvalue)) pv = with_pvalue
      if (pv) then
         lam22 = sum(s1*s1*s2*s2)/real(n,dp)
         lam02 = sum(s2*s2)/real(n,dp)
         lam20 = sum(s1*s1)/real(n,dp)
         z = sqrt(real(n,dp))*sqrt(lam02*lam20/max(tiny(1.0_dp),lam22))*rho
         out(2) = 2.0_dp*(1.0_dp-normal_cdf(abs(z)))
      end if
   end function circular_cor

   function circular_cors(theta, phi, with_pvalue) result(out)
      real(dp), intent(in) :: theta(:), phi(:,:)
      logical, intent(in), optional :: with_pvalue
      real(dp) :: out(size(phi,2),2)
      integer :: j

      do j = 1, size(phi,2)
         if (present(with_pvalue)) then
            out(j,:) = circular_cor(theta,phi(:,j),with_pvalue)
         else
            out(j,:) = circular_cor(theta,phi(:,j))
         end if
      end do
   end function circular_cors

   function kaplan_meier(time, status) result(res)
      real(dp), intent(in) :: time(:)
      integer, intent(in) :: status(:)
      type(km_result) :: res
      integer, allocatable :: ord(:)
      real(dp), allocatable :: ts(:)
      integer, allocatable :: ds(:)
      integer :: n, i, j, k, at_risk, ev
      real(dp) :: surv

      n = size(time)
      allocate(ord(n))
      ord = [(i,i=1,n)]
      call sort_index(time,ord)
      allocate(ts(n),ds(n))
      ts = time(ord)
      ds = status(ord)
      allocate(res%time(max(1,n)),res%risk(max(1,n)),res%events(max(1,n)),res%survival(max(1,n)))
      k = 0
      surv = 1.0_dp
      i = 1
      do while (i <= n)
         j = i
         ev = 0
         do while (j <= n)
            if (abs(ts(j)-ts(i)) > epsilon(1.0_dp)*max(1.0_dp,abs(ts(i)))) exit
            if (ds(j) > 0) ev = ev + 1
            j = j + 1
         end do
         if (ev > 0) then
            k = k + 1
            at_risk = n - i + 1
            surv = surv*(1.0_dp-real(ev,dp)/real(at_risk,dp))
            res%time(k) = ts(i)
            res%risk(k) = at_risk
            res%events(k) = ev
            res%survival(k) = surv
         end if
         i = j
      end do
      if (k < size(res%time)) call shrink_km(res,k)
   end function kaplan_meier

   subroutine sort_index(x, idx)
      real(dp), intent(in) :: x(:)
      integer, intent(inout) :: idx(:)
      integer :: i, j, key

      do i = 2, size(idx)
         key = idx(i)
         j = i-1
         do while (j >= 1)
            if (x(idx(j)) <= x(key)) exit
            idx(j+1) = idx(j)
            j = j-1
         end do
         idx(j+1) = key
      end do
   end subroutine sort_index

   subroutine shrink_km(res,k)
      type(km_result), intent(inout) :: res
      integer, intent(in) :: k
      real(dp), allocatable :: rt(:), rs(:)
      integer, allocatable :: rr(:), re(:)

      allocate(rt(k),rs(k),rr(k),re(k))
      if (k > 0) then
         rt = res%time(1:k)
         rs = res%survival(1:k)
         rr = res%risk(1:k)
         re = res%events(1:k)
      end if
      call move_alloc(rt,res%time)
      call move_alloc(rs,res%survival)
      call move_alloc(rr,res%risk)
      call move_alloc(re,res%events)
   end subroutine shrink_km

   function moran_i(x, w, permutations) result(out)
      real(dp), intent(in) :: x(:), w(:,:)
      integer, intent(in), optional :: permutations
      real(dp) :: out(2), wn(size(w,1),size(w,2)), y(size(x)), xp(size(x)), mx, my, sx, sy, r, obs, rp
      integer, allocatable :: idx(:)
      integer :: i, rnum, exceed

      wn = w
      do i = 1, size(w,1)
         if (sum(abs(wn(i,:))) > tiny(1.0_dp)) wn(i,:) = wn(i,:)/sum(wn(i,:))
      end do
      y = matmul(wn,x)
      mx = mean_r(x)
      my = mean_r(y)
      sx = sqrt(sum((x-mx)**2))
      sy = sqrt(sum((y-my)**2))
      r = sum((x-mx)*(y-my))/max(tiny(1.0_dp),sx*sy)
      obs = r*sqrt(variance_r(y)/max(tiny(1.0_dp),variance_r(x)))
      out = [obs,1.0_dp]
      rnum = 0
      if (present(permutations)) rnum = permutations
      if (rnum > 1) then
         exceed = 0
         do i = 1, rnum
            idx = sample_int(size(x),size(x),.false.)
            xp = x(idx)
            mx = mean_r(xp)
            sx = sqrt(sum((xp-mx)**2))
            rp = sum((xp-mx)*(y-my))/max(tiny(1.0_dp),sx*sy)
            if (abs(rp) >= abs(r)) exceed = exceed + 1
         end do
         out(2) = real(exceed+1,dp)/real(rnum+1,dp)
      end if
   end function moran_i

   function wald_pois_ratio(x, y, alpha) result(res)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), intent(in), optional :: alpha
      type(scalar_test_result) :: res
      real(dp) :: a, lam1, lam2, varat, z

      a = 0.05_dp
      if (present(alpha)) a = alpha
      lam1 = mean_r(x)
      lam2 = mean_r(y)
      res%estimate = lam1/lam2
      varat = lam1/real(size(x),dp)/(lam2*lam2) + lam2/real(size(y),dp)*lam1*lam1/(lam2**4)
      z = normal_quantile(1.0_dp-a/2.0_dp)
      res%lower = res%estimate-z*sqrt(varat)
      res%upper = res%estimate+z*sqrt(varat)
   end function wald_pois_ratio

   function walter_ci(x1,x2,n1,n2,alpha) result(res)
      real(dp), intent(in) :: x1,x2,n1,n2
      real(dp), intent(in), optional :: alpha
      type(scalar_test_result) :: res
      real(dp) :: a,x,y,n,m,se,lr,z

      a = 0.05_dp
      if (present(alpha)) a = alpha
      x = x1+0.5_dp
      y = x2+0.5_dp
      n = n1+0.5_dp
      m = n2+0.5_dp
      se = sqrt(1.0_dp/x+1.0_dp/y-1.0_dp/m-1.0_dp/n)
      lr = log(x/n)-log(y/m)
      z = normal_quantile(1.0_dp-a/2.0_dp)
      res%estimate = exp(lr)
      res%lower = exp(lr-z*se)
      res%upper = exp(lr+z*se)
   end function walter_ci

   function wls_meta(yi,vi) result(out)
      real(dp), intent(in) :: yi(:),vi(:)
      real(dp) :: out(7),w(size(yi)),sw,fe,phi,h,se,tcrit,pv
      integer :: m

      m = size(yi)
      w = 1.0_dp/vi
      sw = sum(w)
      fe = dot_product(yi,w)/sw
      phi = sum((yi-fe)**2/vi)/real(max(1,m-1),dp)
      h = (phi-1.0_dp)/phi
      se = sqrt(phi/sw)
      tcrit = normal_quantile(0.975_dp)
      pv = 2.0_dp*(1.0_dp-student_t_cdf(abs(fe)/se,real(max(1,m-2),dp)))
      out = [fe,se,fe-tcrit*se,fe+tcrit*se,pv,phi,h]
   end function wls_meta

   function ref_meta(yi,vi,tol) result(res)
      real(dp), intent(in) :: yi(:),vi(:)
      real(dp), intent(in), optional :: tol
      type(meta_result) :: res
      real(dp) :: a,b,c,d,fc,fd,eps,tau,w(size(yi)),sw,mfe,v
      integer :: it

      eps = 1.0e-7_dp
      if (present(tol)) eps = tol
      a = 0.0_dp
      b = 100.0_dp
      c = b-(b-a)/1.6180339887498948482_dp
      d = a+(b-a)/1.6180339887498948482_dp
      fc = reml_objective(c,yi,vi)
      fd = reml_objective(d,yi,vi)
      do it = 1, 300
         if (abs(b-a) <= eps) exit
         if (fc < fd) then
            b = d
            d = c
            fd = fc
            c = b-(b-a)/1.6180339887498948482_dp
            fc = reml_objective(c,yi,vi)
         else
            a = c
            c = d
            fc = fd
            d = a+(b-a)/1.6180339887498948482_dp
            fd = reml_objective(d,yi,vi)
         end if
      end do
      tau = 0.5_dp*(a+b)
      w = 1.0_dp/vi
      sw = sum(w)
      mfe = dot_product(yi,w)/sw
      v = 9.0_dp*sw/(sw*sw-sum(w*w))
      res%fixed_mean = mfe
      res%v = v
      res%i2 = tau/(tau+v)
      res%h2 = (tau+v)/v
      res%q = sum((yi-mfe)**2*w)
      res%pvalue = 1.0_dp-chisq_cdf(res%q,real(size(yi)-1,dp))
      res%tau2 = tau
      w = 1.0_dp/(vi+tau)
      res%random_mean = dot_product(yi,w)/sum(w)
   end function ref_meta

   real(dp) function reml_objective(tau,yi,vi) result(v)
      real(dp), intent(in) :: tau,yi(:),vi(:)
      real(dp) :: w(size(yi)),m

      w = 1.0_dp/(vi+tau)
      m = dot_product(w,yi)/sum(w)
      v = -sum(log(w))+log(sum(w))+sum((yi-m)**2*w)
   end function reml_objective

   function silhouette_euclidean(x,cl) result(res)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: cl(:)
      type(silhouette_result) :: res
      integer :: n,g,i,j,c,ni
      real(dp) :: a,b,d,s
      logical :: found_other

      n = size(x,1)
      g = maxval(cl)
      allocate(res%value(n),res%stats(g,4))
      res%value = 0.0_dp
      res%stats = 0.0_dp
      do i = 1, n
         ni = count(cl == cl(i))
         a = 0.0_dp
         if (ni > 1) then
            do j = 1, n
               if (j /= i .and. cl(j) == cl(i)) a = a + sqrt(sum((x(i,:)-x(j,:))**2))
            end do
            a = a/real(ni-1,dp)
         end if
         b = huge(1.0_dp)
         found_other = .false.
         do c = 1, g
            if (c == cl(i)) cycle
            ni = count(cl == c)
            if (ni == 0) cycle
            d = 0.0_dp
            do j = 1, n
               if (cl(j) == c) d = d + sqrt(sum((x(i,:)-x(j,:))**2))
            end do
            b = min(b,d/real(ni,dp))
            found_other = .true.
         end do
         if (.not. found_other) then
            res%value(i) = 0.0_dp
         else
            res%value(i) = (b-a)/max(a,b)
         end if
      end do
      do c = 1, g
         ni = count(cl == c)
         res%stats(c,1) = real(ni,dp)
         if (ni > 0) then
            s = sum(pack(res%value,cl==c))/real(ni,dp)
            res%stats(c,2) = minval(pack(res%value,cl==c))
            res%stats(c,3) = maxval(pack(res%value,cl==c))
            res%stats(c,4) = s
         end if
      end do
   end function silhouette_euclidean

   function permutation_ttest1(x,m,b) result(res)
      real(dp), intent(in) :: x(:),m
      integer, intent(in), optional :: b
      type(scalar_test_result) :: res
      real(dp) :: z(size(x)),stat,pstat
      integer :: r,i,exceed
      real(dp), allocatable :: signs(:)

      r = 999
      if (present(b)) r = b
      z = x-m
      stat = abs(sum(z))
      exceed = 0
      do i = 1, r
         signs = sign_draws(size(x))
         pstat = abs(dot_product(z,signs))
         if (pstat >= stat) exceed = exceed+1
      end do
      res%statistic = stat
      res%pvalue = real(exceed+1,dp)/real(r+1,dp)
   end function permutation_ttest1

   function sign_draws(n) result(s)
      integer, intent(in) :: n
      real(dp), allocatable :: s(:),u(:)
      integer :: i

      u = runif(n)
      allocate(s(n))
      do i = 1, n
         if (u(i) < 0.5_dp) then
            s(i) = -1.0_dp
         else
            s(i) = 1.0_dp
         end if
      end do
   end function sign_draws

   function permutation_ttest2(x,y,b) result(res)
      real(dp), intent(in) :: x(:),y(:)
      integer, intent(in), optional :: b
      type(scalar_test_result) :: res
      real(dp), allocatable :: z(:)
      integer, allocatable :: idx(:)
      integer :: r,i,nx,ny,n,exceed
      real(dp) :: stat,pstat

      nx = size(x)
      ny = size(y)
      n = nx+ny
      r = 999
      if (present(b)) r = b
      allocate(z(n))
      z(1:nx) = x
      z(nx+1:n) = y
      stat = abs(mean_r(x)-mean_r(y))
      exceed = 0
      do i = 1, r
         idx = sample_int(n,n,.false.)
         pstat = abs(sum(z(idx(1:nx)))/real(nx,dp)-sum(z(idx(nx+1:n)))/real(ny,dp))
         if (pstat >= stat) exceed = exceed+1
      end do
      res%statistic = stat
      res%pvalue = real(exceed+1,dp)/real(r+1,dp)
   end function permutation_ttest2

   function bootstrap_ttest1(x,m,b) result(res)
      real(dp), intent(in) :: x(:),m
      integer, intent(in), optional :: b
      type(scalar_test_result) :: res
      real(dp) :: z(size(x)), xb(size(x)), stat, bs, bm
      integer, allocatable :: idx(:)
      integer :: r,i,exceed,n

      n = size(x)
      r = 999
      if (present(b)) r = b
      stat = sqrt(real(n,dp))*(mean_r(x)-m)/sqrt(variance_r(x))
      z = x-mean_r(x)+m
      exceed = 0
      do i = 1, r
         idx = sample_int(n,n,.true.)
         xb = z(idx)
         bm = mean_r(xb)
         bs = sqrt(variance_r(xb))
         if (bs > tiny(1.0_dp)) then
            if (abs(sqrt(real(n,dp))*(bm-m)/bs) >= abs(stat)) exceed = exceed+1
         end if
      end do
      res%statistic = stat
      res%pvalue = real(exceed+1,dp)/real(r+1,dp)
   end function bootstrap_ttest1

   real(dp) function energy_equal_univariate(x,y,permutations) result(pvalue)
      real(dp), intent(in) :: x(:),y(:)
      integer, intent(in), optional :: permutations
      real(dp), allocatable :: z(:),xp(:),yp(:)
      integer, allocatable :: idx(:)
      integer :: nx,ny,n,r,k,exceed
      real(dp) :: stat,pstat

      nx = size(x)
      ny = size(y)
      n = nx+ny
      r = 999
      if (present(permutations)) r = permutations
      stat = energy_stat(x,y)
      allocate(z(n),xp(nx),yp(ny))
      z(1:nx) = x
      z(nx+1:n) = y
      exceed = 0
      do k = 1, r
         idx = sample_int(n,n,.false.)
         xp = z(idx(1:nx))
         yp = z(idx(nx+1:n))
         pstat = energy_stat(xp,yp)
         if (pstat >= stat) exceed = exceed+1
      end do
      pvalue = real(exceed+1,dp)/real(r+1,dp)
   end function energy_equal_univariate

   real(dp) function energy_stat(x,y) result(stat)
      real(dp), intent(in) :: x(:),y(:)
      real(dp) :: xy,xx,yy
      integer :: i,j,nx,ny

      nx = size(x)
      ny = size(y)
      xy = 0.0_dp
      xx = 0.0_dp
      yy = 0.0_dp
      do i = 1, nx
         do j = 1, ny
            xy = xy+abs(x(i)-y(j))
         end do
      end do
      do i = 1, nx
         do j = 1, nx
            xx = xx+abs(x(i)-x(j))
         end do
      end do
      do i = 1, ny
         do j = 1, ny
            yy = yy+abs(y(i)-y(j))
         end do
      end do
      stat = 2.0_dp*xy/real(nx*ny,dp)-xx/real(nx*nx,dp)-yy/real(ny*ny,dp)
      stat = real(nx*ny,dp)/real(nx+ny,dp)*stat
   end function energy_stat

end module rfast2_statistics
