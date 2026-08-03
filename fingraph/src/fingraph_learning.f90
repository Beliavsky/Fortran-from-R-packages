! SPDX-License-Identifier: GPL-3.0-only
module fingraph_learning
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use fingraph_kinds, only : dp
   use fingraph_status, only : fg_ok, fg_invalid_input, fg_size_mismatch, &
      fg_no_convergence
   use fingraph_types, only : fingraph_result
   use fingraph_operators, only : L, A, Lstar, Dstar, Linv, Ainv, vec, vecLmat
   use fingraph_linalg, only : symmetric_pseudoinverse, symmetric_eigen_jacobi, &
      frobenius_norm, vector_norm2, identity_matrix, nnqp_projected_gradient, &
      outer_product, rows_correlation
   implicit none
   private

   public :: learn_connected_graph
   public :: learn_regular_heavytail_graph
   public :: learn_kcomp_heavytail_graph
   public :: compute_student_weight
   public :: compute_augmented_lagrangian_ht
   public :: compute_augmented_lagrangian_kcomp_ht

contains

   subroutine learn_connected_graph(s, result, initial_weights, initialization, &
                                    degree, degrees, rho, maxiter, reltol)
      real(dp), intent(in) :: s(:,:)
      type(fingraph_result), intent(out) :: result
      real(dp), intent(in), optional :: initial_weights(:)
      character(len=*), intent(in), optional :: initialization
      real(dp), intent(in), optional :: degree
      real(dp), intent(in), optional :: degrees(:)
      real(dp), intent(in), optional :: rho, reltol
      integer, intent(in), optional :: maxiter

      real(dp), allocatable :: w(:), wi(:), lw(:,:), lwi(:,:), theta(:,:), thetai(:,:)
      real(dp), allocatable :: ymat(:,:), yvec(:), dvec(:), jmat(:,:)
      real(dp), allocatable :: grad(:), r1(:,:), r2(:), lstar_s(:)
      real(dp) :: rho_now, tol, eta, mu, tau, rnorm, snorm, relative_change
      integer :: p, m, iter, max_it, status
      logical :: converged

      call clear_result(result)
      p = size(s,1)
      if (p < 2 .or. size(s,2) /= p .or. .not. all(ieee_is_finite(s))) then
         result%status = fg_invalid_input
         return
      end if
      m = p*(p-1)/2
      call build_degree_vector(p,degree,degrees,dvec,status)
      if (status /= fg_ok) then
         result%status = status
         return
      end if
      call initialize_weights(s,initial_weights,initialization,w,status)
      if (status /= fg_ok .or. size(w) /= m) then
         result%status = merge(status,fg_size_mismatch,status/=fg_ok)
         return
      end if
      call normalize_initial_weights(w,status)
      if (status /= fg_ok) then
         result%status = status
         return
      end if

      rho_now = 1.0_dp
      if (present(rho)) rho_now = rho
      tol = 1.0e-5_dp
      if (present(reltol)) tol = reltol
      max_it = 10000
      if (present(maxiter)) max_it = maxiter
      if (rho_now <= 0.0_dp .or. tol <= 0.0_dp .or. max_it < 1) then
         result%status = fg_invalid_input
         return
      end if

      allocate(jmat(p,p),ymat(p,p),yvec(p),wi(m),grad(m),r1(p,p),r2(p))
      jmat = 1.0_dp/real(p,dp)
      lw = L(w)
      theta = lw
      ymat = 0.0_dp
      yvec = 0.0_dp
      lstar_s = Lstar(s)
      mu = 2.0_dp
      tau = 2.0_dp
      converged = .false.
      thetai = theta
      wi = w

      do iter = 1,max_it
         grad = lstar_s - Lstar(ymat + rho_now*theta) &
              + Dstar(yvec-rho_now*dvec) &
              + rho_now*(Lstar(lw)+Dstar(diagonal(lw)))
         eta = 1.0_dp/(2.0_dp*rho_now*real(2*p-1,dp))
         wi = max(w-eta*grad,0.0_dp)
         lwi = L(wi)
         call update_theta_connected(lwi,jmat,ymat,rho_now,thetai,status)
         if (status /= fg_ok) then
            result%status = status
            return
         end if
         r1 = thetai-lwi
         ymat = ymat+rho_now*r1
         r2 = diagonal(lwi)-dvec
         yvec = yvec+rho_now*r2

         snorm = rho_now*vector_norm2(Lstar(theta-thetai))
         rnorm = frobenius_norm(r1)
         if (rnorm > mu*snorm) then
            rho_now = rho_now*tau
         else if (snorm > mu*rnorm) then
            rho_now = rho_now/tau
         end if

         relative_change = frobenius_norm(lwi-lw)/max(frobenius_norm(lw),tiny(1.0_dp))
         converged = iter > 1 .and. relative_change < tol
         if (converged) exit
         w = wi
         lw = lwi
         theta = thetai
      end do

      iter = min(iter,max_it)
      result%weights = wi
      result%laplacian = L(wi)
      result%adjacency = A(wi)
      result%theta = thetai
      result%iterations = iter
      result%convergence = converged
      result%rho = rho_now
      result%status = merge(fg_ok,fg_no_convergence,converged)
   end subroutine learn_connected_graph

   subroutine learn_regular_heavytail_graph(x, result, heavy_type, nu, initial_weights, &
                                            initialization, degree, degrees, rho, &
                                            update_rho, maxiter, reltol)
      real(dp), intent(in) :: x(:,:)
      type(fingraph_result), intent(out) :: result
      character(len=*), intent(in), optional :: heavy_type, initialization
      real(dp), intent(in), optional :: nu, initial_weights(:), degree, degrees(:)
      real(dp), intent(in), optional :: rho, reltol
      logical, intent(in), optional :: update_rho
      integer, intent(in), optional :: maxiter

      real(dp), allocatable :: xwork(:,:), cor(:,:), lstar_sq(:,:), w(:), wi(:)
      real(dp), allocatable :: lw(:,:), lwi(:,:), theta(:,:), thetai(:,:), jmat(:,:)
      real(dp), allocatable :: ymat(:,:), yvec(:), dvec(:), grad(:), weighted(:)
      real(dp), allocatable :: r1(:,:), r2(:), primal_lap(:), primal_deg(:)
      real(dp), allocatable :: dual(:), lagr(:), elapsed(:)
      real(dp) :: rho_now, tol, nu_now, eta, mu, tau, rnorm, snorm
      real(dp) :: relative_change, start_time, now
      integer :: n, p, m, iter, max_it, status
      logical :: converged, adapt_rho, student
      character(len=:), allocatable :: distribution

      call clear_result(result)
      n = size(x,1)
      p = size(x,2)
      if (n < 2 .or. p < 2 .or. .not. all(ieee_is_finite(x))) then
         result%status = fg_invalid_input
         return
      end if
      xwork = x
      distribution = "gaussian"
      if (present(heavy_type)) distribution = lowercase(trim(heavy_type))
      student = distribution == "student"
      if (.not. student .and. distribution /= "gaussian") then
         result%status = fg_invalid_input
         return
      end if
      nu_now = huge(1.0_dp)
      if (present(nu)) nu_now = nu
      if (student .and. (.not. present(nu) .or. nu_now <= 2.0_dp)) then
         result%status = fg_invalid_input
         return
      end if

      m = p*(p-1)/2
      call build_degree_vector(p,degree,degrees,dvec,status)
      if (status /= fg_ok) then
         result%status = status
         return
      end if
      call column_correlation(xwork,cor,status)
      if (status /= fg_ok) then
         result%status = status
         return
      end if
      call initialize_weights(cor,initial_weights,initialization,w,status)
      if (status /= fg_ok .or. size(w) /= m) then
         result%status = merge(status,fg_size_mismatch,status/=fg_ok)
         return
      end if
      call normalize_initial_weights(w,status)
      if (status /= fg_ok) then
         result%status = status
         return
      end if
      call build_lstar_scatter(xwork,real(n-1,dp),lstar_sq)

      rho_now = 1.0_dp
      if (present(rho)) rho_now = rho
      tol = 1.0e-5_dp
      if (present(reltol)) tol = reltol
      max_it = 10000
      if (present(maxiter)) max_it = maxiter
      adapt_rho = .true.
      if (present(update_rho)) adapt_rho = update_rho
      if (rho_now <= 0.0_dp .or. tol <= 0.0_dp .or. max_it < 1) then
         result%status = fg_invalid_input
         return
      end if

      allocate(jmat(p,p),ymat(p,p),yvec(p),wi(m),grad(m),weighted(m))
      allocate(r1(p,p),r2(p),primal_lap(max_it),primal_deg(max_it))
      allocate(dual(max_it),lagr(max_it),elapsed(max_it))
      jmat = 1.0_dp/real(p,dp)
      lw = L(w)
      theta = lw
      ymat = 0.0_dp
      yvec = 0.0_dp
      mu = 2.0_dp
      tau = 2.0_dp
      converged = .false.
      thetai = theta
      wi = w
      call cpu_time(start_time)

      do iter = 1,max_it
         call weighted_scatter_adjoint(w,lstar_sq,p,nu_now,student,weighted)
         grad = weighted - Lstar(rho_now*theta+ymat) &
              + Dstar(yvec-rho_now*dvec) &
              + rho_now*(Lstar(lw)+Dstar(diagonal(lw)))
         eta = 1.0_dp/(2.0_dp*rho_now*real(2*p-1,dp))
         wi = max(w-eta*grad,0.0_dp)
         lwi = L(wi)
         call update_theta_connected(lwi,jmat,ymat,rho_now,thetai,status)
         if (status /= fg_ok) then
            result%status = status
            return
         end if
         r1 = thetai-lwi
         ymat = ymat+rho_now*r1
         r2 = diagonal(lwi)-dvec
         yvec = yvec+rho_now*r2

         primal_lap(iter) = frobenius_norm(r1)
         primal_deg(iter) = vector_norm2(r2)
         dual(iter) = rho_now*vector_norm2(Lstar(theta-thetai))
         lagr(iter) = compute_augmented_lagrangian_ht(wi,lstar_sq,thetai,jmat, &
              ymat,yvec,dvec,student,nu_now,rho_now)

         if (adapt_rho) then
            snorm = dual(iter)
            rnorm = primal_lap(iter)
            if (rnorm > mu*snorm) then
               rho_now = rho_now*tau
            else if (snorm > mu*rnorm) then
               rho_now = rho_now/tau
            end if
         end if

         relative_change = frobenius_norm(lwi-lw)/max(frobenius_norm(lw),tiny(1.0_dp))
         converged = iter > 1 .and. relative_change < tol
         call cpu_time(now)
         elapsed(iter) = now-start_time
         if (converged) exit
         w = wi
         lw = lwi
         theta = thetai
      end do

      iter = min(iter,max_it)
      result%weights = wi
      result%laplacian = L(wi)
      result%adjacency = A(wi)
      result%theta = thetai
      result%primal_lap_residual = primal_lap(:iter)
      result%primal_deg_residual = primal_deg(:iter)
      result%dual_residual = dual(:iter)
      result%lagrangian = lagr(:iter)
      result%elapsed_time = elapsed(:iter)
      result%iterations = iter
      result%convergence = converged
      result%rho = rho_now
      result%status = merge(fg_ok,fg_no_convergence,converged)
   end subroutine learn_regular_heavytail_graph

   subroutine learn_kcomp_heavytail_graph(x, result, k, heavy_type, nu, initial_weights, &
                                          initialization, degree, degrees, beta, &
                                          update_beta, early_stopping, rho, update_rho, &
                                          maxiter, reltol, record_objective)
      real(dp), intent(in) :: x(:,:)
      type(fingraph_result), intent(out) :: result
      integer, intent(in), optional :: k, maxiter
      character(len=*), intent(in), optional :: heavy_type, initialization
      real(dp), intent(in), optional :: nu, initial_weights(:), degree, degrees(:)
      real(dp), intent(in), optional :: beta, rho, reltol
      logical, intent(in), optional :: update_beta, early_stopping, update_rho
      logical, intent(in), optional :: record_objective

      real(dp), allocatable :: xs(:,:), cor(:,:), lstar_sq(:,:), w(:), wi(:)
      real(dp), allocatable :: lw(:,:), lwi(:,:), theta(:,:), thetai(:,:)
      real(dp), allocatable :: u(:,:), ymat(:,:), yvec(:), dvec(:), grad(:), weighted(:)
      real(dp), allocatable :: r1(:,:), r2(:), primal_lap(:), primal_deg(:)
      real(dp), allocatable :: dual(:), lagr(:), elapsed(:), beta_hist(:)
      real(dp), allocatable :: values(:), vectors(:,:), uut(:,:)
      real(dp) :: rho_now, beta_now, tol, nu_now, eta, mu, tau, rnorm, snorm
      real(dp) :: relative_change, start_time, now
      integer :: n, p, m, k_now, iter, max_it, status, nzero
      logical :: converged, adapt_rho, adapt_beta, stop_early, student, keep_objective
      character(len=:), allocatable :: distribution

      call clear_result(result)
      n = size(x,1)
      p = size(x,2)
      k_now = 1
      if (present(k)) k_now = k
      if (n < 2 .or. p < 2 .or. k_now < 1 .or. k_now >= p .or. &
          .not. all(ieee_is_finite(x))) then
         result%status = fg_invalid_input
         return
      end if
      call standardize_columns(x,xs,status)
      if (status /= fg_ok) then
         result%status = status
         return
      end if
      distribution = "gaussian"
      if (present(heavy_type)) distribution = lowercase(trim(heavy_type))
      student = distribution == "student"
      if (.not. student .and. distribution /= "gaussian") then
         result%status = fg_invalid_input
         return
      end if
      nu_now = huge(1.0_dp)
      if (present(nu)) nu_now = nu
      if (student .and. (.not. present(nu) .or. nu_now <= 2.0_dp)) then
         result%status = fg_invalid_input
         return
      end if

      m = p*(p-1)/2
      call build_degree_vector(p,degree,degrees,dvec,status)
      if (status /= fg_ok) then
         result%status = status
         return
      end if
      call column_correlation(xs,cor,status)
      if (status /= fg_ok) then
         result%status = status
         return
      end if
      call initialize_weights(cor,initial_weights,initialization,w,status)
      if (status /= fg_ok .or. size(w) /= m) then
         result%status = merge(status,fg_size_mismatch,status/=fg_ok)
         return
      end if
      call normalize_initial_weights(w,status)
      if (status /= fg_ok) then
         result%status = status
         return
      end if
      call build_lstar_scatter(xs,real(n,dp),lstar_sq)

      rho_now = 1.0_dp
      if (present(rho)) rho_now = rho
      beta_now = 1.0e-8_dp
      if (present(beta)) beta_now = beta
      tol = 1.0e-5_dp
      if (present(reltol)) tol = reltol
      max_it = 10000
      if (present(maxiter)) max_it = maxiter
      adapt_rho = .false.
      if (present(update_rho)) adapt_rho = update_rho
      adapt_beta = .true.
      if (present(update_beta)) adapt_beta = update_beta
      stop_early = .false.
      if (present(early_stopping)) stop_early = early_stopping
      keep_objective = .false.
      if (present(record_objective)) keep_objective = record_objective
      if (rho_now <= 0.0_dp .or. beta_now < 0.0_dp .or. tol <= 0.0_dp .or. max_it < 1) then
         result%status = fg_invalid_input
         return
      end if

      allocate(ymat(p,p),yvec(p),wi(m),grad(m),weighted(m),r1(p,p),r2(p))
      allocate(primal_lap(max_it),primal_deg(max_it),dual(max_it),lagr(max_it))
      allocate(elapsed(max_it),beta_hist(max_it))
      lw = L(w)
      theta = lw
      call symmetric_eigen_jacobi(lw,values,vectors,status)
      if (status /= fg_ok) then
         result%status = status
         return
      end if
      u = vectors(:,1:k_now)
      ymat = 0.0_dp
      yvec = 0.0_dp
      mu = 2.0_dp
      tau = 2.0_dp
      converged = .false.
      thetai = theta
      wi = w
      beta_hist = 0.0_dp
      call cpu_time(start_time)

      do iter = 1,max_it
         call weighted_scatter_adjoint(w,lstar_sq,p,nu_now,student,weighted)
         uut = matmul(u,transpose(u))
         grad = weighted + Lstar(beta_now*uut-ymat-rho_now*theta) &
              + Dstar(yvec-rho_now*dvec) &
              + rho_now*(Lstar(lw)+Dstar(diagonal(lw)))
         eta = 1.0_dp/(2.0_dp*rho_now*real(2*p-1,dp))
         wi = max(w-eta*grad,0.0_dp)
         lwi = L(wi)

         call symmetric_eigen_jacobi(lwi,values,vectors,status)
         if (status /= fg_ok) then
            result%status = status
            return
         end if
         u = vectors(:,1:k_now)
         call update_theta_kcomp(lwi,ymat,rho_now,k_now,thetai,status)
         if (status /= fg_ok) then
            result%status = status
            return
         end if
         r1 = thetai-lwi
         ymat = ymat+rho_now*r1
         r2 = diagonal(lwi)-dvec
         yvec = yvec+rho_now*r2

         primal_lap(iter) = frobenius_norm(r1)
         primal_deg(iter) = vector_norm2(r2)
         dual(iter) = rho_now*vector_norm2(Lstar(theta-thetai))
         if (keep_objective) then
            lagr(iter) = compute_augmented_lagrangian_kcomp_ht(wi,lstar_sq,thetai,u, &
                 ymat,yvec,dvec,student,nu_now,rho_now,beta_now,k_now)
         else
            lagr(iter) = 0.0_dp
         end if

         if (adapt_rho) then
            snorm = dual(iter)
            rnorm = primal_lap(iter)
            if (rnorm > mu*snorm) then
               rho_now = rho_now*tau
            else if (snorm > mu*rnorm) then
               rho_now = rho_now/tau
            end if
         end if

         if (adapt_beta) then
            call symmetric_eigen_jacobi(lwi,values,vectors,status)
            if (status /= fg_ok) then
               result%status = status
               return
            end if
            nzero = count(values < 1.0e-9_dp)
            if (k_now < nzero) then
               beta_now = 0.5_dp*beta_now
            else if (k_now > nzero) then
               beta_now = 2.0_dp*beta_now
            else if (stop_early) then
               converged = .true.
            end if
            beta_hist(iter) = beta_now
         end if

         relative_change = frobenius_norm(lwi-lw)/max(frobenius_norm(lw),tiny(1.0_dp))
         if (.not. converged) converged = iter > 1 .and. relative_change < tol
         call cpu_time(now)
         elapsed(iter) = now-start_time
         if (converged) exit
         w = wi
         lw = lwi
         theta = thetai
      end do

      iter = min(iter,max_it)
      result%weights = wi
      result%laplacian = L(wi)
      result%adjacency = A(wi)
      result%theta = thetai
      result%primal_lap_residual = primal_lap(:iter)
      result%primal_deg_residual = primal_deg(:iter)
      result%dual_residual = dual(:iter)
      if (keep_objective) result%lagrangian = lagr(:iter)
      result%elapsed_time = elapsed(:iter)
      if (adapt_beta) result%beta_seq = beta_hist(:iter)
      result%iterations = iter
      result%convergence = converged
      result%rho = rho_now
      result%beta = beta_now
      result%status = merge(fg_ok,fg_no_convergence,converged)
   end subroutine learn_kcomp_heavytail_graph

   pure function compute_student_weight(w,lstar_sq,p,nu) result(weight)
      real(dp), intent(in) :: w(:), lstar_sq(:), nu
      integer, intent(in) :: p
      real(dp) :: weight
      weight = (real(p,dp)+nu)/(dot_product(w,lstar_sq)+nu)
   end function compute_student_weight

   function compute_augmented_lagrangian_ht(w,lstar_sq,theta,jmat,ymat,yvec,dvec, &
                                             student,nu,rho) result(value)
      real(dp), intent(in) :: w(:), lstar_sq(:,:), theta(:,:), jmat(:,:)
      real(dp), intent(in) :: ymat(:,:), yvec(:), dvec(:), nu, rho
      logical, intent(in) :: student
      real(dp) :: value
      real(dp), allocatable :: values(:), vectors(:,:), lw(:,:), dw(:), diff(:,:)
      real(dp) :: ufunc, term
      integer :: q, n, p, status
      n = size(lstar_sq,1)
      p = size(theta,1)
      call symmetric_eigen_jacobi(theta+jmat,values,vectors,status)
      if (status /= fg_ok .or. minval(values) <= 0.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      lw = L(w)
      dw = diagonal(lw)
      ufunc = 0.0_dp
      do q = 1,n
         if (student) then
            term = 1.0_dp+real(n,dp)*dot_product(w,lstar_sq(q,:))/nu
            if (term <= 0.0_dp) then
               value = huge(1.0_dp)
               return
            end if
            ufunc = ufunc+(real(p,dp)+nu)*log(term)
         else
            ufunc = ufunc+real(n,dp)*dot_product(w,lstar_sq(q,:))
         end if
      end do
      ufunc = ufunc/real(n,dp)
      diff = theta-lw
      value = ufunc-sum(log(values))+dot_product(yvec,dw-dvec) &
           + sum(ymat*transpose(diff)) &
           + 0.5_dp*rho*(vector_norm2(dw-dvec)**2+frobenius_norm(lw-theta)**2)
   end function compute_augmented_lagrangian_ht

   function compute_augmented_lagrangian_kcomp_ht(w,lstar_sq,theta,u,ymat,yvec,dvec, &
                                                   student,nu,rho,beta,k) result(value)
      real(dp), intent(in) :: w(:), lstar_sq(:,:), theta(:,:), u(:,:)
      real(dp), intent(in) :: ymat(:,:), yvec(:), dvec(:), nu, rho, beta
      logical, intent(in) :: student
      integer, intent(in) :: k
      real(dp) :: value
      real(dp), allocatable :: values(:), vectors(:,:), lw(:,:), dw(:), diff(:,:), uut(:,:)
      real(dp) :: ufunc, term
      integer :: q, n, p, status
      n = size(lstar_sq,1)
      p = size(theta,1)
      call symmetric_eigen_jacobi(theta,values,vectors,status)
      if (status /= fg_ok .or. minval(values(k+1:p)) <= 0.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      lw = L(w)
      dw = diagonal(lw)
      ufunc = 0.0_dp
      do q = 1,n
         if (student) then
            term = 1.0_dp+real(n,dp)*dot_product(w,lstar_sq(q,:))/nu
            if (term <= 0.0_dp) then
               value = huge(1.0_dp)
               return
            end if
            ufunc = ufunc+(real(p,dp)+nu)*log(term)
         else
            ufunc = ufunc+real(n,dp)*dot_product(w,lstar_sq(q,:))
         end if
      end do
      ufunc = ufunc/real(n,dp)
      diff = theta-lw
      uut = matmul(u,transpose(u))
      value = ufunc-sum(log(values(k+1:p)))+dot_product(yvec,dw-dvec) &
           + sum(ymat*transpose(diff)) &
           + 0.5_dp*rho*(vector_norm2(dw-dvec)**2+frobenius_norm(lw-theta)**2) &
           + beta*dot_product(w,Lstar(uut))
   end function compute_augmented_lagrangian_kcomp_ht

   subroutine update_theta_connected(lw,jmat,ymat,rho,theta,status)
      real(dp), intent(in) :: lw(:,:), jmat(:,:), ymat(:,:), rho
      real(dp), allocatable, intent(out) :: theta(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: values(:), vectors(:,:), scaled(:,:)
      integer :: i, p
      call symmetric_eigen_jacobi(rho*(lw+jmat)-ymat,values,vectors,status)
      if (status /= fg_ok) then
         allocate(theta(0,0))
         return
      end if
      p = size(values)
      allocate(scaled(p,p))
      scaled = vectors
      do i = 1,p
         scaled(:,i) = scaled(:,i)*(values(i)+sqrt(values(i)**2+4.0_dp*rho))/(2.0_dp*rho)
      end do
      theta = matmul(scaled,transpose(vectors))-jmat
      theta = 0.5_dp*(theta+transpose(theta))
   end subroutine update_theta_connected

   subroutine update_theta_kcomp(lw,ymat,rho,k,theta,status)
      real(dp), intent(in) :: lw(:,:), ymat(:,:), rho
      integer, intent(in) :: k
      real(dp), allocatable, intent(out) :: theta(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: values(:), vectors(:,:), selected(:,:), scaled(:,:)
      integer :: i, p, q
      call symmetric_eigen_jacobi(rho*lw-ymat,values,vectors,status)
      if (status /= fg_ok) then
         allocate(theta(0,0))
         return
      end if
      p = size(values)
      q = p-k
      selected = vectors(:,k+1:p)
      allocate(scaled(p,q))
      scaled = selected
      do i = 1,q
         scaled(:,i) = scaled(:,i)*(values(k+i)+sqrt(values(k+i)**2+4.0_dp*rho))/(2.0_dp*rho)
      end do
      theta = matmul(scaled,transpose(selected))
      theta = 0.5_dp*(theta+transpose(theta))
   end subroutine update_theta_kcomp

   subroutine weighted_scatter_adjoint(w,lstar_sq,p,nu,student,weighted)
      real(dp), intent(in) :: w(:), lstar_sq(:,:), nu
      integer, intent(in) :: p
      logical, intent(in) :: student
      real(dp), intent(out) :: weighted(:)
      integer :: q
      weighted = 0.0_dp
      do q = 1,size(lstar_sq,1)
         if (student) then
            weighted = weighted+lstar_sq(q,:)*compute_student_weight(w,lstar_sq(q,:),p,nu)
         else
            weighted = weighted+lstar_sq(q,:)
         end if
      end do
   end subroutine weighted_scatter_adjoint

   subroutine build_lstar_scatter(x,divisor,lstar_sq)
      real(dp), intent(in) :: x(:,:), divisor
      real(dp), allocatable, intent(out) :: lstar_sq(:,:)
      real(dp), allocatable :: row_outer(:,:), temp(:)
      integer :: i, n, p, m
      n = size(x,1)
      p = size(x,2)
      m = p*(p-1)/2
      allocate(lstar_sq(n,m))
      do i = 1,n
         row_outer = outer_product(x(i,:),x(i,:))
         temp = Lstar(row_outer)
         lstar_sq(i,:) = temp/divisor
      end do
   end subroutine build_lstar_scatter

   subroutine initialize_weights(s,initial_weights,initialization,w,status)
      real(dp), intent(in) :: s(:,:)
      real(dp), intent(in), optional :: initial_weights(:)
      character(len=*), intent(in), optional :: initialization
      real(dp), allocatable, intent(out) :: w(:)
      integer, intent(out) :: status
      real(dp), allocatable :: sinv(:,:), r(:,:), h(:,:), b(:)
      integer :: p, m, info
      character(len=:), allocatable :: method
      p = size(s,1)
      m = p*(p-1)/2
      if (present(initial_weights)) then
         if (size(initial_weights) /= m .or. .not. all(ieee_is_finite(initial_weights))) then
            allocate(w(0))
            status = fg_size_mismatch
            return
         end if
         w = max(initial_weights,0.0_dp)
         status = fg_ok
         return
      end if
      method = "naive"
      if (present(initialization)) method = lowercase(trim(initialization))
      call symmetric_pseudoinverse(s,sinv,info)
      if (info /= fg_ok) then
         allocate(w(0))
         status = info
         return
      end if
      select case (method)
      case ("naive")
         w = max(Linv(sinv),0.0_dp)
      case ("qp")
         r = vecLmat(p)
         h = matmul(transpose(r),r)
         b = matmul(transpose(r),vec(sinv))
         allocate(w(m))
         w = max(Linv(sinv),0.0_dp)
         call nnqp_projected_gradient(h,b,w,info,max_iterations=10000,tolerance=1.0e-11_dp)
         if (info /= fg_ok .and. info /= fg_no_convergence) then
            status = info
            return
         end if
      case default
         allocate(w(0))
         status = fg_invalid_input
         return
      end select
      status = fg_ok
   end subroutine initialize_weights

   subroutine normalize_initial_weights(w,status)
      real(dp), intent(inout) :: w(:)
      integer, intent(out) :: status
      real(dp), allocatable :: a0(:,:), rowsum(:)
      integer :: i, p
      a0 = A(w)
      p = size(a0,1)
      if (p < 2) then
         status = fg_invalid_input
         return
      end if
      allocate(rowsum(p))
      rowsum = sum(a0,dim=2)
      do i = 1,p
         if (rowsum(i) > tiny(1.0_dp)) then
            a0(i,:) = a0(i,:)/rowsum(i)
         else
            a0(i,:) = 1.0_dp/real(p-1,dp)
            a0(i,i) = 0.0_dp
         end if
      end do
      w = max(Ainv(a0),0.0_dp)
      status = fg_ok
   end subroutine normalize_initial_weights

   subroutine build_degree_vector(p,degree,degrees,dvec,status)
      integer, intent(in) :: p
      real(dp), intent(in), optional :: degree, degrees(:)
      real(dp), allocatable, intent(out) :: dvec(:)
      integer, intent(out) :: status
      allocate(dvec(p))
      if (present(degrees)) then
         if (size(degrees) /= p .or. any(degrees <= 0.0_dp) .or. &
             .not. all(ieee_is_finite(degrees))) then
            status = fg_size_mismatch
            return
         end if
         dvec = degrees
      else
         dvec = 1.0_dp
         if (present(degree)) dvec = degree
         if (any(dvec <= 0.0_dp) .or. .not. all(ieee_is_finite(dvec))) then
            status = fg_invalid_input
            return
         end if
      end if
      status = fg_ok
   end subroutine build_degree_vector

   subroutine column_correlation(x,cor,status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: cor(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: xt(:,:)
      xt = transpose(x)
      call rows_correlation(xt,cor,status)
   end subroutine column_correlation

   subroutine standardize_columns(x,z,status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: z(:,:)
      integer, intent(out) :: status
      real(dp) :: mean_value, sd_value
      integer :: j, n
      n = size(x,1)
      allocate(z(size(x,1),size(x,2)))
      do j = 1,size(x,2)
         mean_value = sum(x(:,j))/real(n,dp)
         sd_value = sqrt(sum((x(:,j)-mean_value)**2)/real(n-1,dp))
         if (sd_value <= sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(mean_value))) then
            status = fg_invalid_input
            return
         end if
         z(:,j) = (x(:,j)-mean_value)/sd_value
      end do
      status = fg_ok
   end subroutine standardize_columns

   pure function diagonal(a) result(d)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: d(min(size(a,1),size(a,2)))
      integer :: i
      do i = 1,size(d)
         d(i) = a(i,i)
      end do
   end function diagonal

   pure function lowercase(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code
      do i = 1,len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            lower(i:i) = achar(code+iachar('a')-iachar('A'))
         else
            lower(i:i) = text(i:i)
         end if
      end do
   end function lowercase

   subroutine clear_result(result)
      type(fingraph_result), intent(out) :: result
      result%convergence = .false.
      result%iterations = 0
      result%status = fg_ok
      result%rho = 0.0_dp
      result%beta = 0.0_dp
   end subroutine clear_result

end module fingraph_learning
