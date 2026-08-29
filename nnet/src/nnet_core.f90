! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of the numerical core in nnet/src/nnet.c.
module nnet_core
use r_compat, only: dp
use nnet_types, only: nnet_model_t
implicit none
private
public :: build_network, nnet_predict_raw, nnet_objective, nnet_gradient, nnet_objective_gradient
public :: nnet_hessian_exact, nnet_weight_count

real(dp), parameter :: eps_prob = 1.0e-80_dp

contains

pure integer function nnet_weight_count(n_inputs, n_hidden, n_outputs, skip) result(nw)
integer, intent(in) :: n_inputs, n_hidden, n_outputs
logical, intent(in) :: skip
if (n_hidden > 0) then
   nw = n_hidden * (n_inputs + 1) + n_outputs * (n_hidden + 1)
   if (skip) nw = nw + n_outputs * n_inputs
else if (skip) then
   nw = n_outputs * (n_inputs + 1)
else
   nw = 0
end if
end function nnet_weight_count

subroutine build_network(model, n_inputs, n_hidden, n_outputs, linout, entropy, softmax, censored, skip)
type(nnet_model_t), intent(out) :: model
integer, intent(in) :: n_inputs, n_hidden, n_outputs
logical, intent(in), optional :: linout, entropy, softmax, censored, skip
logical :: llin, lent, lsoft, lcens, lskip
integer :: j, pos, k, nw
llin = .false.
if (present(linout)) llin = linout
lent = .false.
if (present(entropy)) lent = entropy
lsoft = .false.
if (present(softmax)) lsoft = softmax
lcens = .false.
if (present(censored)) lcens = censored
lskip = .false.
if (present(skip)) lskip = skip
if (lsoft) then
   llin = .true.
   lent = .false.
end if
if (lcens) then
   llin = .true.
   lent = .false.
   lsoft = .true.
end if
model%n_inputs = n_inputs
model%n_hidden = n_hidden
model%n_outputs = n_outputs
model%n_units = 1 + n_inputs + n_hidden + n_outputs
model%first_hidden = 1 + n_inputs
model%first_output = 1 + n_inputs + n_hidden
model%ns_units = model%n_units
if (llin) model%ns_units = model%n_units - n_outputs
model%entropy = lent
model%softmax = lsoft
model%censored = lcens
model%skip = lskip
nw = nnet_weight_count(n_inputs, n_hidden, n_outputs, lskip)
allocate(model%nconn(0:model%n_units), model%conn(nw), model%wts(nw), model%decay(nw), model%mask(nw))
model%nconn = 0
pos = 0
model%nconn(0) = 0
! destinations 0 through first_hidden-1 have no incoming fitted weights
do j = 0, model%first_hidden - 1
   model%nconn(j + 1) = pos
end do
! hidden units: bias followed by all inputs
do j = model%first_hidden, model%first_output - 1
   pos = pos + 1
   model%conn(pos) = 0
   do k = 1, n_inputs
      pos = pos + 1
      model%conn(pos) = k
   end do
   model%nconn(j + 1) = pos
end do
! output units
do j = model%first_output, model%n_units - 1
   if (n_hidden > 0) then
      pos = pos + 1
      model%conn(pos) = 0
      do k = model%first_hidden, model%first_output - 1
         pos = pos + 1
         model%conn(pos) = k
      end do
      if (lskip) then
         do k = 1, n_inputs
            pos = pos + 1
            model%conn(pos) = k
         end do
      end if
   else if (lskip) then
      pos = pos + 1
      model%conn(pos) = 0
      do k = 1, n_inputs
         pos = pos + 1
         model%conn(pos) = k
      end do
   end if
   model%nconn(j + 1) = pos
end do
model%wts = 0.0_dp
model%decay = 0.0_dp
model%mask = .true.
end subroutine build_network

pure elemental real(dp) function sigmoid(x) result(y)
real(dp), intent(in) :: x
if (x < -15.0_dp) then
   y = 0.0_dp
else if (x > 15.0_dp) then
   y = 1.0_dp
else
   y = 1.0_dp / (1.0_dp + exp(-x))
end if
end function sigmoid

pure elemental real(dp) function sigmoid_prime(value) result(y)
real(dp), intent(in) :: value
y = value * (1.0_dp - value)
end function sigmoid_prime

pure elemental real(dp) function sigmoid_prime_prime(value) result(y)
real(dp), intent(in) :: value
y = value * (1.0_dp - value) * (1.0_dp - 2.0_dp * value)
end function sigmoid_prime_prime

pure subroutine forward_values(model, input, wts, outputs, probs)
type(nnet_model_t), intent(in) :: model
real(dp), intent(in) :: input(:), wts(:)
real(dp), intent(out) :: outputs(0:), probs(0:)
integer :: j, iw
real(dp) :: s, mx, den
outputs = 0.0_dp
probs = 0.0_dp
outputs(0) = 1.0_dp
outputs(1:model%n_inputs) = input(1:model%n_inputs)
do j = model%first_hidden, model%n_units - 1
   s = 0.0_dp
   do iw = model%nconn(j) + 1, model%nconn(j + 1)
      s = s + outputs(model%conn(iw)) * wts(iw)
   end do
   if (j < model%ns_units) s = sigmoid(s)
   outputs(j) = s
end do
if (model%softmax) then
   mx = maxval(outputs(model%first_output:model%n_units - 1))
   den = 0.0_dp
   do j = model%first_output, model%n_units - 1
      probs(j) = exp(outputs(j) - mx)
      den = den + probs(j)
   end do
   probs(model%first_output:model%n_units - 1) = probs(model%first_output:model%n_units - 1) / den
end if
end subroutine forward_values

pure real(dp) function sample_error(model, goal, wx, outputs, probs) result(err)
type(nnet_model_t), intent(in) :: model
real(dp), intent(in) :: goal(:), wx
real(dp), intent(in) :: outputs(0:), probs(0:)
integer :: i, k
real(dp) :: t, d, s
err = 0.0_dp
if (model%softmax) then
   if (model%censored) then
      s = 0.0_dp
      do i = model%first_output, model%n_units - 1
         k = i - model%first_output + 1
         if (goal(k) == 1.0_dp) s = s + probs(i)
      end do
      if (s > 0.0_dp) then
         err = -wx * log(s)
      else
         err = wx * 1000.0_dp
      end if
   else
      do i = model%first_output, model%n_units - 1
         k = i - model%first_output + 1
         t = goal(k)
         if (t > 0.0_dp) then
            if (probs(i) > 0.0_dp) then
               err = err - wx * t * log(probs(i))
            else
               err = err + wx * 1000.0_dp
            end if
         end if
      end do
   end if
else
   do i = model%first_output, model%n_units - 1
      k = i - model%first_output + 1
      t = goal(k)
      if (model%entropy) then
         if (t > 0.0_dp) err = err - wx * t * log((outputs(i) + eps_prob) / t)
         if (t < 1.0_dp) err = err - wx * (1.0_dp - t) * log((1.0_dp - outputs(i) + eps_prob) / (1.0_dp - t))
      else
         d = outputs(i) - t
         err = err + wx * d * d
      end if
   end do
end if
end function sample_error

pure subroutine backward_values(model, goal, wx, outputs, probs, wts, slopes, errors)
type(nnet_model_t), intent(in) :: model
real(dp), intent(in) :: goal(:), wx, outputs(0:), probs(0:), wts(:)
real(dp), intent(inout) :: slopes(:)
real(dp), intent(out) :: errors(0:)
real(dp) :: error_sums(0:model%n_units-1)
integer :: i, j, iw, src, k
real(dp) :: s, denom
error_sums = 0.0_dp
errors = 0.0_dp
if (model%softmax) then
   if (model%censored) then
      denom = 0.0_dp
      do i = model%first_output, model%n_units - 1
         k = i - model%first_output + 1
         if (goal(k) == 1.0_dp) denom = denom + probs(i)
      end do
      do i = model%first_output, model%n_units - 1
         k = i - model%first_output + 1
         error_sums(i) = probs(i)
         if (goal(k) == 1.0_dp) error_sums(i) = error_sums(i) - probs(i) / denom
      end do
   else
      s = sum(goal)
      do i = model%first_output, model%n_units - 1
         k = i - model%first_output + 1
         error_sums(i) = s * probs(i) - goal(k)
      end do
   end if
else if (model%entropy) then
   do i = model%first_output, model%n_units - 1
      k = i - model%first_output + 1
      error_sums(i) = outputs(i) - goal(k)
   end do
else
   do i = model%first_output, model%n_units - 1
      k = i - model%first_output + 1
      error_sums(i) = 2.0_dp * (outputs(i) - goal(k))
      if (i < model%ns_units) error_sums(i) = error_sums(i) * sigmoid_prime(outputs(i))
   end do
end if
error_sums(model%first_hidden:model%first_output-1) = 0.0_dp
do j = model%n_units - 1, model%first_hidden, -1
   errors(j) = error_sums(j)
   if (j < model%first_output) errors(j) = errors(j) * sigmoid_prime(outputs(j))
   do iw = model%nconn(j) + 1, model%nconn(j + 1)
      src = model%conn(iw)
      error_sums(src) = error_sums(src) + errors(j) * wts(iw)
      slopes(iw) = slopes(iw) + wx * errors(j) * outputs(src)
   end do
end do
end subroutine backward_values

pure function nnet_objective(model, x, y, case_weights, wts) result(value)
type(nnet_model_t), intent(in) :: model
real(dp), intent(in) :: x(:,:), y(:,:), case_weights(:), wts(:)
real(dp) :: value
real(dp) :: outputs(0:model%n_units-1), probs(0:model%n_units-1)
integer :: i
value = 0.0_dp
do i = 1, size(x,1)
   call forward_values(model, x(i,:), wts, outputs, probs)
   value = value + sample_error(model, y(i,:), case_weights(i), outputs, probs)
end do
value = value + sum(model%decay * wts * wts)
end function nnet_objective

pure function nnet_gradient(model, x, y, case_weights, wts) result(grad)
type(nnet_model_t), intent(in) :: model
real(dp), intent(in) :: x(:,:), y(:,:), case_weights(:), wts(:)
real(dp), allocatable :: grad(:)
real(dp) :: outputs(0:model%n_units-1), probs(0:model%n_units-1), errors(0:model%n_units-1)
integer :: i
allocate(grad(size(wts)))
grad = 2.0_dp * model%decay * wts
do i = 1, size(x,1)
   call forward_values(model, x(i,:), wts, outputs, probs)
   call backward_values(model, y(i,:), case_weights(i), outputs, probs, wts, grad, errors)
end do
end function nnet_gradient

pure subroutine nnet_objective_gradient(model, x, y, case_weights, wts, value, grad)
type(nnet_model_t), intent(in) :: model
real(dp), intent(in) :: x(:,:), y(:,:), case_weights(:), wts(:)
real(dp), intent(out) :: value
real(dp), allocatable, intent(out) :: grad(:)
value = nnet_objective(model, x, y, case_weights, wts)
grad = nnet_gradient(model, x, y, case_weights, wts)
end subroutine nnet_objective_gradient

pure function nnet_predict_raw(model, x, wts) result(pred)
type(nnet_model_t), intent(in) :: model
real(dp), intent(in) :: x(:,:), wts(:)
real(dp), allocatable :: pred(:,:)
real(dp) :: outputs(0:model%n_units-1), probs(0:model%n_units-1)
integer :: i, j, k
allocate(pred(size(x,1), model%n_outputs))
do i = 1, size(x,1)
   call forward_values(model, x(i,:), wts, outputs, probs)
   do j = model%first_output, model%n_units - 1
      k = j - model%first_output + 1
      if (model%softmax) then
         pred(i,k) = probs(j)
      else
         pred(i,k) = outputs(j)
      end if
   end do
end do
end function nnet_predict_raw

pure function nnet_hessian_exact(model, x, y, case_weights, wts) result(hess)
type(nnet_model_t), intent(in) :: model
real(dp), intent(in) :: x(:,:), y(:,:), case_weights(:), wts(:)
real(dp), allocatable :: hess(:,:)
real(dp) :: outputs(0:model%n_units-1), probs(0:model%n_units-1), errors(0:model%n_units-1)
real(dp) :: dummy_slopes(size(wts)), hvec(0:model%n_units-1), h1(0:model%n_units-1)
real(dp) :: wmat(0:model%n_units-1,0:model%n_units-1)
real(dp) :: out, s, sum1, sum2, t, tmp, tot, pden
integer :: obs, i, j, iw, to1, to2, from1, from2, j1, j2, k
logical :: first1, first2
allocate(hess(size(wts), size(wts)))
hess = 0.0_dp
wmat = 0.0_dp
do j = model%first_output, model%n_units - 1
   do iw = model%nconn(j) + 1, model%nconn(j+1)
      wmat(model%conn(iw),j) = wts(iw)
   end do
end do
do obs = 1, size(x,1)
   call forward_values(model, x(obs,:), wts, outputs, probs)
   dummy_slopes = 0.0_dp
   call backward_values(model, y(obs,:), 1.0_dp, outputs, probs, wts, dummy_slopes, errors)
   hvec = 0.0_dp
   h1 = 0.0_dp
   if (model%softmax) then
      tot = 0.0_dp
      pden = 0.0_dp
      do i = 0, model%n_units - 1
         sum1 = 0.0_dp
         sum2 = 0.0_dp
         tot = 0.0_dp
         pden = 0.0_dp
         do j = model%first_output, model%n_units - 1
            sum1 = sum1 + wmat(i,j) * probs(j)
            k = j - model%first_output + 1
            t = y(obs,k)
            pden = pden + t * probs(j)
            sum2 = sum2 + wmat(i,j) * probs(j) * t
            tot = tot + t
         end do
         hvec(i) = sum1
         if (pden /= 0.0_dp) h1(i) = sum2 / pden
      end do
      if (model%censored) tot = 1.0_dp
      do to1 = 0, model%n_units - 1
         do j1 = model%nconn(to1)+1, model%nconn(to1+1)
            from1 = model%conn(j1)
            first1 = to1 < model%first_output
            do to2 = 0, model%n_units - 1
               do j2 = model%nconn(to2)+1, model%nconn(to2+1)
                  if (j2 > j1) cycle
                  from2 = model%conn(j2)
                  first2 = to2 < model%first_output
                  if ((.not. first1) .and. (.not. first2)) then
                     if (model%censored) then
                        tmp = -probs(to1)*probs(to2) * (1.0_dp - &
                           y(obs,to1-model%first_output+1)*y(obs,to2-model%first_output+1)/(pden*pden))
                        if (to1 == to2) tmp = tmp + probs(to1) * &
                           (1.0_dp - y(obs,to1-model%first_output+1)/pden)
                        hess(j1,j2) = hess(j1,j2) + case_weights(obs)*tmp*outputs(from1)*outputs(from2)
                     else
                        tmp = -probs(to1)*probs(to2)
                        if (to1 == to2) tmp = tmp + probs(to1)
                        hess(j1,j2) = hess(j1,j2) + case_weights(obs)*tot*tmp*outputs(from1)*outputs(from2)
                     end if
                  else if (first1 .and. first2) then
                     sum1 = 0.0_dp
                     sum2 = 0.0_dp
                     do i = model%first_output, model%n_units - 1
                        sum1 = sum1 + errors(i)*wmat(to1,i)
                        tmp = wmat(to1,i)*wmat(to2,i)*probs(i)
                        if (model%censored) tmp = tmp * (1.0_dp - y(obs,i-model%first_output+1)/pden)
                        sum2 = sum2 + tmp
                     end do
                     if (model%censored) then
                        sum2 = sum2 - hvec(to1)*hvec(to2) + h1(to1)*h1(to2)
                        s = sigmoid_prime(outputs(to1))*sigmoid_prime(outputs(to2))*sum2
                     else
                        sum2 = sum2 - hvec(to1)*hvec(to2)
                        s = sigmoid_prime(outputs(to1))*sigmoid_prime(outputs(to2))*tot*sum2
                     end if
                     if (to1 == to2) s = s + sigmoid_prime_prime(outputs(to1))*sum1
                     hess(j1,j2) = hess(j1,j2) + case_weights(obs)*s*outputs(from1)*outputs(from2)
                  else
                     if (to1 < to2) then
                        tmp = wmat(to1,to2) - hvec(to1)
                        if (model%censored) tmp = tmp + y(obs,to2-model%first_output+1)/pden * &
                           (h1(to1)-wmat(to1,to2))
                        hess(j1,j2) = hess(j1,j2) + case_weights(obs)*outputs(from1)*sigmoid_prime(outputs(to1))* &
                           (outputs(from2)*probs(to2)*tmp*tot + merge(errors(to2),0.0_dp,to1==from2))
                     else
                        tmp = wmat(to2,to1) - hvec(to2)
                        if (model%censored) tmp = tmp + y(obs,to1-model%first_output+1)/pden * &
                           (h1(to2)-wmat(to2,to1))
                        hess(j1,j2) = hess(j1,j2) + case_weights(obs)*outputs(from2)*sigmoid_prime(outputs(to2))* &
                           (outputs(from1)*probs(to1)*tmp*tot + merge(errors(to1),0.0_dp,to2==from1))
                     end if
                  end if
               end do
            end do
         end do
      end do
   else
      do i = model%first_output, model%n_units - 1
         out = outputs(i)
         s = sigmoid_prime(out)
         t = y(obs,i-model%first_output+1)
         if (model%ns_units < model%n_units) then
            hvec(i) = 2.0_dp
         else if (model%entropy) then
            hvec(i) = out*(1.0_dp-out)
         else
            hvec(i) = sigmoid_prime_prime(out)*2.0_dp*(out-t) + 2.0_dp*s*s
         end if
      end do
      do to1 = 0, model%n_units - 1
         do j1 = model%nconn(to1)+1, model%nconn(to1+1)
            from1 = model%conn(j1)
            first1 = to1 < model%first_output
            do to2 = 0, model%n_units - 1
               do j2 = model%nconn(to2)+1, model%nconn(to2+1)
                  if (j2 > j1) cycle
                  from2 = model%conn(j2)
                  first2 = to2 < model%first_output
                  if ((.not.first1) .and. (.not.first2)) then
                     if (to1 == to2) hess(j1,j2) = hess(j1,j2) + case_weights(obs)*hvec(to1)*outputs(from1)*outputs(from2)
                  else if (first1 .and. first2) then
                     sum1 = 0.0_dp
                     sum2 = 0.0_dp
                     do i = model%first_output, model%n_units - 1
                        sum1 = sum1 + errors(i)*wmat(to1,i)
                        sum2 = sum2 + wmat(to1,i)*wmat(to2,i)*hvec(i)
                     end do
                     s = sigmoid_prime(outputs(to1))*sigmoid_prime(outputs(to2))*sum2
                     if (to1 == to2) s = s + sigmoid_prime_prime(outputs(to1))*sum1
                     hess(j1,j2) = hess(j1,j2) + case_weights(obs)*s*outputs(from1)*outputs(from2)
                  else if (to1 < to2) then
                     hess(j1,j2) = hess(j1,j2) + case_weights(obs)*outputs(from1)*sigmoid_prime(outputs(to1))* &
                        (outputs(from2)*wmat(to1,to2)*hvec(to2) + merge(errors(to2),0.0_dp,to1==from2))
                  else
                     hess(j1,j2) = hess(j1,j2) + case_weights(obs)*outputs(from2)*sigmoid_prime(outputs(to2))* &
                        (outputs(from1)*wmat(to2,to1)*hvec(to1) + merge(errors(to1),0.0_dp,to2==from1))
                  end if
               end do
            end do
         end do
      end do
   end if
end do
! lower triangle was accumulated; make symmetric and add penalty curvature
do j1 = 1, size(wts)
   do j2 = 1, j1-1
      hess(j2,j1) = hess(j1,j2)
   end do
   hess(j1,j1) = hess(j1,j1) + 2.0_dp*model%decay(j1)
end do
end function nnet_hessian_exact

end module nnet_core
