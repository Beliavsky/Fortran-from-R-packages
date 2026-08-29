! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module nnet_fit_mod
use r_compat, only: dp, optim_result_t, optim_bfgs, runif_vec
use nnet_types, only: nnet_model_t
use nnet_core, only: build_network, nnet_objective, nnet_gradient, nnet_predict_raw, nnet_hessian_exact
implicit none
private
public :: nnet_fit, nnet_refit, nnet_predict, nnet_predict_class

type(nnet_model_t), save :: ctx_model
real(dp), allocatable, save :: ctx_x(:,:), ctx_y(:,:), ctx_case_weights(:), ctx_base_wts(:)
integer, allocatable, save :: ctx_var_idx(:)

contains

subroutine nnet_fit(model, x, y, hidden_size, case_weights, initial_wts, mask, linout, entropy, softmax, censored, skip, &
   rang, decay, maxit, max_nwts, abstol, reltol, hessian)
type(nnet_model_t), intent(out) :: model
real(dp), intent(in) :: x(:,:), y(:,:)
integer, intent(in) :: hidden_size
real(dp), intent(in), optional :: case_weights(:), initial_wts(:), decay(:)
logical, intent(in), optional :: mask(:), linout, entropy, softmax, censored, skip, hessian
real(dp), intent(in), optional :: rang, abstol, reltol
integer, intent(in), optional :: maxit, max_nwts
real(dp), allocatable :: cw(:), u(:)
real(dp) :: rr
integer :: nw, mnw
logical :: llin, lent, lsoft, lcens, lskip
llin=.false.
if(present(linout)) llin=linout
lent=.false.
if(present(entropy)) lent=entropy
lsoft=.false.
if(present(softmax)) lsoft=softmax
lcens=.false.
if(present(censored)) lcens=censored
lskip=.false.
if(present(skip)) lskip=skip
if (size(x,1) /= size(y,1)) error stop "nnet_fit: x and y row counts differ"
if (size(x,1) < 1 .or. size(x,2) < 1 .or. size(y,2) < 1) error stop "nnet_fit: empty data"
if (llin .and. lent) error stop "nnet_fit: entropy fit only for logistic output units"
if (lsoft .and. size(y,2) < 2) error stop "nnet_fit: softmax requires at least two outputs"
call build_network(model, size(x,2), hidden_size, size(y,2), llin, lent, lsoft, lcens, lskip)
nw = size(model%wts)
if (nw == 0) error stop "nnet_fit: no weights to fit"
mnw=1000
if(present(max_nwts)) mnw=max_nwts
if(nw>mnw) error stop "nnet_fit: too many weights"
rr = 0.7_dp
if (present(rang)) rr = rang
if (present(initial_wts)) then
   if (size(initial_wts) /= nw) error stop "nnet_fit: incorrect initial_wts length"
   model%wts = initial_wts
else if (rr > 0.0_dp) then
   u = runif_vec(nw)
   model%wts = (2.0_dp*u - 1.0_dp) * rr
else
   model%wts = 0.0_dp
end if
if (present(mask)) then
   if (size(mask) /= nw) error stop "nnet_fit: incorrect mask length"
   model%mask = mask
else
   model%mask = .true.
end if
if (present(decay)) then
   if (size(decay) == 1) then
      model%decay = decay(1)
   else if (size(decay) == nw) then
      model%decay = decay
   else
      error stop "nnet_fit: decay must have length 1 or number of weights"
   end if
else
   model%decay = 0.0_dp
end if
allocate(cw(size(x,1)))
if (present(case_weights)) then
   if (size(case_weights) /= size(x,1) .or. any(case_weights < 0.0_dp)) error stop "nnet_fit: invalid case weights"
   cw = case_weights
else
   cw = 1.0_dp
end if
call fit_existing(model, x, y, cw, maxit, abstol, reltol)
model%fitted = nnet_predict_raw(model, x, model%wts)
allocate(model%residuals(size(y,1),size(y,2)))
model%residuals = y - model%fitted
if (present(hessian)) then
   if (hessian) model%hessian = nnet_hessian_exact(model, x, y, cw, model%wts)
end if
end subroutine nnet_fit

subroutine nnet_refit(model, x, y, case_weights, maxit, abstol, reltol, hessian)
type(nnet_model_t), intent(inout) :: model
real(dp), intent(in) :: x(:,:), y(:,:)
real(dp), intent(in), optional :: case_weights(:)
integer, intent(in), optional :: maxit
real(dp), intent(in), optional :: abstol, reltol
logical, intent(in), optional :: hessian
real(dp), allocatable :: cw(:)
allocate(cw(size(x,1)))
cw=1.0_dp
if(present(case_weights)) cw=case_weights
call fit_existing(model,x,y,cw,maxit,abstol,reltol)
model%fitted = nnet_predict_raw(model,x,model%wts)
if (allocated(model%residuals)) deallocate(model%residuals)
allocate(model%residuals(size(y,1),size(y,2)))
model%residuals=y-model%fitted
if (present(hessian)) then
   if(hessian) model%hessian=nnet_hessian_exact(model,x,y,cw,model%wts)
end if
end subroutine nnet_refit

subroutine fit_existing(model, x, y, cw, maxit, abstol, reltol)
type(nnet_model_t), intent(inout) :: model
real(dp), intent(in) :: x(:,:), y(:,:), cw(:)
integer, intent(in), optional :: maxit
real(dp), intent(in), optional :: abstol, reltol
type(optim_result_t) :: opt
real(dp), allocatable :: p0(:)
integer :: i, k, mi
real(dp) :: rt, at
mi=100
if(present(maxit)) mi=maxit
rt=1.0e-8_dp
if(present(reltol)) rt=reltol
at=1.0e-4_dp
if(present(abstol)) at=abstol
ctx_model=model
ctx_x=x
ctx_y=y
ctx_case_weights=cw
ctx_base_wts=model%wts
ctx_var_idx=pack([(i,i=1,size(model%wts))],model%mask)
allocate(p0(size(ctx_var_idx)))
do k=1,size(ctx_var_idx)
p0(k)=model%wts(ctx_var_idx(k))
end do
if(size(p0)==0) then
   model%value=nnet_objective(model,x,y,cw,model%wts)
   model%convergence=0
   model%counts=0
   return
end if
if (nnet_objective(model,x,y,cw,model%wts) <= at) then
   model%value=nnet_objective(model,x,y,cw,model%wts)
   model%convergence=0
   model%counts=[1,0]
   return
end if
opt=optim_bfgs(fit_objective_var,p0,maxit=mi,reltol=rt,gr=fit_gradient_var)
do k=1,size(ctx_var_idx)
model%wts(ctx_var_idx(k))=opt%par(k)
end do
model%value=opt%value
model%convergence=opt%convergence
model%counts=opt%counts
end subroutine fit_existing

pure function fit_objective_var(par) result(value)
real(dp), intent(in) :: par(:)
real(dp) :: value
real(dp), allocatable :: full(:)
integer :: k
full=ctx_base_wts
do k=1,size(ctx_var_idx)
full(ctx_var_idx(k))=par(k)
end do
value=nnet_objective(ctx_model,ctx_x,ctx_y,ctx_case_weights,full)
end function fit_objective_var

pure function fit_gradient_var(par) result(g)
real(dp), intent(in) :: par(:)
real(dp), allocatable :: g(:), full(:), gf(:)
integer :: k
full=ctx_base_wts
do k=1,size(ctx_var_idx)
full(ctx_var_idx(k))=par(k)
end do
gf=nnet_gradient(ctx_model,ctx_x,ctx_y,ctx_case_weights,full)
allocate(g(size(par)))
do k=1,size(ctx_var_idx)
g(k)=gf(ctx_var_idx(k))
end do
end function fit_gradient_var

pure function nnet_predict(model,x) result(pred)
type(nnet_model_t), intent(in) :: model
real(dp), intent(in) :: x(:,:)
real(dp), allocatable :: pred(:,:)
if(size(x,2)/=model%n_inputs) error stop "nnet_predict: wrong number of input columns"
pred=nnet_predict_raw(model,x,model%wts)
end function nnet_predict

function nnet_predict_class(model,x) result(cls)
type(nnet_model_t), intent(in) :: model
real(dp), intent(in) :: x(:,:)
integer, allocatable :: cls(:)
real(dp), allocatable :: p(:,:)
integer :: i
p=nnet_predict(model,x)
allocate(cls(size(x,1)))
if(size(p,2)==1) then
   cls=1+merge(1,0,p(:,1)>0.5_dp)
else
   do i=1,size(x,1)
   cls(i)=maxloc(p(i,:),dim=1)
   end do
end if
end function nnet_predict_class

end module nnet_fit_mod
