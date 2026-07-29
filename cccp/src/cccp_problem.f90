! SPDX-License-Identifier: GPL-3.0-or-later
module cccp_problem
   use cccp_kinds, only : dp
   use cccp_types
   use cccp_solver, only : solve_lp, solve_qp, solve_dnl, solve_dcp
   implicit none
   private

   type, public :: dlp_problem
      real(dp), allocatable :: q(:), a(:,:), b(:)
      type(cone_constraint), allocatable :: cones(:)
   end type dlp_problem
   type, public :: dqp_problem
      real(dp), allocatable :: p(:,:), q(:), a(:,:), b(:)
      type(cone_constraint), allocatable :: cones(:)
   end type dqp_problem

   public :: dlp, dqp, cps
   interface cps
      module procedure cps_lp
      module procedure cps_qp
   end interface cps

contains
   function dlp(q,a,b,cones) result(prob)
      real(dp),intent(in)::q(:)
      real(dp),intent(in),optional::a(:,:),b(:)
      type(cone_constraint),intent(in),optional::cones(:)
      type(dlp_problem)::prob
      integer::n
      n=size(q);allocate(prob%q(n));prob%q=q
      if(present(a))then;allocate(prob%a(size(a,1),n),prob%b(size(a,1)));prob%a=a;prob%b=b
      else;allocate(prob%a(0,n),prob%b(0));end if
      if(present(cones))then;allocate(prob%cones(size(cones)));prob%cones=cones
      else;allocate(prob%cones(0));end if
   end function dlp
   function dqp(p,q,a,b,cones) result(prob)
      real(dp),intent(in)::p(:,:),q(:)
      real(dp),intent(in),optional::a(:,:),b(:)
      type(cone_constraint),intent(in),optional::cones(:)
      type(dqp_problem)::prob
      integer::n
      n=size(q);allocate(prob%p(n,n),prob%q(n));prob%p=p;prob%q=q
      if(present(a))then;allocate(prob%a(size(a,1),n),prob%b(size(a,1)));prob%a=a;prob%b=b
      else;allocate(prob%a(0,n),prob%b(0));end if
      if(present(cones))then;allocate(prob%cones(size(cones)));prob%cones=cones
      else;allocate(prob%cones(0));end if
   end function dqp
   function cps_lp(prob,control) result(sol)
      type(dlp_problem),intent(in)::prob
      type(cccp_control),intent(in),optional::control
      type(cccp_solution)::sol
      call solve_lp(prob%q,prob%a,prob%b,prob%cones,control,sol)
   end function cps_lp
   function cps_qp(prob,control) result(sol)
      type(dqp_problem),intent(in)::prob
      type(cccp_control),intent(in),optional::control
      type(cccp_solution)::sol
      call solve_qp(prob%p,prob%q,prob%a,prob%b,prob%cones,control,sol)
   end function cps_qp
end module cccp_problem
