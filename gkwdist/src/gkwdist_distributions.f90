! SPDX-License-Identifier: MIT
module gkwdist_distributions
   use gkwdist_kinds, only : dp
   use gkwdist_core, only : fam_gkw,fam_bkw,fam_kkw,fam_ekw,fam_mc,fam_kw,fam_beta, &
      dgkw_scalar,pgkw_scalar,qgkw_scalar,rgkw_scalar,family_nll,family_derivatives
   implicit none
   private
   public :: dgkw,pgkw,qgkw,rgkw,llgkw,grgkw,hsgkw
   public :: dbkw,pbkw,qbkw,rbkw,llbkw,grbkw,hsbkw
   public :: dkkw,pkkw,qkkw,rkkw,llkkw,grkkw,hskkw
   public :: dekw,pekw,qekw,rekw,llekw,grekw,hsekw
   public :: dmc,pmc,qmc,rmc,llmc,grmc,hsmc
   public :: dkw,pkw,qkw,rkw,llkw,grkw,hskw
   public :: dbeta_,pbeta_,qbeta_,rbeta_,llbeta,grbeta,hsbeta

contains

   pure elemental function dgkw(x,alpha,beta,gamma,delta,lambda,log_prob) result(v)
      real(dp),intent(in)::x,alpha,beta,gamma,delta,lambda
      logical,intent(in),optional::log_prob
      real(dp)::v
      v=dgkw_scalar(x,alpha,beta,gamma,delta,lambda,log_prob)
   end function
   pure elemental function pgkw(q,alpha,beta,gamma,delta,lambda,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,alpha,beta,gamma,delta,lambda
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::v
      v=pgkw_scalar(q,alpha,beta,gamma,delta,lambda,lower_tail,log_p)
   end function
   pure elemental function qgkw(p,alpha,beta,gamma,delta,lambda,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,alpha,beta,gamma,delta,lambda
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::v
      v=qgkw_scalar(p,alpha,beta,gamma,delta,lambda,lower_tail,log_p)
   end function
   function rgkw(n,alpha,beta,gamma,delta,lambda) result(x)
      integer,intent(in)::n
      real(dp),intent(in)::alpha,beta,gamma,delta,lambda
      real(dp),allocatable::x(:); integer::i
      allocate(x(max(0,n)))
      do i=1,size(x); x(i)=rgkw_scalar(alpha,beta,gamma,delta,lambda); end do
   end function

   pure elemental function dbkw(x,alpha,beta,gamma,delta,log_prob) result(v)
      real(dp),intent(in)::x,alpha,beta,gamma,delta; logical,intent(in),optional::log_prob; real(dp)::v
      v=dgkw_scalar(x,alpha,beta,gamma,delta,1.0_dp,log_prob)
   end function
   pure elemental function pbkw(q,alpha,beta,gamma,delta,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,alpha,beta,gamma,delta; logical,intent(in),optional::lower_tail,log_p; real(dp)::v
      v=pgkw_scalar(q,alpha,beta,gamma,delta,1.0_dp,lower_tail,log_p)
   end function
   pure elemental function qbkw(p,alpha,beta,gamma,delta,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,alpha,beta,gamma,delta; logical,intent(in),optional::lower_tail,log_p; real(dp)::v
      v=qgkw_scalar(p,alpha,beta,gamma,delta,1.0_dp,lower_tail,log_p)
   end function
   function rbkw(n,alpha,beta,gamma,delta) result(x)
      integer,intent(in)::n; real(dp),intent(in)::alpha,beta,gamma,delta; real(dp),allocatable::x(:); integer::i
      allocate(x(max(0,n))); do i=1,size(x); x(i)=rgkw_scalar(alpha,beta,gamma,delta,1.0_dp); end do
   end function

   pure elemental function dkkw(x,alpha,beta,delta,lambda,log_prob) result(v)
      real(dp),intent(in)::x,alpha,beta,delta,lambda; logical,intent(in),optional::log_prob; real(dp)::v
      v=dgkw_scalar(x,alpha,beta,1.0_dp,delta,lambda,log_prob)
   end function
   pure elemental function pkkw(q,alpha,beta,delta,lambda,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,alpha,beta,delta,lambda; logical,intent(in),optional::lower_tail,log_p; real(dp)::v
      v=pgkw_scalar(q,alpha,beta,1.0_dp,delta,lambda,lower_tail,log_p)
   end function
   pure elemental function qkkw(p,alpha,beta,delta,lambda,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,alpha,beta,delta,lambda; logical,intent(in),optional::lower_tail,log_p; real(dp)::v
      v=qgkw_scalar(p,alpha,beta,1.0_dp,delta,lambda,lower_tail,log_p)
   end function
   function rkkw(n,alpha,beta,delta,lambda) result(x)
      integer,intent(in)::n; real(dp),intent(in)::alpha,beta,delta,lambda; real(dp),allocatable::x(:); integer::i
      allocate(x(max(0,n))); do i=1,size(x); x(i)=rgkw_scalar(alpha,beta,1.0_dp,delta,lambda); end do
   end function

   pure elemental function dekw(x,alpha,beta,lambda,log_prob) result(v)
      real(dp),intent(in)::x,alpha,beta,lambda; logical,intent(in),optional::log_prob; real(dp)::v
      v=dgkw_scalar(x,alpha,beta,1.0_dp,0.0_dp,lambda,log_prob)
   end function
   pure elemental function pekw(q,alpha,beta,lambda,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,alpha,beta,lambda; logical,intent(in),optional::lower_tail,log_p; real(dp)::v
      v=pgkw_scalar(q,alpha,beta,1.0_dp,0.0_dp,lambda,lower_tail,log_p)
   end function
   pure elemental function qekw(p,alpha,beta,lambda,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,alpha,beta,lambda; logical,intent(in),optional::lower_tail,log_p; real(dp)::v
      v=qgkw_scalar(p,alpha,beta,1.0_dp,0.0_dp,lambda,lower_tail,log_p)
   end function
   function rekw(n,alpha,beta,lambda) result(x)
      integer,intent(in)::n; real(dp),intent(in)::alpha,beta,lambda; real(dp),allocatable::x(:); integer::i
      allocate(x(max(0,n))); do i=1,size(x); x(i)=rgkw_scalar(alpha,beta,1.0_dp,0.0_dp,lambda); end do
   end function

   pure elemental function dmc(x,gamma,delta,lambda,log_prob) result(v)
      real(dp),intent(in)::x,gamma,delta,lambda; logical,intent(in),optional::log_prob; real(dp)::v
      v=dgkw_scalar(x,1.0_dp,1.0_dp,gamma,delta,lambda,log_prob)
   end function
   pure elemental function pmc(q,gamma,delta,lambda,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,gamma,delta,lambda; logical,intent(in),optional::lower_tail,log_p; real(dp)::v
      v=pgkw_scalar(q,1.0_dp,1.0_dp,gamma,delta,lambda,lower_tail,log_p)
   end function
   pure elemental function qmc(p,gamma,delta,lambda,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,gamma,delta,lambda; logical,intent(in),optional::lower_tail,log_p; real(dp)::v
      v=qgkw_scalar(p,1.0_dp,1.0_dp,gamma,delta,lambda,lower_tail,log_p)
   end function
   function rmc(n,gamma,delta,lambda) result(x)
      integer,intent(in)::n; real(dp),intent(in)::gamma,delta,lambda; real(dp),allocatable::x(:); integer::i
      allocate(x(max(0,n))); do i=1,size(x); x(i)=rgkw_scalar(1.0_dp,1.0_dp,gamma,delta,lambda); end do
   end function

   pure elemental function dkw(x,alpha,beta,log_prob) result(v)
      real(dp),intent(in)::x,alpha,beta; logical,intent(in),optional::log_prob; real(dp)::v
      v=dgkw_scalar(x,alpha,beta,1.0_dp,0.0_dp,1.0_dp,log_prob)
   end function
   pure elemental function pkw(q,alpha,beta,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,alpha,beta; logical,intent(in),optional::lower_tail,log_p; real(dp)::v
      v=pgkw_scalar(q,alpha,beta,1.0_dp,0.0_dp,1.0_dp,lower_tail,log_p)
   end function
   pure elemental function qkw(p,alpha,beta,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,alpha,beta; logical,intent(in),optional::lower_tail,log_p; real(dp)::v
      v=qgkw_scalar(p,alpha,beta,1.0_dp,0.0_dp,1.0_dp,lower_tail,log_p)
   end function
   function rkw(n,alpha,beta) result(x)
      integer,intent(in)::n; real(dp),intent(in)::alpha,beta; real(dp),allocatable::x(:); integer::i
      allocate(x(max(0,n))); do i=1,size(x); x(i)=rgkw_scalar(alpha,beta,1.0_dp,0.0_dp,1.0_dp); end do
   end function

   pure elemental function dbeta_(x,gamma,delta,log_prob) result(v)
      real(dp),intent(in)::x,gamma,delta; logical,intent(in),optional::log_prob; real(dp)::v
      v=dgkw_scalar(x,1.0_dp,1.0_dp,gamma,delta,1.0_dp,log_prob)
   end function
   pure elemental function pbeta_(q,gamma,delta,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,gamma,delta; logical,intent(in),optional::lower_tail,log_p; real(dp)::v
      v=pgkw_scalar(q,1.0_dp,1.0_dp,gamma,delta,1.0_dp,lower_tail,log_p)
   end function
   pure elemental function qbeta_(p,gamma,delta,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,gamma,delta; logical,intent(in),optional::lower_tail,log_p; real(dp)::v
      v=qgkw_scalar(p,1.0_dp,1.0_dp,gamma,delta,1.0_dp,lower_tail,log_p)
   end function
   function rbeta_(n,gamma,delta) result(x)
      integer,intent(in)::n; real(dp),intent(in)::gamma,delta; real(dp),allocatable::x(:); integer::i
      allocate(x(max(0,n))); do i=1,size(x); x(i)=rgkw_scalar(1.0_dp,1.0_dp,gamma,delta,1.0_dp); end do
   end function

   function llgkw(par,data) result(v)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::v
      v=family_nll(fam_gkw,par,data)
   end function llgkw

   function llbkw(par,data) result(v)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::v
      v=family_nll(fam_bkw,par,data)
   end function llbkw

   function llkkw(par,data) result(v)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::v
      v=family_nll(fam_kkw,par,data)
   end function llkkw

   function llekw(par,data) result(v)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::v
      v=family_nll(fam_ekw,par,data)
   end function llekw

   function llmc(par,data) result(v)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::v
      v=family_nll(fam_mc,par,data)
   end function llmc

   function llkw(par,data) result(v)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::v
      v=family_nll(fam_kw,par,data)
   end function llkw

   function llbeta(par,data) result(v)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::v
      v=family_nll(fam_beta,par,data)
   end function llbeta

   function grgkw(par,data) result(g)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::g(5),h(5,5),v
      call family_derivatives(fam_gkw,par,data,v,g,h)
   end function grgkw

   function grbkw(par,data) result(g)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::g(4),h(4,4),v
      call family_derivatives(fam_bkw,par,data,v,g,h)
   end function grbkw

   function grkkw(par,data) result(g)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::g(4),h(4,4),v
      call family_derivatives(fam_kkw,par,data,v,g,h)
   end function grkkw

   function grekw(par,data) result(g)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::g(3),h(3,3),v
      call family_derivatives(fam_ekw,par,data,v,g,h)
   end function grekw

   function grmc(par,data) result(g)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::g(3),h(3,3),v
      call family_derivatives(fam_mc,par,data,v,g,h)
   end function grmc

   function grkw(par,data) result(g)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::g(2),h(2,2),v
      call family_derivatives(fam_kw,par,data,v,g,h)
   end function grkw

   function grbeta(par,data) result(g)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::g(2),h(2,2),v
      call family_derivatives(fam_beta,par,data,v,g,h)
   end function grbeta

   function hsgkw(par,data) result(h)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::h(5,5),g(5),v
      call family_derivatives(fam_gkw,par,data,v,g,h)
   end function hsgkw

   function hsbkw(par,data) result(h)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::h(4,4),g(4),v
      call family_derivatives(fam_bkw,par,data,v,g,h)
   end function hsbkw

   function hskkw(par,data) result(h)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::h(4,4),g(4),v
      call family_derivatives(fam_kkw,par,data,v,g,h)
   end function hskkw

   function hsekw(par,data) result(h)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::h(3,3),g(3),v
      call family_derivatives(fam_ekw,par,data,v,g,h)
   end function hsekw

   function hsmc(par,data) result(h)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::h(3,3),g(3),v
      call family_derivatives(fam_mc,par,data,v,g,h)
   end function hsmc

   function hskw(par,data) result(h)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::h(2,2),g(2),v
      call family_derivatives(fam_kw,par,data,v,g,h)
   end function hskw

   function hsbeta(par,data) result(h)
      real(dp),intent(in)::par(:),data(:)
      real(dp)::h(2,2),g(2),v
      call family_derivatives(fam_beta,par,data,v,g,h)
   end function hsbeta

end module gkwdist_distributions
