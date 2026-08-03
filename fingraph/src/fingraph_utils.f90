! SPDX-License-Identifier: GPL-3.0-only
module fingraph_utils
   use fingraph_kinds, only : dp
   use fingraph_types, only : graph_metrics
   use fingraph_linalg, only : frobenius_norm
   implicit none
   private
   public :: block_diag, block_diag2, block_diag3
   public :: relative_error, metrics, fscore, recall, specificity, fdr, npv, accuracy
   public :: pairwise_matrix_rownorm2, upper_view_vec, prial

   interface block_diag
      module procedure block_diag2
      module procedure block_diag3
   end interface block_diag
contains
   pure function block_diag2(a,b) result(c)
      real(dp), intent(in) :: a(:,:),b(:,:)
      real(dp), allocatable :: c(:,:)
      integer :: n1,n2
      n1=size(a,1); n2=size(b,1)
      if (size(a,2)/=n1 .or. size(b,2)/=n2) then
         allocate(c(0,0)); return
      end if
      allocate(c(n1+n2,n1+n2)); c=0.0_dp
      c(1:n1,1:n1)=a; c(n1+1:n1+n2,n1+1:n1+n2)=b
   end function block_diag2

   pure function block_diag3(a,b,d) result(c)
      real(dp), intent(in) :: a(:,:),b(:,:),d(:,:)
      real(dp), allocatable :: c(:,:)
      real(dp), allocatable :: ab(:,:)
      ab=block_diag2(a,b); c=block_diag2(ab,d)
   end function block_diag3

   pure function relative_error(west,wtrue) result(value)
      real(dp), intent(in) :: west(:,:),wtrue(:,:)
      real(dp) :: value,den
      den=frobenius_norm(wtrue)
      if (den<=tiny(1.0_dp)) then
         value=frobenius_norm(west-wtrue)
      else
         value=frobenius_norm(west-wtrue)/den
      end if
   end function relative_error

   pure function metrics(wtrue,west,eps) result(m)
      real(dp), intent(in) :: wtrue(:,:),west(:,:)
      real(dp), intent(in), optional :: eps
      type(graph_metrics) :: m
      real(dp) :: threshold,den
      logical :: edge,est_edge
      integer :: n,i,j
      threshold=1e-4_dp; if (present(eps)) threshold=eps
      n=min(size(wtrue,1),size(west,1))
      do i=1,n-1
         do j=i+1,n
            edge=abs(wtrue(i,j))>threshold
            est_edge=abs(west(i,j))>threshold
            if (edge .and. est_edge) then
               m%true_positive=m%true_positive+1.0_dp
            else if (.not.edge .and. est_edge) then
               m%false_positive=m%false_positive+1.0_dp
            else if (edge .and. .not.est_edge) then
               m%false_negative=m%false_negative+1.0_dp
            else
               m%true_negative=m%true_negative+1.0_dp
            end if
         end do
      end do
      den=2.0_dp*m%true_positive+m%false_negative+m%false_positive
      if (den>0.0_dp) m%fscore=2.0_dp*m%true_positive/den
      den=m%true_positive+m%false_negative
      if (den>0.0_dp) m%recall=m%true_positive/den
      den=m%true_negative+m%false_positive
      if (den>0.0_dp) m%specificity=m%true_negative/den
      den=m%true_positive+m%true_negative+m%false_positive+m%false_negative
      if (den>0.0_dp) m%accuracy=(m%true_positive+m%true_negative)/den
      den=m%true_negative+m%false_negative
      if (den>0.0_dp) m%npv=m%true_negative/den
      den=m%false_positive+m%true_positive
      if (den>0.0_dp) m%fdr=m%false_positive/den
   end function metrics

   pure function fscore(wtrue,west,eps) result(value)
      real(dp), intent(in) :: wtrue(:,:),west(:,:)
      real(dp), intent(in), optional :: eps
      real(dp) :: value
      type(graph_metrics) :: m
      m=metrics(wtrue,west,eps); value=m%fscore
   end function fscore
   pure function recall(wtrue,west,eps) result(value)
      real(dp), intent(in) :: wtrue(:,:),west(:,:)
      real(dp), intent(in), optional :: eps
      real(dp) :: value
      type(graph_metrics) :: m
      m=metrics(wtrue,west,eps); value=m%recall
   end function recall
   pure function specificity(wtrue,west,eps) result(value)
      real(dp), intent(in) :: wtrue(:,:),west(:,:)
      real(dp), intent(in), optional :: eps
      real(dp) :: value
      type(graph_metrics) :: m
      m=metrics(wtrue,west,eps); value=m%specificity
   end function specificity
   pure function fdr(wtrue,west,eps) result(value)
      real(dp), intent(in) :: wtrue(:,:),west(:,:)
      real(dp), intent(in), optional :: eps
      real(dp) :: value
      type(graph_metrics) :: m
      m=metrics(wtrue,west,eps); value=m%fdr
   end function fdr
   pure function npv(wtrue,west,eps) result(value)
      real(dp), intent(in) :: wtrue(:,:),west(:,:)
      real(dp), intent(in), optional :: eps
      real(dp) :: value
      type(graph_metrics) :: m
      m=metrics(wtrue,west,eps); value=m%npv
   end function npv
   pure function accuracy(wtrue,west,eps) result(value)
      real(dp), intent(in) :: wtrue(:,:),west(:,:)
      real(dp), intent(in), optional :: eps
      real(dp) :: value
      type(graph_metrics) :: m
      m=metrics(wtrue,west,eps); value=m%accuracy
   end function accuracy

   pure function pairwise_matrix_rownorm2(m) result(v)
      real(dp), intent(in) :: m(:,:)
      real(dp), allocatable :: v(:,:)
      integer :: n,i,j
      n=size(m,1); allocate(v(n,n)); v=0.0_dp
      do i=1,n-1
         do j=i+1,n
            v(i,j)=sum((m(i,:)-m(j,:))**2); v(j,i)=v(i,j)
         end do
      end do
   end function pairwise_matrix_rownorm2

   pure function upper_view_vec(m) result(v)
      real(dp), intent(in) :: m(:,:)
      real(dp), allocatable :: v(:)
      integer :: p,i,j,k
      p=size(m,1); allocate(v(p*(p-1)/2)); k=0
      do i=1,p-1
         do j=i+1,p
            k=k+1; v(k)=m(i,j)
         end do
      end do
   end function upper_view_vec

   pure function prial(ltrue,lest,lscm) result(value)
      real(dp), intent(in) :: ltrue(:,:),lest(:,:),lscm(:,:)
      real(dp) :: value,den
      den=frobenius_norm(lscm-ltrue)
      if (den<=tiny(1.0_dp)) then
         value=0.0_dp
      else
         value=100.0_dp*(1.0_dp-(frobenius_norm(lest-ltrue)/den)**2)
      end if
   end function prial
end module fingraph_utils
