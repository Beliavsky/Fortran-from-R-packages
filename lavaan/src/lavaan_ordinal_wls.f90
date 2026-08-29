module lavaan_ordinal_wls
   use lavaan_kinds, only : dp
   use lavaan_ordinal, only : polychoric_matrix
   use lavaan_linalg, only : inverse_general
   implicit none
   private

   type, public :: ordinal_wls_result
      real(dp), allocatable :: correlation(:, :)
      real(dp), allocatable :: gamma(:, :), weight(:, :), dwls_weight(:)
      integer :: status = 0, n_jackknife = 0
   end type ordinal_wls_result

   public :: ordinal_wls_correlation_weights

contains

   subroutine ordinal_wls_correlation_weights(data,result)
      integer,intent(in)::data(:,:)
      type(ordinal_wls_result),intent(out)::result
      real(dp),allocatable::off(:),jk(:,:),meanjk(:),gamma_off(:,:),w_off(:,:),corj(:,:),ridge(:,:)
      integer,allocatable::offpos(:)
      integer::n,p,q,qo,i,j,r,info,m
      integer,allocatable::sub(:,:)

      n=size(data,1)
      p=size(data,2)
      q=p*(p+1)/2
      qo=p*(p-1)/2
      if(n<4 .or. p<2) then
      result%status=-1
      return
      end if
      call polychoric_matrix(data,result%correlation,info)
      if(info/=0) then
      result%status=info
      return
      end if
      allocate(off(qo),offpos(qo))
      call offdiag_vector(result%correlation,off,offpos)
      allocate(jk(n,qo))
      jk=0.0_dp
      m=0
      allocate(sub(n-1,p))
      do r=1,n
         if(r>1) sub(1:r-1,:)=data(1:r-1,:)
         if(r<n) sub(r:n-1,:)=data(r+1:n,:)
         call polychoric_matrix(sub,corj,info)
         if(info==0) then
            m=m+1
            call offdiag_vector(corj,jk(m,:),offpos)
         end if
      end do
      result%n_jackknife=m
      if(m<max(3,qo+1)) then
         result%status=-2
         return
      end if
      meanjk=sum(jk(1:m,:),dim=1)/real(m,dp)
      allocate(gamma_off(qo,qo))
      gamma_off=0.0_dp
      do r=1,m
         off=jk(r,:)-meanjk
         gamma_off=gamma_off+spread(off,2,qo)*spread(off,1,qo)
      end do
      ! Jackknife covariance estimates Var(theta); multiply by n for asymptotic Gamma.
      gamma_off=real(n,dp)*real(m-1,dp)/real(m,dp)*gamma_off
      allocate(ridge(qo,qo))
      ridge=gamma_off
      do i=1,qo
         ridge(i,i)=ridge(i,i)+1.0e-8_dp*max(1.0_dp,abs(ridge(i,i)))
      end do
      call inverse_general(ridge,w_off,info)
      if(info/=0) then
         w_off=0.0_dp
         do i=1,qo
            if(gamma_off(i,i)>1.0e-12_dp) w_off(i,i)=1.0_dp/gamma_off(i,i)
         end do
      end if
      allocate(result%gamma(q,q),result%weight(q,q),result%dwls_weight(q))
      result%gamma=0.0_dp
      result%weight=0.0_dp
      result%dwls_weight=0.0_dp
      do i=1,qo
         result%dwls_weight(offpos(i))=w_off(i,i)
         do j=1,qo
            result%gamma(offpos(i),offpos(j))=gamma_off(i,j)
            result%weight(offpos(i),offpos(j))=w_off(i,j)
         end do
      end do
      result%status=0
   contains
      subroutine offdiag_vector(cor,v,pos)
         real(dp),intent(in)::cor(:,:)
         real(dp),intent(out)::v(:)
         integer,intent(out)::pos(:)
         integer::a,b,kk,vpos
         kk=0
         vpos=0
         do b=1,size(cor,1)
            do a=b,size(cor,1)
               vpos=vpos+1
               if(a/=b) then
               kk=kk+1
               v(kk)=cor(a,b)
               pos(kk)=vpos
               end if
            end do
         end do
      end subroutine offdiag_vector
   end subroutine ordinal_wls_correlation_weights
end module lavaan_ordinal_wls
