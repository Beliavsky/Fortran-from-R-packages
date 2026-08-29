! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_bivariate_fit
   use r_compat, only: dp
   use matrixdist_linalg
   use matrixdist_ph, only: ph_exit_rates
   use matrixdist_dph, only: dph_exit_probs
   use matrixdist_multivariate, only: bivph_density, bivdph_density
   implicit none
   private
   public :: bivph_loglik, emstep_bivph, bivdph_loglik, emstep_bivdph
contains
   function outer(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::c(size(a),size(b))
      integer::i,j
      do i=1,size(a)
      do j=1,size(b)
      c(i,j)=a(i)*b(j)
      end do
      end do
   end function

   function bivph_loglik(alpha,s11,s12,s22,obs,weight) result(ll)
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:),obs(:,:)
      real(dp),intent(in),optional::weight(:)
      real(dp)::ll,w,f
      integer::k
      ll=0.0_dp
      do k=1,size(obs,1)
         w=1.0_dp
         if(present(weight))w=weight(k)
         f=bivph_density(obs(k,1),obs(k,2),alpha,s11,s12,s22)
         if(f<=0.0_dp)then
         ll=-huge(1.0_dp)
         return
         end if
         ll=ll+w*log(f)
      end do
   end function

   subroutine emstep_bivph(alpha,s11,s12,s22,obs,weight)
      real(dp),intent(inout)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp),intent(in)::obs(:,:)
      real(dp),intent(in),optional::weight(:)
      integer::p1,p2,p,k,i,j
      real(dp)::w,total,density,aux,rate
      real(dp),allocatable::exitv(:),bmean(:),zmean(:),nmean(:,:)
      real(dp),allocatable::e1(:,:),e2(:,:),g1(:,:),g2(:,:),c1(:,:),c2(:,:)
      real(dp),allocatable::bm1(:,:),bm2(:,:),v1(:),v2(:),r1(:),r2(:)
      p1=size(s11,1)
      p2=size(s22,1)
      p=p1+p2
      allocate(bmean(p1),zmean(p),nmean(p,p+1),c1(p1,p1),c2(p2,p2))
      bmean=0.0_dp
      zmean=0.0_dp
      nmean=0.0_dp
      total=0.0_dp
      exitv=ph_exit_rates(s22)
      do k=1,size(obs,1)
         w=1.0_dp
         if(present(weight))w=weight(k)
         total=total+w
         e2=matrix_exponential(s22*obs(k,2))
         v2=matmul(e2,exitv)
         bm1=outer(matmul(s12,v2),alpha)
         g1=matrix_exponential(matrix_vanloan(s11,s11,bm1)*obs(k,1))
         e1=g1(1:p1,1:p1)
         c1=g1(1:p1,p1+1:2*p1)
         r1=matmul(transpose(e1),alpha)
         r2=matmul(transpose(s12),r1)
         bm2=outer(exitv,r2)
         g2=matrix_exponential(matrix_vanloan(s22,s22,bm2)*obs(k,2))
         c2=g2(1:p2,p2+1:2*p2)
         v1=matmul(e1,matmul(s12,v2))
         density=dot_product(alpha,v1)
         if(density<=0.0_dp)cycle
         ! r2 below becomes alpha^T e1 S12 e2 as a column transpose.
         r2=matmul(transpose(e2),matmul(transpose(s12),r1))
         do i=1,p1
            bmean(i)=bmean(i)+alpha(i)*v1(i)*w/density
            zmean(i)=zmean(i)+c1(i,i)*w/density
            do j=1,p1
            nmean(i,j)=nmean(i,j)+s11(i,j)*c1(j,i)*w/density
            end do
            do j=1,p2
               aux=v2(j)*r1(i)
               nmean(i,p1+j)=nmean(i,p1+j)+s12(i,j)*aux*w/density
            end do
         end do
         do i=1,p2
            zmean(p1+i)=zmean(p1+i)+c2(i,i)*w/density
            nmean(p1+i,p+1)=nmean(p1+i,p+1)+r2(i)*exitv(i)*w/density
            do j=1,p2
            nmean(p1+i,p1+j)=nmean(p1+i,p1+j)+s22(i,j)*c2(j,i)*w/density
            end do
         end do
      end do
      if(total<=0.0_dp)return
      alpha=max(0.0_dp,bmean/total)
      do i=1,p1
         if(zmean(i)<=0.0_dp)cycle
         s11(i,:)=0.0_dp
         s11(i,i)=0.0_dp
         do j=1,p1
            if(j/=i)then
            s11(i,j)=max(0.0_dp,nmean(i,j)/zmean(i))
            s11(i,i)=s11(i,i)-s11(i,j)
            end if
         end do
         do j=1,p2
         s12(i,j)=max(0.0_dp,nmean(i,p1+j)/zmean(i))
         s11(i,i)=s11(i,i)-s12(i,j)
         end do
      end do
      do i=1,p2
         if(zmean(p1+i)<=0.0_dp)cycle
         rate=max(0.0_dp,nmean(p1+i,p+1)/zmean(p1+i))
         s22(i,:)=0.0_dp
         s22(i,i)=-rate
         do j=1,p2
            if(j/=i)then
            s22(i,j)=max(0.0_dp,nmean(p1+i,p1+j)/zmean(p1+i))
            s22(i,i)=s22(i,i)-s22(i,j)
            end if
         end do
      end do
   end subroutine emstep_bivph

   function bivdph_loglik(alpha,s11,s12,s22,obs,weight) result(ll)
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      integer,intent(in)::obs(:,:)
      real(dp),intent(in),optional::weight(:)
      real(dp)::ll,w,f
      integer::k
      ll=0.0_dp
      do k=1,size(obs,1)
         w=1.0_dp
         if(present(weight))w=weight(k)
         f=bivdph_density(obs(k,1),obs(k,2),alpha,s11,s12,s22)
         if(f<=0.0_dp)then
         ll=-huge(1.0_dp)
         return
         end if
         ll=ll+w*log(f)
      end do
   end function

   subroutine emstep_bivdph(alpha,s11,s12,s22,obs,weight)
      real(dp),intent(inout)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      integer,intent(in)::obs(:,:)
      real(dp),intent(in),optional::weight(:)
      integer::p1,p2,p,k,i,j,m,k1,k2
      real(dp)::w,total,density,rowsum
      real(dp),allocatable::exitv(:),bmean(:),nmean(:,:),a0(:,:),a2(:,:),a3(:)
      real(dp),allocatable::p1m(:,:),p2m(:,:),c1(:,:),c2(:,:),left(:),right(:)
      p1=size(s11,1)
      p2=size(s22,1)
      p=p1+p2
      allocate(bmean(p1),nmean(p,p+1),c1(p1,p1),c2(p2,p2))
      bmean=0.0_dp
      nmean=0.0_dp
      total=0.0_dp
      exitv=dph_exit_probs(s22)
      do k=1,size(obs,1)
         k1=obs(k,1)
         k2=obs(k,2)
         if(k1<1 .or. k2<1)cycle
         w=1.0_dp
         if(present(weight))w=weight(k)
         total=total+w
         p1m=matrix_power(k1-1,s11)
         p2m=matrix_power(k2-1,s22)
         a0=matmul(p1m,matmul(s12,p2m))
         a3=matmul(transpose(a0),alpha)
         right=matmul(a0,exitv)
         density=dot_product(alpha,right)
         if(density<=0.0_dp)cycle
         left=matmul(transpose(p1m),alpha)
         do i=1,p1
         bmean(i)=bmean(i)+alpha(i)*right(i)*w/density
         end do
         if(k1>1)then
            c1=0.0_dp
            do m=0,k1-2
               c1=c1+matmul(matrix_power(k1-m-2,s11), &
                    matmul(s12,matmul(p2m,outer(exitv,matmul(transpose(matrix_power(m,s11)),alpha)))))
            end do
            do i=1,p1
            do j=1,p1
            nmean(i,j)=nmean(i,j)+s11(i,j)*c1(j,i)*w/density
            end do
            end do
         end if
         ! transitions from first to second block
         a2=outer(matmul(p2m,exitv),left)
         do i=1,p1
         do j=1,p2
         nmean(i,p1+j)=nmean(i,p1+j)+s12(i,j)*a2(j,i)*w/density
         end do
         end do
         do i=1,p2
         nmean(p1+i,p+1)=nmean(p1+i,p+1)+a3(i)*exitv(i)*w/density
         end do
         if(k2>1)then
            c2=0.0_dp
            do m=0,k2-2
               c2=c2+matmul(matrix_power(k2-m-2,s22), &
                    outer(exitv,matmul(transpose(matmul(p1m,matmul(s12,matrix_power(m,s22)))),alpha)))
            end do
            do i=1,p2
            do j=1,p2
            nmean(p1+i,p1+j)=nmean(p1+i,p1+j)+s22(i,j)*c2(j,i)*w/density
            end do
            end do
         end if
      end do
      if(total<=0.0_dp)return
      alpha=max(0.0_dp,bmean/total)
      do i=1,p1
         rowsum=sum(nmean(i,:))
         if(rowsum<=0.0_dp)cycle
         do j=1,p1
         s11(i,j)=max(0.0_dp,nmean(i,j)/rowsum)
         end do
         do j=1,p2
         s12(i,j)=max(0.0_dp,nmean(i,p1+j)/rowsum)
         end do
      end do
      do i=1,p2
         rowsum=sum(nmean(p1+i,:))
         if(rowsum<=0.0_dp)cycle
         do j=1,p2
         s22(i,j)=max(0.0_dp,nmean(p1+i,p1+j)/rowsum)
         end do
      end do
   end subroutine emstep_bivdph
end module matrixdist_bivariate_fit
