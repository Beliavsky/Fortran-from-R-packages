program test_parity_targets
   use SpatialExtremes
   implicit none
   integer :: info,i
   integer,allocatable :: parts(:,:),sizes(:)
   real(dp) :: acond(2,3),data_cond(2),xcond(40,2)
   real(dp),allocatable :: lines(:,:),tbm(:,:),contrib(:,:)
   real(dp) :: coord2(3,2),data(8,3)
   real(dp) :: theta(2),ll
   type(composite_se_t) :: se
   type(dic_result_t) :: dic
   real(dp) :: chain_loc(3,2),chain_scale(3,2),chain_shape(3,2),dgev_data(4,2)

   call check(abs(stirling_second_kind(5,2)-15.0_dp)<1e-12_dp,'Stirling S(5,2)')
   call check(bell_number(4)==15,'Bell B4')
   call list_set_partitions(4,parts,sizes,info)
   call check(info==0.and.size(parts,2)==15,'partition enumeration')
   call check(all(sizes>=1).and.all(sizes<=4),'partition block counts')

   acond=reshape([0.5_dp,0.0_dp,0.0_dp,0.5_dp,0.5_dp,0.5_dp],[2,3])
   data_cond=[2.0_dp,3.0_dp]
   xcond=simulate_conditional_max_linear(data_cond,acond,acond,40,info)
   call check(info==0,'conditional max-linear status')
   do i=1,40
      call check(maxval(abs(xcond(i,:)-data_cond))<1e-11_dp,'conditional max-linear honors data')
   end do

   lines=van_der_corput_lines(20)
   call check(maxval(abs(sum(lines*lines,dim=2)-1.0_dp))<1e-12_dp,'Van der Corput unit lines')
   coord2=reshape([0.0_dp,0.5_dp,1.0_dp,0.0_dp,0.25_dp,0.8_dp],[3,2])
   tbm=simulate_gaussian_process_tbm(80,coord2,TBM_POWEREXP,0.1_dp,0.9_dp,1.2_dp,1.0_dp,80,.false.,info)
   call check(info==0.and.all(tbm==tbm).and.size(tbm,2)==3,'turning-band point simulation')
   tbm=simulate_gaussian_process_tbm(3,coord2,TBM_POWEREXP,0.0_dp,1.0_dp,1.0_dp,1.0_dp,20,.true.,info)
   call check(info==0.and.size(tbm,2)==9,'turning-band grid simulation')
   tbm=simulate_schlather_tbm(2,coord2,TBM_POWEREXP,0.0_dp,1.0_dp,1.0_dp,20,.false.,0.1_dp,info)
   call check(info==0.and.size(tbm,2)==3.and.all(tbm>0.0_dp),'Schlather turning-band max-stable simulation')
   tbm=simulate_geomgauss_tbm(2,coord2,TBM_POWEREXP,1.0_dp,0.0_dp,1.0_dp,1.0_dp,20,.false.,0.1_dp,info)
   call check(info==0.and.size(tbm,2)==3.and.all(tbm>0.0_dp),'geometric turning-band max-stable simulation')
   tbm=simulate_extremalt_tbm(2,coord2,TBM_POWEREXP,0.0_dp,1.0_dp,1.0_dp,2.0_dp,20,.false.,0.1_dp,info)
   call check(info==0.and.size(tbm,2)==3.and.all(tbm>=0.0_dp),'extremal-t turning-band max-stable simulation')

   theta=[0.2_dp,-0.1_dp]
   call composite_sandwich(test_contributions,theta,se)
   call check(se%info==0,'generic sandwich status')
   call check(all(se%stderr>0.0_dp),'generic sandwich positive SE')
   call check(maxval(abs(se%covariance-transpose(se%covariance)))<1e-10_dp,'sandwich covariance symmetry')

   data=reshape([1.1_dp,1.7_dp,2.2_dp,3.0_dp,1.3_dp,2.0_dp,2.8_dp,4.1_dp, &
                 1.4_dp,1.9_dp,2.5_dp,3.4_dp,1.2_dp,1.8_dp,2.9_dp,3.8_dp, &
                 1.6_dp,2.1_dp,2.6_dp,3.5_dp,1.5_dp,2.3_dp,3.1_dp,4.3_dp],[8,3])
   se=brownresnick_frechet_standard_errors(data,coord2,1.3_dp,1.0_dp)
   call check(se%info==0,'Brown-Resnick standard errors')
   call check(all(se%stderr>0.0_dp),'Brown-Resnick positive SE')
   contrib=brown_resnick_loglik_contributions(data,coord2,[0.0_dp,0.0_dp,0.0_dp], &
      [1.0_dp,1.0_dp,1.0_dp],[0.0_dp,0.0_dp,0.0_dp],1.3_dp,1.0_dp)
   ll=brown_resnick_loglik(data,coord2,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
      [0.0_dp,0.0_dp,0.0_dp],1.3_dp,1.0_dp)
   call check(abs(sum(contrib)-ll)<1e-10_dp,'full-margin contribution decomposition')

   dgev_data=reshape([0.0_dp,0.2_dp,-0.1_dp,0.3_dp,1.0_dp,1.2_dp,0.8_dp,1.1_dp],[4,2])
   chain_loc(:,1)=0.0_dp
   chain_loc(:,2)=1.0_dp
   chain_scale=1.0_dp
   chain_shape=0.0_dp
   dic=latent_dic(dgev_data,chain_loc,chain_scale,chain_shape)
   call check(abs(dic%effective_npar)<1e-12_dp,'latent DIC identical-chain effective npar')
   call check(abs(dic%dic-dic%dbar)<1e-12_dp,'latent DIC identical-chain identity')

   print *, 'test_parity_targets: PASS'
contains
   subroutine test_contributions(x,c)
      real(dp),intent(in)::x(:)
      real(dp),allocatable,intent(out)::c(:,:)
      real(dp)::b1(6),b2(6),q1(6),q2(6)
      integer::r,t,k
      b1=[1.0_dp,-0.4_dp,0.7_dp,1.4_dp,-1.1_dp,0.2_dp]
      b2=[0.3_dp,1.2_dp,-0.8_dp,0.5_dp,0.9_dp,-1.3_dp]
      q1=[0.2_dp,0.5_dp,0.7_dp,0.4_dp,0.9_dp,0.3_dp]
      q2=[0.6_dp,0.3_dp,0.4_dp,0.8_dp,0.2_dp,0.7_dp]
      allocate(c(3,2))
      r=0
      do t=1,3
         do k=1,2
            r=r+1
            c(t,k)=b1(r)*x(1)+b2(r)*x(2)-0.5_dp*(q1(r)*x(1)**2+q2(r)*x(2)**2)
         end do
      end do
   end subroutine test_contributions

   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
      print *, 'FAIL: ',trim(msg)
      error stop 1
      end if
   end subroutine check
end program test_parity_targets
