program test_cholperm
 use truncated_normal
 implicit none
 real(dp) :: sigma(3,3), lo(3), hi(3), err
 type(cholperm_result) :: cp

 sigma = reshape([1.0_dp,0.3_dp,0.2_dp, &
                  0.3_dp,1.5_dp,0.4_dp, &
                  0.2_dp,0.4_dp,1.2_dp], [3,3])
 lo = [-1.0_dp,0.0_dp,-2.0_dp]
 hi = [2.0_dp,1.5_dp,0.5_dp]

 call cholperm(sigma,lo,hi,cp,'GGE')
 if(cp%status/=0) error stop 'GGE cholperm status'
 err=maxval(abs(matmul(cp%lmat,transpose(cp%lmat))-sigma(cp%perm,cp%perm)))
 if(err>2.0e-12_dp) error stop 'GGE cholperm reconstruction'

 call cholperm(sigma,lo,hi,cp,'GB')
 if(cp%status/=0) error stop 'GB cholperm status'
 err=maxval(abs(matmul(cp%lmat,transpose(cp%lmat))-sigma(cp%perm,cp%perm)))
 if(err>2.0e-12_dp) error stop 'GB cholperm reconstruction'

 print '(a)', 'test_cholperm: PASS'
end program test_cholperm
