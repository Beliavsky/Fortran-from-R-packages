module boot_glm_diag
    use boot_kinds, only : dp
    implicit none
    private
    public :: glm_diagnostics, gaussian_scale, gamma_scale
contains
    elemental real(dp) function gaussian_scale(deviance,df_residual) result(scale)
        real(dp),intent(in)::deviance
        integer,intent(in)::df_residual
        scale=sqrt(deviance/real(df_residual,dp))
    end function gaussian_scale

    real(dp) function gamma_scale(y,fitted,weights,df_residual) result(scale)
        real(dp),intent(in)::y(:),fitted(:),weights(:)
        integer,intent(in)::df_residual
        if(size(y)/=size(fitted).or.size(y)/=size(weights))error stop "gamma_scale: size mismatch"
        scale=sqrt(sum(weights*(y/fitted-1.0_dp)**2)/real(df_residual,dp))
    end function gamma_scale

    subroutine glm_diagnostics(deviance_resid,pearson_resid,hat,rank,scale,res,rd,rp,cook)
        real(dp),intent(in)::deviance_resid(:),pearson_resid(:),hat(:),scale
        integer,intent(in)::rank
        real(dp),intent(out)::res(size(hat)),rd(size(hat)),rp(size(hat)),cook(size(hat))
        real(dp)::dev(size(hat)),pear(size(hat))
        integer::i
        if(size(deviance_resid)/=size(hat).or.size(pearson_resid)/=size(hat))error stop "glm_diagnostics: mismatch"
        dev=deviance_resid/scale
        pear=pearson_resid/scale
        do i=1,size(hat)
            rp(i)=pear(i)/sqrt(max(tiny(1.0_dp),1.0_dp-hat(i)))
            rd(i)=dev(i)/sqrt(max(tiny(1.0_dp),1.0_dp-hat(i)))
            cook(i)=hat(i)*rp(i)**2/(max(tiny(1.0_dp),1.0_dp-hat(i))*real(rank,dp))
            res(i)=sign(1.0_dp,dev(i))*sqrt(dev(i)**2+hat(i)*rp(i)**2)
        end do
    end subroutine glm_diagnostics
end module boot_glm_diag
