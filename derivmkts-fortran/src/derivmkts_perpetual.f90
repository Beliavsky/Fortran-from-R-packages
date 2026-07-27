! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_perpetual
    use derivmkts_kinds, only: dp
    use derivmkts_barriers, only: ur, dr
    use derivmkts_types, only: perpetual_result
    implicit none
    private
    public :: callperpetual, putperpetual
contains
    pure function callperpetual(s,k,v,r,d) result(out)
        real(dp), intent(in) :: s,k,v,r,d
        type(perpetual_result) :: out
        real(dp) :: g,h1
        if (d < 1.0e-13_dp .and. r > 0.0_dp) then
            out%price = s
            out%barrier = huge(1.0_dp)
            return
        end if
        g = sqrt(((r-d)/(v*v)-0.5_dp)**2+2.0_dp*r/(v*v))
        h1 = 0.5_dp-(r-d)/(v*v)+g
        out%barrier = k*h1/(h1-1.0_dp)
        out%price = (max(out%barrier,s)-k)*ur(s,v,r,1.0_dp,d,out%barrier,.true.)
    end function callperpetual

    pure function putperpetual(s,k,v,r,d) result(out)
        real(dp), intent(in) :: s,k,v,r,d
        type(perpetual_result) :: out
        real(dp) :: g,h2
        if (r < 1.0e-13_dp .and. d > 0.0_dp) then
            out%price = k
            out%barrier = 0.0_dp
            return
        end if
        g = sqrt(((r-d)/(v*v)-0.5_dp)**2+2.0_dp*r/(v*v))
        h2 = 0.5_dp-(r-d)/(v*v)-g
        out%barrier = k*h2/(h2-1.0_dp)
        out%price = (k-min(s,out%barrier))*dr(s,v,r,1.0_dp,d,out%barrier,.true.)
    end function putperpetual
end module derivmkts_perpetual
