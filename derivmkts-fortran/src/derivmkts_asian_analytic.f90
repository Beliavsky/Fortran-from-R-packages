! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_asian_analytic
    use derivmkts_kinds, only: dp
    use derivmkts_black_scholes, only: bscall, bsput
    use derivmkts_types, only: option_pair
    implicit none
    private
    public :: geomavgprice, geomavgpricecall, geomavgpriceput
    public :: geomavgstrike, geomavgstrikecall, geomavgstrikeput
contains

    pure real(dp) function asian_sigma(v,m,continuous) result(siga)
        real(dp), intent(in) :: v
        integer, intent(in) :: m
        logical, intent(in) :: continuous
        if (continuous) then
            siga = v/sqrt(3.0_dp)
        else
            siga = v*sqrt(real((m+1)*(2*m+1),dp)/6.0_dp)/real(m,dp)
        end if
    end function asian_sigma

    pure real(dp) function asian_dividend(r,d,v,m,continuous) result(da)
        real(dp), intent(in) :: r,d,v
        integer, intent(in) :: m
        logical, intent(in) :: continuous
        if (continuous) then
            da = 0.5_dp*(r+d+v*v/6.0_dp)
        else
            da = 0.5_dp*(r*real(m-1,dp)/real(m,dp) + &
                (d+0.5_dp*v*v)*real(m+1,dp)/real(m,dp) - &
                (v/real(m,dp))**2*real((m+1)*(2*m+1),dp)/6.0_dp)
        end if
    end function asian_dividend

    pure real(dp) function asian_rho(m,continuous) result(rho)
        integer, intent(in) :: m
        logical, intent(in) :: continuous
        if (continuous) then
            rho = 0.5_dp*sqrt(3.0_dp)
        else
            rho = 0.5_dp*sqrt(6.0_dp*real(m+1,dp)/real(2*m+1,dp))
        end if
    end function asian_rho

    pure function geomavgprice(s,k,v,r,tt,d,m,continuous) result(out)
        real(dp), intent(in) :: s,k,v,r,tt,d
        integer, intent(in) :: m
        logical, intent(in), optional :: continuous
        type(option_pair) :: out
        logical :: cont
        cont = .false.
        if (present(continuous)) cont = continuous
        out%call = geomavgpricecall(s,k,v,r,tt,d,m,cont)
        out%put = geomavgpriceput(s,k,v,r,tt,d,m,cont)
    end function geomavgprice

    pure real(dp) function geomavgpricecall(s,k,v,r,tt,d,m,continuous) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d
        integer, intent(in) :: m
        logical, intent(in), optional :: continuous
        logical :: cont
        real(dp) :: siga, da
        cont = .false.
        if (present(continuous)) cont = continuous
        siga = asian_sigma(v,m,cont)
        da = asian_dividend(r,d,v,m,cont)
        price = bscall(s,k,siga,r,tt,da)
    end function geomavgpricecall

    pure real(dp) function geomavgpriceput(s,k,v,r,tt,d,m,continuous) result(price)
        real(dp), intent(in) :: s,k,v,r,tt,d
        integer, intent(in) :: m
        logical, intent(in), optional :: continuous
        logical :: cont
        real(dp) :: siga, da
        cont = .false.
        if (present(continuous)) cont = continuous
        siga = asian_sigma(v,m,cont)
        da = asian_dividend(r,d,v,m,cont)
        price = bsput(s,k,siga,r,tt,da)
    end function geomavgpriceput

    pure function geomavgstrike(s,km,v,r,tt,d,m,continuous) result(out)
        real(dp), intent(in) :: s,km,v,r,tt,d
        integer, intent(in) :: m
        logical, intent(in), optional :: continuous
        type(option_pair) :: out
        logical :: cont
        cont = .false.
        if (present(continuous)) cont = continuous
        out%call = geomavgstrikecall(s,km,v,r,tt,d,m,cont)
        out%put = geomavgstrikeput(s,km,v,r,tt,d,m,cont)
    end function geomavgstrike

    pure real(dp) function geomavgstrikecall(s,km,v,r,tt,d,m,continuous) result(price)
        real(dp), intent(in) :: s,km,v,r,tt,d
        integer, intent(in) :: m
        logical, intent(in), optional :: continuous
        logical :: cont
        real(dp) :: siga, da, rho, vol
        cont = .false.
        if (present(continuous)) cont = continuous
        siga = asian_sigma(v,m,cont)
        da = asian_dividend(r,d,v,m,cont)
        rho = asian_rho(m,cont)
        vol = sqrt(max(siga*siga+v*v-2.0_dp*rho*siga*v,0.0_dp))
        price = bscall(s,km,vol,da,tt,d)
    end function geomavgstrikecall

    pure real(dp) function geomavgstrikeput(s,km,v,r,tt,d,m,continuous) result(price)
        real(dp), intent(in) :: s,km,v,r,tt,d
        integer, intent(in) :: m
        logical, intent(in), optional :: continuous
        logical :: cont
        real(dp) :: siga, da, rho, vol
        cont = .false.
        if (present(continuous)) cont = continuous
        siga = asian_sigma(v,m,cont)
        da = asian_dividend(r,d,v,m,cont)
        rho = asian_rho(m,cont)
        vol = sqrt(max(siga*siga+v*v-2.0_dp*rho*siga*v,0.0_dp))
        price = bsput(s,km,vol,da,tt,d)
    end function geomavgstrikeput

end module derivmkts_asian_analytic
