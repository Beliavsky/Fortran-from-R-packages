! SPDX-License-Identifier: GPL-2.0-or-later
module dykstra
    use dykstra_kinds, only : dp
    use dykstra_solver, only : dykstra_result, dykstra_solve
    implicit none
    private
    public :: dp, dykstra_result, dykstra_solve
end module dykstra
