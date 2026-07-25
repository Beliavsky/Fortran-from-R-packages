! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2 of the License, or any later version.
module deoptimr
    use deoptimr_kinds, only: dp
    use deoptimr_interfaces, only: objective_function, constraint_function
    use deoptimr_rng, only: seed_rng
    use deoptimr_types, only: jde_control, ncde_control, de_result, ncde_result
    use deoptimr_jde, only: jde_optimize, spjde_optimize
    use deoptimr_ncde, only: ncde_optimize
    implicit none
    public
end module deoptimr
