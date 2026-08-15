! rmutil computational translation umbrella module
! Copyright (C) 1998-2001 J.K. Lindsey
! Copyright (C) 2026 OpenAI (modern Fortran translation)
! SPDX-License-Identifier: GPL-2.0-or-later
module rmutil
   use rmutil_kinds
   use rmutil_special
   use rmutil_integrate
   use rmutil_quadrature
   use rmutil_linalg
   use rmutil_ode
   use rmutil_continuous
   use rmutil_discrete
   use rmutil_covariate
   use rmutil_utils
   use rmutil_pkpd
   implicit none
   public
end module rmutil
