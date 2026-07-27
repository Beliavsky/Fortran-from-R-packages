! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of backtest-modern-fortran.
! It is free software: you may redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 2 of the License, or any later version.

module backtest_kinds
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)

end module backtest_kinds
