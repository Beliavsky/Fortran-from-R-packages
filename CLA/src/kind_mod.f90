! SPDX-License-Identifier: GPL-3.0-or-later
module kind_mod
   ! dp follows the project convention kind(1.0d0)
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
end module kind_mod
