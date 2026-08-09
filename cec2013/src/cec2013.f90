! SPDX-License-Identifier: GPL-3.0-or-later
module cec2013
   use cec2013_kinds, only : dp
   use cec2013_data, only : cec2013_context, cec2013_dimension_supported, &
      CEC2013_OK, CEC2013_BAD_DIMENSION, CEC2013_IO_ERROR, CEC2013_BAD_PROBLEM, CEC2013_BAD_SHAPE
   use cec2013_functions, only : cec2013_evaluate, cec2013_evaluate_batch, cec2013_optimum_value
   implicit none
   public
end module cec2013
