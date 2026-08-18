! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet
  use degreenet_kinds
  use degreenet_rng
  use degreenet_math
  use degreenet_distributions
  use degreenet_compound
  use degreenet_models
  use degreenet_fit
  use degreenet_observation
  use degreenet_observed_fit
  use degreenet_diagnostics
  use degreenet_simulation
  use degreenet_graph
  implicit none
  public
  character(len=*), parameter :: statnet_attribution = &
    "Based on 'statnet' project software (statnet.org)."
contains
  subroutine degreenet_print_attribution()
    print '(a)', statnet_attribution
    print '(a)', 'For license and citation information see statnet.org/attribution'
  end subroutine degreenet_print_attribution
end module degreenet
