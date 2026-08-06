! SPDX-License-Identifier: GPL-2.0-only
module tsgarch_benchmarks
  use ghyp_kinds, only : dp
  implicit none
  private
  public :: fcp_benchmark_data, laurent_benchmark_data, log_relative_error
contains
  subroutine fcp_benchmark_data(coefficient,se_h,se_opg,se_qmle)
    real(dp),intent(out)::coefficient(4),se_h(4),se_opg(4),se_qmle(4)
    coefficient=[-0.619041e-2_dp,0.107613e-1_dp,0.153134_dp,0.805974_dp]
    se_h=[0.846212e-2_dp,0.285271e-2_dp,0.265228e-1_dp,0.335527e-1_dp]
    se_opg=[0.843359e-2_dp,0.132298e-2_dp,0.139737e-1_dp,0.165604e-1_dp]
    se_qmle=[0.918935e-2_dp,0.649319e-2_dp,0.535317e-1_dp,0.724614e-1_dp]
  end subroutine fcp_benchmark_data

  subroutine laurent_benchmark_data(coefficient,se_h)
    real(dp),intent(out)::coefficient(6),se_h(6)
    coefficient=[0.04016_dp,0.04028_dp,0.15189_dp,0.46892_dp,0.84713_dp,1.33403_dp]
    se_h=[0.01408_dp,0.00558_dp,0.01188_dp,0.04969_dp,0.01096_dp,0.13814_dp]
  end subroutine laurent_benchmark_data

  elemental real(dp) function log_relative_error(x,benchmark) result(value)
    real(dp),intent(in)::x,benchmark
    if(abs(x-benchmark)<=epsilon(1.0_dp)*max(1.0_dp,abs(benchmark)))then
      value=huge(1.0_dp)
    else if(abs(benchmark)<=tiny(1.0_dp))then
      value=-huge(1.0_dp)
    else
      value=-log10(abs(x-benchmark)/abs(benchmark))
    end if
  end function log_relative_error
end module tsgarch_benchmarks
