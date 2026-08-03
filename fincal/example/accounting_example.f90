! SPDX-License-Identifier: GPL-2.0-or-later
program accounting_example
   use fincal
   implicit none
   type(inventory_result) :: result_value
   real(dp), allocatable :: depreciation(:)

   result_value = cogs(2.0_dp, 2.0_dp, [3.0_dp, 5.0_dp], [3.0_dp, 5.0_dp], 7.0_dp, 'LIFO')
   print '(a,f8.2)', 'LIFO COGS: ', result_value%cost_of_goods
   print '(a,f8.2)', 'ending inventory: ', result_value%ending_inventory

   depreciation = double_declining_balance(1200.0_dp, 200.0_dp, 5)
   print '(a,*(f10.2))', 'DDB schedule: ', depreciation
end program accounting_example
