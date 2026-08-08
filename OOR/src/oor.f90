! Upstream OOR license declaration: LGPL (version unspecified).
module oor
   use oor_kinds, only : dp
   use oor_interfaces, only : scalar_objective, vector_objective
   use oor_random, only : set_random_seed
   use oor_test_functions, only : guirland, sin1, difficult, difficult2, double_sine
   use oor_poo, only : poo, poo_result, poo_leaf, poo_tree
   use oor_stosoo, only : stosoo, stosoo_options, stosoo_result, soo_node, soo_level
   implicit none
   public
end module oor
