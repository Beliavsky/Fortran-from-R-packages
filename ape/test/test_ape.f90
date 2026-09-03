program test_ape
   use ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use ape
   implicit none
   integer :: failures

   failures = 0
   call test_tree_operations(failures)
   call test_tree_comparison(failures)
   call test_split_collections(failures)
   call test_tree_editing(failures)
   call test_reconstruction(failures)
   call test_fastme(failures)
   call test_incomplete_reconstruction(failures)
   call test_completion(failures)
   call test_statistics(failures)
   call test_dna_models(failures)
   call test_graph_statistics(failures)
   call test_misc_statistics(failures)
   call test_ordination(failures)
   call test_continuous_ace(failures)
   call test_discrete_ace(failures)
   call test_pgls(failures)
   call test_chronopl(failures)
   call test_chronos_clock(failures)
   call test_chronos_models(failures)
   call test_compar_ou(failures)
   call test_compar_lynch(failures)
   call test_corphylo(failures)
   call test_binary_pglmm(failures)
   call test_reconstruct(failures)
   call test_skyline(failures)
   call test_birthdeath(failures)
   call test_birthdeath_extended(failures)

   if (failures /= 0) then
      error stop 'ape tests failed'
   end if
   print '(a)', 'All ape deterministic tests passed.'

contains

   subroutine test_tree_operations(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed tree-operation check.
      type(phylo_tree) :: tree
      type(phylo_tree) :: rescaled
      type(phylo_tree) :: dated
      type(phylo_tree) :: mpl_tree
      integer :: edge(6, 2)
      real(dp) :: edge_length(6)
      real(dp), allocatable :: depth(:)
      real(dp), allocatable :: node_distance(:, :)
      real(dp), allocatable :: contrasts(:)
      real(dp), allocatable :: variances(:)
      real(dp), allocatable :: ancestral(:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: p_value(:)
      integer, allocatable :: counts(:)
      integer, allocatable :: balance(:, :)
      integer, allocatable :: path(:)
      integer, allocatable :: topo_depth(:)
      real(dp), allocatable :: covariance(:, :)
      logical, allocatable :: descendant_matrix(:, :)
      integer :: info

      edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
      edge_length = 1.0_dp
      tree = make_phylo_tree(4, edge, edge_length)
      call check(tree%valid(), 'tree validity', failures)
      call check(tree%root() == 5, 'tree root', failures)
      call check(mrca(tree, 1, 2) == 6, 'MRCA siblings', failures)
      call check(mrca(tree, 1, 4) == 5, 'MRCA across clades', failures)

      call node_depth_edgelength(tree, depth, info)
      call check(info == 0, 'edge-length depth status', failures)
      call check_close(depth(1), 2.0_dp, 1.0e-12_dp, 'tip depth', failures)
      call check_close(depth(5), 0.0_dp, 1.0e-12_dp, 'root depth', failures)

      call descendant_tip_counts(tree, counts, info)
      call check(info == 0, 'descendant counts status', failures)
      call check(all(counts == [1, 1, 1, 1, 4, 2, 2]), 'descendant tip counts', failures)
      call check(cherry_count(tree) == 2, 'cherry count', failures)
      call check(is_binary_tree(tree), 'binary tree test', failures)
      call check(is_ultrametric_tree(tree), 'ultrametric tree test', failures)
      call check(is_monophyletic(tree, [1, 2]), 'monophyletic clade', failures)
      call check(.not. is_monophyletic(tree, [1, 3]), 'non-monophyletic tips', failures)

      call balance_counts(tree, balance, info)
      call check(info == 0, 'balance status', failures)
      call check(all(balance(1, :) == [2, 2]), 'root balance', failures)
      call check(all(balance(2, :) == [1, 1]), 'left balance', failures)

      call node_path(tree, 1, 4, path, info)
      call check(info == 0, 'node path status', failures)
      call check(all(path == [1, 6, 5, 7, 4]), 'node path values', failures)

      call dist_nodes(tree, node_distance, info)
      call check(info == 0, 'dist.nodes status', failures)
      call check_close(node_distance(1, 2), 2.0_dp, 1.0e-12_dp, 'sibling distance', failures)
      call check_close(node_distance(1, 4), 4.0_dp, 1.0e-12_dp, 'cross-clade distance', failures)

      call node_depth_count(tree, 2, topo_depth, info)
      call check(info == 0, 'topological depth status', failures)
      call check(all(topo_depth == [1, 1, 1, 1, 3, 2, 2]), 'topological depth values', failures)
      call phylogenetic_vcv(tree, covariance, info)
      call check(info == 0, 'phylogenetic VCV status', failures)
      call check_close(covariance(1, 1), 2.0_dp, 1.0e-12_dp, 'phylogenetic variance', failures)
      call check_close(covariance(1, 2), 1.0_dp, 1.0e-12_dp, 'phylogenetic covariance', failures)
      call check_close(covariance(1, 3), 0.0_dp, 1.0e-12_dp, 'cross-clade covariance', failures)
      call tip_descendant_matrix(tree, descendant_matrix, info)
      call check(info == 0, 'descendant matrix status', failures)
      call check(all(descendant_matrix(1, :)), 'root descendant matrix row', failures)
      call check(all(descendant_matrix(2, :) .eqv. [.true., .true., .false., .false.]), &
         'clade descendant matrix row', failures)

      call pic(tree, [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp], contrasts, variances, .true., rescaled, info)
      call check(info == 0, 'PIC status', failures)
      call check_close(contrasts(2), -sqrt(2.0_dp), 1.0e-12_dp, 'PIC left contrast', failures)
      call check_close(contrasts(3), -sqrt(2.0_dp), 1.0e-12_dp, 'PIC right contrast', failures)
      call check_close(contrasts(1), -4.0_dp / sqrt(3.0_dp), 1.0e-12_dp, 'PIC root contrast', failures)
      call check_close(variances(1), 3.0_dp, 1.0e-12_dp, 'PIC root variance', failures)

      call ace_pic(tree, [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp], ancestral, variances, .true., info)
      call check(info == 0, 'ace PIC status', failures)
      call check_close(ancestral(1), 4.0_dp, 1.0e-12_dp, 'ace PIC root estimate', failures)
      call check_close(ancestral(2), 2.0_dp, 1.0e-12_dp, 'ace PIC left estimate', failures)
      call check_close(ancestral(3), 6.0_dp, 1.0e-12_dp, 'ace PIC right estimate', failures)
      call check_close(variances(1), 3.0_dp, 1.0e-12_dp, 'ace PIC root variance', failures)

      mpl_tree = tree
      mpl_tree%edge_length = [4.0_dp, 6.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 2.0_dp]
      call chrono_mpl(mpl_tree, dated, info, standard_error, p_value)
      call check(info == 0, 'chronoMPL status', failures)
      call check(all(abs(dated%edge_length - [5.0_dp, 5.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp]) < 1.0e-12_dp), &
         'chronoMPL dated lengths', failures)
      call check_close(standard_error(1), sqrt(3.0_dp), 1.0e-12_dp, 'chronoMPL root SE', failures)
      call check_close(standard_error(2), 1.0_dp, 1.0e-12_dp, 'chronoMPL child SE', failures)
      call check_close(p_value(1), erfc(1.0_dp / sqrt(6.0_dp)), 1.0e-12_dp, 'chronoMPL root test', failures)
      call check_close(p_value(2), erfc(1.0_dp / sqrt(2.0_dp)), 1.0e-12_dp, 'chronoMPL child test', failures)
      call check_close(p_value(3), 1.0_dp, 1.0e-12_dp, 'chronoMPL equal-path test', failures)

      call compute_brtime(tree, [3.0_dp, 1.0_dp, 1.0_dp], dated, info)
      call check(info == 0, 'compute.brtime numeric status', failures)
      call check(all(abs(dated%edge_length - [2.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]) < 1.0e-12_dp), &
         'compute.brtime numeric lengths', failures)
   end subroutine test_tree_operations

   subroutine test_tree_comparison(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed tree-comparison check.
      type(phylo_tree) :: tree_a
      type(phylo_tree) :: tree_b
      integer :: edge_a(5, 2)
      integer :: edge_b(5, 2)
      real(dp) :: length_a(5)
      real(dp) :: length_b(5)
      real(dp) :: score
      integer :: distance
      integer :: info

      edge_a = reshape([5, 5, 5, 6, 6, 1, 2, 6, 3, 4], [5, 2])
      edge_b = reshape([5, 5, 5, 6, 6, 1, 3, 6, 2, 4], [5, 2])
      length_a = [0.5_dp, 0.5_dp, 1.0_dp, 0.5_dp, 0.5_dp]
      length_b = [0.5_dp, 0.5_dp, 2.0_dp, 0.5_dp, 0.5_dp]
      tree_a = make_phylo_tree(4, edge_a, length_a)
      tree_b = make_phylo_tree(4, edge_b, length_b)

      call topological_distance_ph85(tree_a, tree_a, distance, info)
      call check(info == 0 .and. distance == 0, 'PH85 identical-tree distance', failures)
      call topological_distance_ph85(tree_a, tree_b, distance, info)
      call check(info == 0 .and. distance == 2, 'PH85 quartet distance', failures)
      call branch_score_distance(tree_a, tree_a, score, info)
      call check_close(score, 0.0_dp, 1.0e-12_dp, 'branch score identical trees', failures)
      call branch_score_distance(tree_a, tree_b, score, info)
      call check_close(score, sqrt(5.0_dp), 1.0e-12_dp, 'branch score unmatched quartet splits', failures)
   end subroutine test_tree_comparison

   subroutine test_split_collections(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed split-collection check.
      type(phylo_tree) :: tree_a
      type(phylo_tree) :: tree_b
      type(phylo_tree) :: rooted
      type(phylo_tree) :: consensus
      type(phylo_tree) :: trees(3)
      type(phylo_tree) :: rooted_trees(2)
      type(split_collection) :: collection
      logical, allocatable :: split(:, :)
      real(dp), allocatable :: support(:)
      integer, allocatable :: counts(:)
      integer :: edge_a(5, 2)
      integer :: edge_b(5, 2)
      integer :: rooted_edge(6, 2)
      integer :: info
      integer :: i

      edge_a = reshape([5, 5, 5, 6, 6, 1, 2, 6, 3, 4], [5, 2])
      edge_b = reshape([5, 5, 5, 6, 6, 1, 3, 6, 2, 4], [5, 2])
      tree_a = make_phylo_tree(4, edge_a)
      tree_b = make_phylo_tree(4, edge_b)
      call tree_bipartitions(tree_a, split, info)
      call check(info == 0 .and. size(split, 1) == 1, 'quartet split count', failures)
      call check(all(split(1, :) .eqv. [.true., .true., .false., .false.]), 'quartet canonical split', failures)

      trees = [tree_a, tree_a, tree_b]
      call bitsplits(trees, collection, info)
      call check(info == 0 .and. collection%n_split == 2, 'bitsplits unique count', failures)
      call check(sum(collection%frequency) == 3, 'bitsplits total frequency', failures)
      call check(any(collection%frequency == 2), 'bitsplits repeated frequency', failures)
      call count_bipartitions(tree_a, trees, counts, info)
      call check(info == 0 .and. size(counts) == 1 .and. counts(1) == 2, 'count bipartitions', failures)

      call consensus_tree(trees, 0.5_dp, consensus, support, info)
      call check(info == 0 .and. consensus%valid(), 'majority consensus status', failures)
      call tree_bipartitions(consensus, split, info)
      call check(info == 0 .and. size(split, 1) == 1, 'majority consensus split count', failures)
      call check(all(split(1, :) .eqv. [.true., .true., .false., .false.]), &
         'majority consensus topology', failures)
      call check(size(support) == 2, 'majority consensus support count', failures)
      call check_close(support(1), 1.0_dp, 1.0e-12_dp, 'majority consensus root support', failures)
      call check_close(support(2), 2.0_dp / 3.0_dp, 1.0e-12_dp, 'majority consensus split support', failures)
      call consensus_tree(trees, 1.0_dp, consensus, support, info)
      call check(info == 0 .and. consensus%valid(), 'strict consensus status', failures)
      call tree_bipartitions(consensus, split, info)
      call check(info == 0 .and. size(split, 1) == 0, 'strict consensus star topology', failures)

      rooted_edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
      rooted = make_phylo_tree(4, rooted_edge)
      rooted_trees = [rooted, rooted]
      call prop_part(rooted_trees, collection, info)
      call check(info == 0 .and. collection%n_split == 3, 'prop.part unique clades', failures)
      call check(all(collection%frequency == 2), 'prop.part clade frequencies', failures)
      do i = 1, collection%n_split
         call check(count(collection%split(i, :)) == 2 .or. count(collection%split(i, :)) == 4, &
            'prop.part clade membership', failures)
      end do
      call prop_clades(rooted, rooted_trees, counts, info, rooted=.true.)
      call check(info == 0 .and. all(counts == 2), 'prop.clades rooted support', failures)
      call prop_clades(rooted, rooted_trees, counts, info)
      call check(info == 0 .and. all(counts == 2), 'prop.clades SHORTwise support', failures)
      call consensus_tree(rooted_trees, 1.0_dp, consensus, support, info, rooted=.true.)
      call check(info == 0 .and. consensus%valid(), 'rooted strict consensus status', failures)
      call prop_clades(consensus, rooted_trees, counts, info, rooted=.true.)
      call check(info == 0 .and. all(counts == 2), 'rooted strict consensus support', failures)
   end subroutine test_split_collections

   subroutine test_tree_editing(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed tree-editing check.
      type(phylo_tree) :: tree
      type(phylo_tree) :: pruned
      type(phylo_tree) :: clade
      type(phylo_tree) :: edited
      type(phylo_tree) :: star
      integer :: edge(8, 2)
      integer :: star_edge(4, 2)
      integer :: unary_edge(6, 2)
      real(dp) :: edge_length(8)
      real(dp) :: star_length(4)
      real(dp) :: unary_length(6)
      real(dp), allocatable :: original_distance(:, :)
      real(dp), allocatable :: edited_distance(:, :)
      integer, allocatable :: old_tips(:)
      integer :: info
      integer :: i
      integer :: j

      edge = reshape([6, 6, 7, 7, 8, 8, 9, 9, 7, 9, 1, 8, 2, 3, 4, 5], [8, 2])
      edge_length = 1.0_dp
      tree = make_phylo_tree(5, edge, edge_length)
      call check(tree%valid(), 'tree-edit source validity', failures)
      call dist_nodes(tree, original_distance, info)
      call check(info == 0, 'tree-edit source distances', failures)

      call drop_tips(tree, [3], pruned, old_tips, info)
      call check(info == 0 .and. pruned%valid(), 'drop tips status', failures)
      call check(all(old_tips == [1, 2, 4, 5]), 'drop tips mapping', failures)
      call check(.not. has_singles(pruned), 'drop tips singleton collapse', failures)
      call dist_nodes(pruned, edited_distance, info)
      call check(info == 0, 'drop tips distances status', failures)
      do i = 1, size(old_tips)
         do j = 1, size(old_tips)
            call check_close(edited_distance(i, j), original_distance(old_tips(i), old_tips(j)), 1.0e-12_dp, &
               'drop tips patristic preservation', failures)
         end do
      end do

      call keep_tips(tree, [5, 1, 4], pruned, old_tips, info)
      call check(info == 0 .and. all(old_tips == [1, 4, 5]), 'keep tips mapping', failures)
      call dist_nodes(pruned, edited_distance, info)
      do i = 1, size(old_tips)
         do j = 1, size(old_tips)
            call check_close(edited_distance(i, j), original_distance(old_tips(i), old_tips(j)), 1.0e-12_dp, &
               'keep tips patristic preservation', failures)
         end do
      end do

      call extract_clade(tree, 7, clade, old_tips, info)
      call check(info == 0 .and. clade%valid(), 'extract clade status', failures)
      call check(all(old_tips == [1, 2, 3]), 'extract clade mapping', failures)
      call dist_nodes(clade, edited_distance, info)
      do i = 1, size(old_tips)
         do j = 1, size(old_tips)
            call check_close(edited_distance(i, j), original_distance(old_tips(i), old_tips(j)), 1.0e-12_dp, &
               'extract clade patristic preservation', failures)
         end do
      end do

      call reroot_node(tree, 7, edited, info)
      call check(info == 0 .and. edited%valid(), 'reroot node status', failures)
      call check(edited%root() == 6, 'reroot node numbering', failures)
      call dist_nodes(edited, edited_distance, info)
      call check(maxval(abs(edited_distance(1:5, 1:5) - original_distance(1:5, 1:5))) < 1.0e-12_dp, &
         'reroot patristic preservation', failures)

      call unroot_tree(tree, edited, info)
      call check(info == 0 .and. edited%valid(), 'unroot status', failures)
      call check(.not. is_rooted_tree(edited), 'unroot structural status', failures)
      call check(edited%n_node == tree%n_node - 1, 'unroot node count', failures)
      call dist_nodes(edited, edited_distance, info)
      call check(maxval(abs(edited_distance(1:5, 1:5) - original_distance(1:5, 1:5))) < 1.0e-12_dp, &
         'unroot patristic preservation', failures)

      call root_outgroup(tree, [1], edited, info, .true.)
      call check(info == 0 .and. edited%valid(), 'single-outgroup root status', failures)
      call check(is_rooted_tree(edited), 'single-outgroup resolved root status', failures)
      call check(is_monophyletic(edited, [1]), 'single-outgroup monophyly', failures)
      call dist_nodes(edited, edited_distance, info)
      call check(maxval(abs(edited_distance(1:5, 1:5) - original_distance(1:5, 1:5))) < 1.0e-12_dp, &
         'single-outgroup rooting preserves distances', failures)

      call root_outgroup(tree, [4, 5], edited, info, .true.)
      call check(info == 0 .and. edited%valid(), 'multi-outgroup root status', failures)
      call check(is_rooted_tree(edited), 'multi-outgroup resolved root status', failures)
      call check(is_monophyletic(edited, [4, 5]), 'multi-outgroup monophyly', failures)
      call dist_nodes(edited, edited_distance, info)
      call check(maxval(abs(edited_distance(1:5, 1:5) - original_distance(1:5, 1:5))) < 1.0e-12_dp, &
         'multi-outgroup rooting preserves distances', failures)

      unary_edge = reshape([4, 5, 5, 6, 7, 7, 5, 1, 6, 7, 2, 3], [6, 2])
      unary_length = 1.0_dp
      edited = make_phylo_tree(3, unary_edge, unary_length)
      call check(has_singles(edited), 'has singles detection', failures)
      call dist_nodes(edited, original_distance, info)
      call collapse_singles(edited, pruned, info)
      call check(info == 0 .and. .not. has_singles(pruned), 'collapse singles status', failures)
      call dist_nodes(pruned, edited_distance, info)
      call check(maxval(abs(edited_distance(1:3, 1:3) - original_distance(1:3, 1:3))) < 1.0e-12_dp, &
         'collapse singles patristic preservation', failures)

      star_edge = reshape([5, 5, 5, 5, 1, 2, 3, 4], [4, 2])
      star_length = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      star = make_phylo_tree(4, star_edge, star_length)
      call dist_nodes(star, original_distance, info)
      call check(.not. is_rooted_tree(star), 'star unrooted structural status', failures)
      call multi2di(star, edited, info)
      call check(info == 0 .and. edited%valid(), 'multi2di status', failures)
      call check(is_binary_tree(edited), 'multi2di binary result', failures)
      call dist_nodes(edited, edited_distance, info)
      call check(maxval(abs(edited_distance(1:4, 1:4) - original_distance(1:4, 1:4))) < 1.0e-12_dp, &
         'multi2di patristic preservation', failures)
      call di2multi(edited, 1.0e-12_dp, pruned, info)
      call check(info == 0 .and. pruned%valid(), 'di2multi status', failures)
      call check(pruned%n_node == 1, 'di2multi contracted zero branches', failures)
      call dist_nodes(pruned, edited_distance, info)
      call check(maxval(abs(edited_distance(1:4, 1:4) - original_distance(1:4, 1:4))) < 1.0e-12_dp, &
         'di2multi zero-edge patristic preservation', failures)
   end subroutine test_tree_editing

   subroutine test_reconstruction(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed reconstruction check.
      real(dp) :: distance(5, 5)
      real(dp) :: variance_matrix(5, 5)
      real(dp), allocatable :: reconstructed(:, :)
      type(phylo_tree) :: tree
      integer :: info

      distance = reshape([ &
         0.0_dp, 5.0_dp, 9.0_dp, 9.0_dp, 8.0_dp, &
         5.0_dp, 0.0_dp, 10.0_dp, 10.0_dp, 9.0_dp, &
         9.0_dp, 10.0_dp, 0.0_dp, 8.0_dp, 7.0_dp, &
         9.0_dp, 10.0_dp, 8.0_dp, 0.0_dp, 3.0_dp, &
         8.0_dp, 9.0_dp, 7.0_dp, 3.0_dp, 0.0_dp], [5, 5])

      call nj(distance, tree, info)
      call check(info == 0 .and. tree%valid(), 'NJ reconstruction status', failures)
      call dist_nodes(tree, reconstructed, info)
      call check(info == 0, 'NJ patristic status', failures)
      call check(maxval(abs(reconstructed(1:5, 1:5) - distance)) < 1.0e-12_dp, &
         'NJ additive-matrix parity', failures)

      call bionj(distance, tree, info)
      call check(info == 0 .and. tree%valid(), 'BIONJ reconstruction status', failures)
      call dist_nodes(tree, reconstructed, info)
      call check(info == 0, 'BIONJ patristic status', failures)
      call check(maxval(abs(reconstructed(1:5, 1:5) - distance)) < 1.0e-12_dp, &
         'BIONJ additive-matrix parity', failures)


      call triang_mtd(distance, tree, info)
      call check(info == 0 .and. tree%valid(), 'triangMtd reconstruction status', failures)
      call dist_nodes(tree, reconstructed, info)
      call check(info == 0, 'triangMtd patristic status', failures)
      call check(maxval(abs(reconstructed(1:5, 1:5) - distance)) < 1.0e-12_dp, &
         'triangMtd additive-matrix parity', failures)

      variance_matrix = 1.0_dp
      variance_matrix(1, 1) = 0.0_dp
      variance_matrix(2, 2) = 0.0_dp
      variance_matrix(3, 3) = 0.0_dp
      variance_matrix(4, 4) = 0.0_dp
      variance_matrix(5, 5) = 0.0_dp
      call mvr(distance, variance_matrix, tree, info)
      call check(info == 0 .and. tree%valid(), 'MVR reconstruction status', failures)
      call dist_nodes(tree, reconstructed, info)
      call check(info == 0, 'MVR patristic status', failures)
      call check(maxval(abs(reconstructed(1:5, 1:5) - distance)) < 1.0e-12_dp, &
         'MVR unit-variance NJ parity', failures)


      call njs(distance, tree, info)
      call check(info == 0 .and. tree%valid(), 'NJ* complete reconstruction status', failures)
      call dist_nodes(tree, reconstructed, info)
      call check(info == 0, 'NJ* complete patristic status', failures)
      call check(maxval(abs(reconstructed(1:5, 1:5) - distance)) < 1.0e-12_dp, &
         'NJ* complete-matrix NJ parity', failures)

      distance(1, 3) = -1.0_dp
      distance(3, 1) = -1.0_dp
      call njs(distance, tree, info)
      call check(info == 0 .and. tree%valid(), 'NJ* incomplete reconstruction status', failures)
      call dist_nodes(tree, reconstructed, info)
      call check(info == 0, 'NJ* incomplete patristic status', failures)
      distance(1, 3) = 9.0_dp
      distance(3, 1) = 9.0_dp
      call check(maxval(abs(reconstructed(1:5, 1:5) - distance)) < 1.0e-12_dp, &
         'NJ* missing additive distance recovery', failures)


      distance(1, 3) = -1.0_dp
      distance(3, 1) = -1.0_dp
      call triang_mtds(distance, tree, info)
      call check(info == 0 .and. tree%valid(), 'triangMtds incomplete reconstruction status', failures)
      call dist_nodes(tree, reconstructed, info)
      distance(1, 3) = 9.0_dp
      distance(3, 1) = 9.0_dp
      call check(maxval(abs(reconstructed(1:5, 1:5) - distance)) < 1.0e-12_dp, &
         'triangMtds missing additive distance recovery', failures)
   end subroutine test_reconstruction

   subroutine test_incomplete_reconstruction(failures)
      integer, intent(inout) :: failures !! Running count incremented for every NJ*/BIONJ*/MVR* parity check.
      real(dp) :: distance(6, 6)
      real(dp) :: variance_matrix(6, 6)
      type(phylo_tree) :: tree
      integer :: i
      integer :: info
      integer :: j
      logical :: found
      real(dp) :: value

      distance = reshape([ &
         0.0_dp, 2.0_dp, 7.0_dp, -1.0_dp, 5.0_dp, 9.0_dp, &
         2.0_dp, 0.0_dp, 6.0_dp, 8.0_dp, 4.0_dp, -1.0_dp, &
         7.0_dp, 6.0_dp, 0.0_dp, 3.0_dp, 7.0_dp, 5.0_dp, &
         -1.0_dp, 8.0_dp, 3.0_dp, 0.0_dp, 6.0_dp, 4.0_dp, &
         5.0_dp, 4.0_dp, 7.0_dp, 6.0_dp, 0.0_dp, 3.0_dp, &
         9.0_dp, -1.0_dp, 5.0_dp, 4.0_dp, 3.0_dp, 0.0_dp], [6, 6])
      variance_matrix = 0.0_dp
      do i = 1, 5
         do j = i + 1, 6
            variance_matrix(i, j) = 0.7_dp + 0.17_dp * real(i, dp) + 0.11_dp * real(j, dp)
            variance_matrix(j, i) = variance_matrix(i, j)
         end do
      end do

      call bionjs(distance, tree, info)
      call check(info == 0 .and. tree%valid(), 'BIONJ* incomplete reconstruction status', failures)
      value = edge_length_to_child(tree, 9, found)
      call check(found, 'BIONJ* internal edge lookup', failures)
      call check_close(value, 1.6532258064516128_dp, 1.0e-12_dp, 'BIONJ* upstream C branch parity', failures)
      value = edge_length_to_child(tree, 4, found)
      call check_close(value, 1.0910362893177818_dp, 1.0e-12_dp, 'BIONJ* final branch parity', failures)

      call mvrs(distance, variance_matrix, tree, info)
      call check(info == 0 .and. tree%valid(), 'MVR* incomplete reconstruction status', failures)
      value = edge_length_to_child(tree, 10, found)
      call check(found, 'MVR* internal edge lookup', failures)
      call check_close(value, 2.2707982928129691_dp, 1.0e-12_dp, 'MVR* upstream C branch parity', failures)
      value = edge_length_to_child(tree, 8, found)
      call check_close(value, 1.0568424003927852_dp, 1.0e-12_dp, 'MVR* final branch parity', failures)
   end subroutine test_incomplete_reconstruction

   subroutine test_fastme(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed FastME parity check.
      type(phylo_tree) :: tree
      type(phylo_tree) :: reference
      real(dp) :: d6(6, 6)
      real(dp) :: d7(7, 7)
      real(dp), allocatable :: observed_distance(:, :)
      real(dp), allocatable :: reference_distance(:, :)
      integer :: edge6(9, 2)
      integer :: edge7(11, 2)
      real(dp) :: length6(9)
      integer :: info
      integer :: topo_distance

      d6 = reshape([ &
         0.0_dp, 5.0_dp, 9.0_dp, 9.0_dp, 8.0_dp, 7.0_dp, &
         5.0_dp, 0.0_dp, 10.0_dp, 10.0_dp, 9.0_dp, 8.0_dp, &
         9.0_dp, 10.0_dp, 0.0_dp, 8.0_dp, 7.0_dp, 9.0_dp, &
         9.0_dp, 10.0_dp, 8.0_dp, 0.0_dp, 3.0_dp, 6.0_dp, &
         8.0_dp, 9.0_dp, 7.0_dp, 3.0_dp, 0.0_dp, 5.0_dp, &
         7.0_dp, 8.0_dp, 9.0_dp, 6.0_dp, 5.0_dp, 0.0_dp], [6, 6])
      edge6(:, 1) = [7, 7, 7, 8, 9, 9, 10, 10, 8]
      edge6(:, 2) = [1, 2, 8, 9, 3, 10, 4, 5, 6]

      call fastme_ols(d6, tree, info, .false.)
      call check(info == 0 .and. tree%valid(), 'FastME OLS insertion status', failures)
      length6 = [2.0_dp, 3.0_dp, 2.5_dp, 2.0_dp / 3.0_dp, 4.5_dp, 1.5_dp, 2.0_dp, 1.0_dp, 2.5_dp]
      reference = make_phylo_tree(6, edge6, length6)
      call dist_nodes(tree, observed_distance, info)
      call dist_nodes(reference, reference_distance, info)
      call check(maxval(abs(observed_distance(1:6, 1:6) - reference_distance(1:6, 1:6))) < 1.0e-11_dp, &
         'FastME OLS branch-length parity', failures)

      call fastme_bal(d6, tree, info, .false., .false.)
      call check(info == 0 .and. tree%valid(), 'FastME BME insertion status', failures)
      length6 = [2.0_dp, 3.0_dp, 2.25_dp, 0.75_dp, 4.75_dp, 1.25_dp, 2.0_dp, 1.0_dp, 2.75_dp]
      reference = make_phylo_tree(6, edge6, length6)
      call dist_nodes(tree, observed_distance, info)
      call dist_nodes(reference, reference_distance, info)
      call check(maxval(abs(observed_distance(1:6, 1:6) - reference_distance(1:6, 1:6))) < 1.0e-11_dp, &
         'FastME balanced branch-length parity', failures)

      d7 = reshape([ &
         0.0_dp, 12.792331380261_dp, 18.519955721532_dp, 17.221386350524_dp, 10.487526584573_dp, &
            14.782352437918_dp, 15.039955783440_dp, &
         12.792331380261_dp, 0.0_dp, 7.678567398831_dp, 9.945547669852_dp, 4.653896686348_dp, &
            14.362303476978_dp, 10.565172080548_dp, &
         18.519955721532_dp, 7.678567398831_dp, 0.0_dp, 9.848775342027_dp, 3.308133362474_dp, &
            10.039560921698_dp, 12.418077900912_dp, &
         17.221386350524_dp, 9.945547669852_dp, 9.848775342027_dp, 0.0_dp, 12.993711461408_dp, &
            8.223811952507_dp, 10.343109141175_dp, &
         10.487526584573_dp, 4.653896686348_dp, 3.308133362474_dp, 12.993711461408_dp, 0.0_dp, &
            14.550626493169_dp, 4.126728243822_dp, &
         14.782352437918_dp, 14.362303476978_dp, 10.039560921698_dp, 8.223811952507_dp, 14.550626493169_dp, &
            0.0_dp, 8.251579413600_dp, &
         15.039955783440_dp, 10.565172080548_dp, 12.418077900912_dp, 10.343109141175_dp, 4.126728243822_dp, &
            8.251579413600_dp, 0.0_dp], [7, 7])

      call fastme_ols(d7, tree, info, .false.)
      call check(info == 0, 'FastME OLS nonadditive insertion status', failures)
      edge7(:, 1) = [8, 8, 9, 9, 10, 10, 11, 12, 12, 11, 8]
      edge7(:, 2) = [1, 9, 2, 10, 3, 11, 12, 4, 6, 7, 5]
      reference = make_phylo_tree(7, edge7)
      call topological_distance_ph85(tree, reference, topo_distance, info)
      call check(info == 0 .and. topo_distance == 0, 'FastME OLS insertion topology parity', failures)
      call check_close(sum(tree%edge_length), 34.2517579607373_dp, 2.0e-11_dp, &
         'FastME OLS insertion total length', failures)

      call fastme_ols(d7, tree, info, .true.)
      call check(info == 0, 'FastME OLS NNI status', failures)
      edge7(:, 1) = [8, 8, 9, 10, 10, 9, 8, 11, 11, 12, 12]
      edge7(:, 2) = [1, 9, 10, 4, 6, 7, 11, 2, 12, 3, 5]
      reference = make_phylo_tree(7, edge7)
      call topological_distance_ph85(tree, reference, topo_distance, info)
      call check(info == 0 .and. topo_distance == 0, 'FastME OLS NNI topology parity', failures)
      call check_close(sum(tree%edge_length), 33.44048497130686_dp, 2.0e-11_dp, &
         'FastME OLS NNI total length', failures)

      call fastme_bal(d7, tree, info, .false., .false.)
      call check(info == 0, 'FastME BME nonadditive insertion status', failures)
      edge7(:, 1) = [8, 8, 9, 9, 10, 10, 11, 11, 8, 12, 12]
      edge7(:, 2) = [1, 9, 2, 10, 3, 11, 4, 6, 12, 5, 7]
      reference = make_phylo_tree(7, edge7)
      call topological_distance_ph85(tree, reference, topo_distance, info)
      call check(info == 0 .and. topo_distance == 0, 'FastME BME insertion topology parity', failures)
      call check_close(sum(tree%edge_length), 34.327775361935295_dp, 2.0e-11_dp, &
         'FastME BME insertion total length', failures)

      call fastme_bal(d7, tree, info, .true., .false.)
      call check(info == 0, 'FastME balanced NNI status', failures)
      edge7(:, 1) = [8, 8, 9, 10, 10, 9, 11, 11, 12, 12, 8]
      edge7(:, 2) = [1, 9, 10, 5, 7, 11, 3, 12, 4, 6, 2]
      reference = make_phylo_tree(7, edge7)
      call topological_distance_ph85(tree, reference, topo_distance, info)
      call check(info == 0 .and. topo_distance == 0, 'FastME balanced NNI topology parity', failures)
      call check_close(sum(tree%edge_length), 34.28059407184715_dp, 2.0e-11_dp, &
         'FastME balanced NNI total length', failures)

      call fastme_bal(d7, tree, info, .true., .true.)
      call check(info == 0, 'FastME balanced SPR status', failures)
      edge7(:, 1) = [8, 8, 9, 10, 10, 9, 8, 11, 12, 12, 11]
      edge7(:, 2) = [1, 9, 10, 5, 3, 2, 11, 12, 4, 6, 7]
      reference = make_phylo_tree(7, edge7)
      call topological_distance_ph85(tree, reference, topo_distance, info)
      call check(info == 0 .and. topo_distance == 0, 'FastME balanced SPR topology parity', failures)
      call check_close(sum(tree%edge_length), 33.43844390533529_dp, 2.0e-11_dp, &
         'FastME balanced SPR total length', failures)
   end subroutine test_fastme

   subroutine test_completion(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed missing-distance check.
      real(dp) :: distance(4, 4)
      real(dp), allocatable :: completed(:, :)
      real(dp) :: nan_value
      integer :: unresolved

      distance = reshape([ &
         0.0_dp, 2.0_dp, 4.0_dp, 4.0_dp, &
         2.0_dp, 0.0_dp, 4.0_dp, 4.0_dp, &
         4.0_dp, 4.0_dp, 0.0_dp, 2.0_dp, &
         4.0_dp, 4.0_dp, 2.0_dp, 0.0_dp], [4, 4])
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      distance(1, 3) = nan_value
      distance(3, 1) = nan_value

      call ultrametric_completion(distance, completed, unresolved)
      call check(unresolved == 0, 'ultrametric completion status', failures)
      call check_close(completed(1, 3), 4.0_dp, 1.0e-12_dp, 'ultrametric completion value', failures)

      call additive_completion(distance, completed, unresolved)
      call check(unresolved == 0, 'additive completion status', failures)
      call check_close(completed(1, 3), 4.0_dp, 1.0e-12_dp, 'additive completion value', failures)
   end subroutine test_completion

   subroutine test_statistics(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed phylogenetic-statistics check.
      type(phylo_tree) :: tree
      type(yule_result) :: fit
      integer :: edge(6, 2)
      real(dp) :: edge_length(6)
      real(dp), allocatable :: times(:)
      real(dp), allocatable :: widths(:)
      integer, allocatable :: lineages(:)
      real(dp) :: gamma
      real(dp) :: total_depth
      integer :: info

      edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
      edge_length = 1.0_dp
      tree = make_phylo_tree(4, edge, edge_length)

      call branching_times(tree, times, info)
      call check(info == 0, 'branching times status', failures)
      call check(all(abs(times - [2.0_dp, 1.0_dp, 1.0_dp]) < 1.0e-12_dp), &
         'branching times values', failures)

      call gamma_stat(tree, gamma, info)
      call check(info == 0, 'gamma statistic status', failures)
      call check_close(gamma, -sqrt(2.0_dp / 3.0_dp), 1.0e-12_dp, 'gamma statistic value', failures)

      call yule_fit(tree, fit, info)
      call check(info == 0, 'Yule fit status', failures)
      call check_close(fit%lambda, 1.0_dp / 3.0_dp, 1.0e-12_dp, 'Yule lambda', failures)
      call check_close(fit%standard_error, 1.0_dp / (3.0_dp * sqrt(2.0_dp)), 1.0e-12_dp, &
         'Yule standard error', failures)
      call check_close(fit%log_likelihood, -2.0_dp + log(2.0_dp / 3.0_dp), 1.0e-12_dp, &
         'Yule log likelihood', failures)

      call coalescent_intervals(tree, lineages, widths, total_depth, info)
      call check(info == 0, 'coalescent intervals status', failures)
      call check(all(lineages == [4, 3, 2]), 'coalescent lineages', failures)
      call check(all(abs(widths - [1.0_dp, 0.0_dp, 1.0_dp]) < 1.0e-12_dp), &
         'coalescent widths', failures)
      call check_close(total_depth, 2.0_dp, 1.0e-12_dp, 'coalescent total depth', failures)

      call ltt_coordinates(tree, times, lineages, info)
      call check(info == 0, 'LTT coordinate status', failures)
      call check(all(abs(times - [-2.0_dp, -1.0_dp, -1.0_dp, 0.0_dp]) < 1.0e-12_dp), &
         'LTT backward times', failures)
      call check(all(lineages == [1, 2, 3, 4]), 'LTT lineage counts', failures)
      call ltt_coordinates(tree, times, lineages, info, backward=.false., step_type='s')
      call check(info == 0, 'LTT forward lower-step status', failures)
      call check(all(abs(times - [0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]) < 1.0e-12_dp), &
         'LTT forward times', failures)
      call check(all(lineages == [2, 3, 4, 4]), 'LTT lower-step lineage counts', failures)
   end subroutine test_statistics

   subroutine test_dna_models(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed nucleotide-distance check.
      integer :: a(4)
      integer :: b(4)
      integer :: long_a(20)
      integer :: long_b(20)
      integer :: variance_b(20)
      integer :: alignment(3, 6)
      integer :: alignment2(2, 20)
      integer :: alignment3(3, 20)
      integer :: ambiguity_alignment(3, 3)
      integer :: terminal_alignment(1, 7)
      integer, allocatable :: converted_alignment(:, :)
      real(dp) :: state_proportions(17)
      real(dp) :: base_frequency(4)
      real(dp) :: distance
      real(dp) :: variance_value
      real(dp), allocatable :: matrix_distance(:, :)
      real(dp), allocatable :: matrix_variance(:, :)
      real(dp), allocatable :: table(:, :)
      logical, allocatable :: keep(:)
      logical, allocatable :: segregating(:)
      integer, allocatable :: positions(:)
      character(len=1), allocatable :: amino(:)
      integer :: comparable
      integer :: info
      integer :: i

      a = [dna_a, dna_c, dna_g, dna_t]
      b = [dna_a, dna_t, dna_g, dna_t]
      call dna_distance(a, b, 'RAW', distance, comparable, info)
      call check(info == 0 .and. comparable == 4, 'RAW DNA status', failures)
      call check_close(distance, 0.25_dp, 1.0e-12_dp, 'RAW DNA distance', failures)

      call dna_distance(a, b, 'JC69', distance, comparable, info)
      call check(info == 0, 'JC69 DNA status', failures)
      call check_close(distance, -0.75_dp * log(2.0_dp / 3.0_dp), 1.0e-12_dp, &
         'JC69 DNA distance', failures)

      call dna_distance(a, b, 'K80', distance, comparable, info)
      call check(info == 0, 'K80 DNA status', failures)
      call check_close(distance, -0.5_dp * log(0.5_dp), 1.0e-12_dp, 'K80 DNA distance', failures)
      call dna_distance(a, b, 'TS', distance, comparable, info)
      call check_close(distance, 1.0_dp, 1.0e-12_dp, 'TS DNA count', failures)
      call dna_distance(a, b, 'TV', distance, comparable, info)
      call check_close(distance, 0.0_dp, 1.0e-12_dp, 'TV DNA count', failures)

      do i = 1, 5
         long_a(4 * i - 3:4 * i) = [dna_a, dna_c, dna_g, dna_t]
      end do
      long_b = long_a
      long_b(1) = dna_g
      base_frequency = 0.25_dp
      call dna_distance(long_a, long_b, 'F81', distance, comparable, info, base_frequency=base_frequency)
      call check(info == 0, 'F81 DNA status', failures)
      call check_close(distance, -0.75_dp * log(14.0_dp / 15.0_dp), 1.0e-12_dp, &
         'F81 DNA distance', failures)
      call dna_distance(long_a, long_b, 'K81', distance, comparable, info)
      call check(info == 0, 'K81 DNA status', failures)
      call check_close(distance, -0.5_dp * log(0.9_dp), 1.0e-12_dp, 'K81 DNA distance', failures)
      call dna_distance(long_a, long_b, 'F84', distance, comparable, info, base_frequency=base_frequency)
      call check(info == 0, 'F84 DNA status', failures)
      call check_close(distance, -0.5_dp * log(0.9_dp), 1.0e-12_dp, 'F84 DNA distance', failures)
      call dna_distance(long_a, long_b, 'T92', distance, comparable, info, base_frequency=base_frequency)
      call check(info == 0, 'T92 DNA status', failures)
      call check_close(distance, -0.5_dp * log(0.9_dp), 1.0e-12_dp, 'T92 DNA distance', failures)
      call dna_distance(long_a, long_b, 'TN93', distance, comparable, info, base_frequency=base_frequency)
      call check(info == 0, 'TN93 DNA status', failures)
      call check_close(distance, -0.25_dp * log(0.8_dp), 1.0e-12_dp, 'TN93 DNA distance', failures)
      call dna_distance(long_a, long_b, 'LOGDET', distance, comparable, info)
      call check(info == 0 .and. distance > 0.0_dp, 'LOGDET DNA distance', failures)
      call dna_distance(long_a, long_b, 'PARALIN', distance, comparable, info)
      call check(info == 0 .and. distance > 0.0_dp, 'PARALIN DNA distance', failures)

      variance_b = long_a
      variance_b(1) = dna_g
      variance_b(2) = dna_t
      variance_b(3) = dna_c
      variance_b(4) = dna_g
      call dna_distance_with_variance(long_a, variance_b, 'JC69', distance, variance_value, comparable, info)
      call check(info == 0 .and. comparable == 20, 'JC69 variance status', failures)
      call check_close(variance_value, 0.01487603305785124_dp, 1.0e-14_dp, 'JC69 variance', failures)
      call dna_distance_with_variance(long_a, variance_b, 'JC69', distance, variance_value, comparable, info, &
         gamma_shape=2.0_dp)
      call check(info == 0, 'JC69 gamma variance status', failures)
      call check_close(variance_value, 0.006505653639058279_dp, 1.0e-14_dp, 'JC69 gamma variance', failures)
      call dna_distance_with_variance(long_a, variance_b, 'K80', distance, variance_value, comparable, info)
      call check(info == 0, 'K80 variance status', failures)
      call check_close(variance_value, 0.015341996173469385_dp, 1.0e-14_dp, 'K80 variance', failures)
      call dna_distance_with_variance(long_a, variance_b, 'K80', distance, variance_value, comparable, info, &
         gamma_shape=2.0_dp)
      call check(info == 0, 'K80 gamma variance status', failures)
      call check_close(variance_value, 0.021314913491980753_dp, 1.0e-14_dp, 'K80 gamma variance', failures)
      call dna_distance_with_variance(long_a, variance_b, 'F81', distance, variance_value, comparable, info, &
         base_frequency=base_frequency)
      call check(info == 0, 'F81 variance status', failures)
      call check_close(variance_value, 0.01487603305785124_dp, 1.0e-14_dp, 'F81 variance', failures)
      call dna_distance_with_variance(long_a, variance_b, 'K81', distance, variance_value, comparable, info)
      call check(info == 0, 'K81 variance status', failures)
      call check_close(variance_value, 0.15341996173469385_dp, 1.0e-13_dp, 'K81 upstream variance', failures)
      call dna_distance_with_variance(long_a, variance_b, 'F84', distance, variance_value, comparable, info, &
         base_frequency=base_frequency)
      call check(info == 0, 'F84 variance status', failures)
      call check_close(variance_value, 0.015341996173469385_dp, 1.0e-14_dp, 'F84 variance', failures)
      call dna_distance_with_variance(long_a, variance_b, 'T92', distance, variance_value, comparable, info, &
         base_frequency=base_frequency)
      call check(info == 0, 'T92 variance status', failures)
      call check_close(variance_value, 0.015341996173469385_dp, 1.0e-14_dp, 'T92 variance', failures)
      call dna_distance_with_variance(long_a, variance_b, 'TN93', distance, variance_value, comparable, info, &
         base_frequency=base_frequency)
      call check(info == 0, 'TN93 variance status', failures)
      call check_close(variance_value, 0.015341996173469385_dp, 1.0e-14_dp, 'TN93 variance', failures)
      call dna_distance_with_variance(long_a, variance_b, 'TN93', distance, variance_value, comparable, info, &
         gamma_shape=2.0_dp, base_frequency=base_frequency)
      call check(info == 0, 'TN93 gamma variance status', failures)
      call check_close(variance_value, 0.021314913491980753_dp, 1.0e-14_dp, 'TN93 gamma variance', failures)
      call dna_distance_with_variance(long_a, variance_b, 'LOGDET', distance, variance_value, comparable, info)
      call check(info == 0, 'LOGDET variance status', failures)
      call check_close(variance_value, 0.011409023668639051_dp, 1.0e-13_dp, 'LOGDET upstream variance', failures)
      call dna_distance_with_variance(long_a, variance_b, 'PARALIN', distance, variance_value, comparable, info)
      call check(info == 0, 'PARALIN variance status', failures)
      call check_close(variance_value, -0.14013622222761954_dp, 1.0e-13_dp, 'PARALIN upstream variance', failures)

      alignment(1, :) = [dna_a, dna_a, dna_a, dna_t, dna_t, dna_t]
      alignment(2, :) = [dna_a, dna_a, dna_g, dna_t, dna_t, dna_t]
      alignment(3, :) = [dna_a, dna_unknown, dna_g, dna_t, dna_t, dna_t]
      call dna_base_frequencies(alignment, base_frequency)
      call check_close(sum(base_frequency), 1.0_dp, 1.0e-12_dp, 'DNA base-frequency normalization', failures)
      call dna_global_deletion_mask(alignment, keep)
      call check(count(keep) == 5, 'DNA global deletion mask', failures)
      call dna_segregating_sites(alignment, segregating)
      call check(segregating(3), 'DNA segregating site', failures)
      call dna_pattern_positions([dna_a, dna_t, dna_g, dna_a, dna_t, dna_g], &
         [dna_a, dna_t, dna_g], positions)
      call check(all(positions == [1, 4]), 'DNA pattern positions', failures)
      call translate_dna([dna_a, dna_t, dna_g, dna_g, dna_c, dna_c], amino, info)
      call check(info == 0 .and. all(amino == ['M', 'A']), 'DNA translation', failures)
      call translate_dna([dna_a, dna_a, dna_r, dna_t, dna_g, dna_a], amino, info)
      call check(info == 0 .and. all(amino == ['K', '*']), 'DNAbin ambiguous standard translation', failures)
      call translate_dna([dna_t, dna_g, dna_a], amino, info, genetic_code=2)
      call check(info == 0 .and. amino(1) == 'W', 'DNAbin genetic code 2', failures)
      call translate_dna([dna_c, dna_t, dna_n], amino, info, genetic_code=3)
      call check(info == 0 .and. amino(1) == 'L', 'DNAbin genetic code 3 upstream table', failures)
      call translate_dna([dna_a, dna_g, dna_r], amino, info, genetic_code=5)
      call check(info == 0 .and. amino(1) == 'S', 'DNAbin genetic code 5 ambiguity', failures)
      call translate_dna([dna_t, dna_a, dna_r], amino, info, genetic_code=6)
      call check(info == 0 .and. amino(1) == '*', 'DNAbin genetic code 6 upstream table', failures)

      ambiguity_alignment(:, 1) = [dna_a, dna_r, dna_r]
      ambiguity_alignment(:, 2) = [dna_a, dna_y, dna_a]
      ambiguity_alignment(:, 3) = [dna_a, dna_a, dna_unknown]
      call dna_segregating_sites(ambiguity_alignment, segregating, trailing_gaps_as_n=.false.)
      call check(.not. segregating(1), 'DNAbin compatible ambiguous site', failures)
      call check(segregating(2), 'DNAbin incompatible ambiguous site', failures)
      call check(segregating(3), 'DNAbin historical unknown-site behavior', failures)

      terminal_alignment(1, :) = [dna_gap, dna_gap, dna_a, dna_gap, dna_c, dna_gap, dna_gap]
      call dna_leading_trailing_gaps_to_n(terminal_alignment, converted_alignment)
      call check(all(converted_alignment(1, :) == [dna_n, dna_n, dna_a, dna_gap, dna_c, dna_n, dna_n]), &
         'DNAbin terminal gaps to N', failures)
      call dna_base_proportions(reshape([dna_a, dna_c, dna_g, dna_t, dna_r, dna_m, dna_w, dna_s, dna_k, dna_y, &
         dna_v, dna_h, dna_d, dna_b, dna_n, dna_gap, dna_unknown], [1, 17]), state_proportions)
      call check(all(abs(state_proportions - 1.0_dp / 17.0_dp) < 1.0e-12_dp), &
         'DNAbin 17-state proportions', failures)
      call dna_distance([dna_a, dna_gap, dna_gap, dna_t], [dna_a, dna_gap, dna_t, dna_t], &
         'INDEL', distance, comparable, info)
      call check_close(distance, 1.0_dp, 1.0e-12_dp, 'INDEL distance', failures)

      alignment2(1, :) = long_a
      alignment2(2, :) = long_b
      call dna_distance_matrix(alignment2, 'K80', matrix_distance, info, pairwise_deletion=.true.)
      call check(info == 0 .and. matrix_distance(1, 2) > 0.0_dp, 'DNA distance matrix', failures)
      call dna_distance_matrix(alignment2, 'BH87', matrix_distance, info, pairwise_deletion=.true.)
      call check(info == 0 .and. matrix_distance(1, 2) > 0.0_dp, 'BH87 distance matrix', failures)

      alignment3(1, :) = long_a
      alignment3(2, :) = long_a
      alignment3(3, :) = long_a
      alignment3(2, 1) = dna_g
      alignment3(2, 2) = dna_a
      alignment3(3, 1) = dna_g
      alignment3(3, 3) = dna_t
      call dna_distance_matrix(alignment3, 'GG95', matrix_distance, info)
      call check(info == 0 .and. all(matrix_distance >= 0.0_dp), 'GG95 distance matrix', failures)
      call dna_distance_matrix_with_variance(alignment3, 'GG95', matrix_distance, matrix_variance, info)
      call check(info == 0 .and. all(matrix_variance >= 0.0_dp), 'GG95 variance matrix', failures)
      alignment2(1, :) = long_a
      alignment2(2, :) = variance_b
      call dna_distance_matrix_with_variance(alignment2, 'K80', matrix_distance, matrix_variance, info, &
         pairwise_deletion=.true.)
      call check(info == 0, 'K80 variance matrix status', failures)
      call check_close(matrix_variance(1, 2), 0.015341996173469385_dp, 1.0e-14_dp, &
         'K80 variance matrix', failures)

      call dna_distance_matrix_with_variance(alignment2, 'LOGDET', matrix_distance, matrix_variance, info)
      call check(info == 0, 'LOGDET variance matrix status', failures)
      call check_close(matrix_variance(1, 2), 0.011409023668639051_dp, 1.0e-13_dp, &
         'LOGDET variance matrix', failures)
      call dna_distance_matrix_with_variance(alignment2, 'PARALIN', matrix_distance, matrix_variance, info)
      call check(info == 0, 'PARALIN variance matrix status', failures)
      call check_close(matrix_variance(1, 2), -0.14013622222761954_dp, 1.0e-13_dp, &
         'PARALIN variance matrix', failures)

      alignment2 = dna_a
      alignment2(1, 1) = dna_gap
      call dna_distance_matrix(alignment2, 'INDEL', matrix_distance, info)
      call check(info == 0, 'INDEL matrix status', failures)
      call check_close(matrix_distance(1, 2), 1.0_dp, 1.0e-12_dp, 'INDEL matrix gap preservation', failures)
      allocate(table(4, 4))
      call dna_contingency_table(long_a, long_b, table)
      call check_close(sum(table), 20.0_dp, 1.0e-12_dp, 'DNA contingency table', failures)
   end subroutine test_dna_models

   subroutine test_graph_statistics(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed graph/statistics check.
      real(dp) :: distance(5, 5)
      real(dp) :: weight(5, 5)
      real(dp) :: x(5)
      integer, allocatable :: adjacency(:, :)
      integer, allocatable :: bins(:)
      real(dp), allocatable :: delta_bar(:)
      type(moran_result) :: moran
      real(dp) :: total_weight
      integer :: info
      integer :: i

      distance = reshape([ &
         0.0_dp, 1.0_dp, 4.0_dp, 7.0_dp, 8.0_dp, &
         1.0_dp, 0.0_dp, 2.0_dp, 6.0_dp, 7.0_dp, &
         4.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 5.0_dp, &
         7.0_dp, 6.0_dp, 3.0_dp, 0.0_dp, 4.0_dp, &
         8.0_dp, 7.0_dp, 5.0_dp, 4.0_dp, 0.0_dp], [5, 5])
      call minimum_spanning_tree(distance, adjacency, total_weight, info)
      call check(info == 0, 'MST status', failures)
      call check(sum(adjacency) == 8, 'MST edge count', failures)
      call check_close(total_weight, 10.0_dp, 1.0e-12_dp, 'MST total weight', failures)

      weight = 0.0_dp
      do i = 1, 4
         weight(i, i + 1) = 1.0_dp
         weight(i + 1, i) = 1.0_dp
      end do
      x = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp]
      call moran_i(x, weight, moran, info)
      call check(info == 0, 'Moran I status', failures)
      call check_close(moran%observed, 0.43346774193548376_dp, 1.0e-12_dp, 'Moran I observed', failures)
      call check_close(moran%expected, -0.25_dp, 1.0e-12_dp, 'Moran I expected', failures)
      call check_close(moran%standard_deviation, 0.34685285012751044_dp, 1.0e-12_dp, &
         'Moran I standard deviation', failures)

      call delta_plot_statistics(distance(1:4, 1:4), 10, bins, delta_bar, info)
      call check(info == 0, 'delta plot status', failures)
      call check(sum(bins) == 1, 'delta plot quartet count', failures)
      call check(all(delta_bar >= 0.0_dp .and. delta_bar <= 1.0_dp), 'delta plot range', failures)
   end subroutine test_graph_statistics

   subroutine test_misc_statistics(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed miscellaneous-statistics check.
      type(tree_count_result) :: count_result
      type(diversification_gof_result) :: gof
      type(diversification_time_result) :: time_fit
      type(chi_square_result) :: chi
      integer :: genes(3, 5)
      real(dp), allocatable :: distance(:, :)
      real(dp), allocatable :: variance(:, :)
      real(dp), allocatable :: individual_p(:)
      real(dp), allocatable :: contrasts(:)
      logical :: split_a(4)
      logical :: split_b(4)
      logical :: split_c(4)
      logical :: split_matrix(2, 4)
      integer :: info

      call howmanytrees(4, count_result, info)
      call check(info == 0, 'howmanytrees rooted binary status', failures)
      call check_close(count_result%value, 15.0_dp, 1.0e-12_dp, 'howmanytrees rooted binary count', failures)
      call howmanytrees(5, count_result, info, rooted=.false.)
      call check_close(count_result%value, 15.0_dp, 1.0e-12_dp, 'howmanytrees unrooted binary count', failures)
      call howmanytrees(4, count_result, info, binary=.false.)
      call check_close(count_result%value, 26.0_dp, 1.0e-12_dp, 'howmanytrees rooted multifurcating count', failures)
      call howmanytrees(4, count_result, info, labeled=.false.)
      call check_close(count_result%value, 2.0_dp, 1.0e-12_dp, 'howmanytrees unlabeled rooted count', failures)

      genes(1, :) = [1, 1, 2, 2, -9]
      genes(2, :) = [1, 2, 2, 2, 3]
      genes(3, :) = [2, 2, 2, 1, 3]
      call gene_distance_matrix(genes, distance, info, variance=variance, missing_value=-9)
      call check(info == 0, 'dist.gene count status', failures)
      call check_close(distance(1, 2), 1.0_dp, 1.0e-12_dp, 'dist.gene count 1-2', failures)
      call check_close(distance(1, 3), 3.0_dp, 1.0e-12_dp, 'dist.gene count 1-3', failures)
      call check_close(variance(1, 3), 0.75_dp, 1.0e-12_dp, 'dist.gene count variance', failures)
      call gene_distance_matrix(genes, distance, info, method='percentage', pairwise_deletion=.true., &
         variance=variance, missing_value=-9)
      call check_close(distance(1, 3), 0.75_dp, 1.0e-12_dp, 'dist.gene percentage 1-3', failures)
      call check_close(distance(2, 3), 0.4_dp, 1.0e-12_dp, 'dist.gene pairwise deletion', failures)
      call check_close(variance(1, 3), 0.046875_dp, 1.0e-12_dp, 'dist.gene percentage variance', failures)

      split_a = [.true., .true., .false., .false.]
      split_b = [.true., .false., .true., .false.]
      split_c = [.true., .true., .true., .false.]
      call check(.not. splits_compatible(split_a, split_b), 'incompatible split pair', failures)
      call check(splits_compatible(split_a, split_c), 'compatible split pair', failures)
      split_matrix(1, :) = split_a
      split_matrix(2, :) = split_c
      call check(all_splits_compatible(split_matrix), 'all compatible splits', failures)

      call diversification_gof([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], gof, info)
      call check(info == 0, 'diversi.gof status', failures)
      call check_close(gof%cramer_von_mises, 0.10074998262842674_dp, 1.0e-12_dp, &
         'diversi.gof Cramer-von Mises', failures)
      call check_close(gof%anderson_darling, 0.6526299309545924_dp, 1.0e-12_dp, &
         'diversi.gof Anderson-Darling', failures)

      call diversification_time([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], &
         time_fit, info, census=[1, 1, 1, 0, 1, 0], breakpoint=3.5_dp)
      call check(info == 0, 'diversi.time status', failures)
      call check(time_fit%n_event == 4 .and. time_fit%n_censored == 2, &
         'diversi.time event counts', failures)
      call check_close(time_fit%constant_rate, 0.19047619047619047_dp, 1.0e-14_dp, &
         'diversi.time constant rate', failures)
      call check_close(time_fit%log_likelihood_a, -10.63291230641413_dp, 1.0e-12_dp, &
         'diversi.time model A likelihood', failures)
      call check_close(time_fit%weibull_shape, 1.5820959702741118_dp, 1.0e-11_dp, &
         'diversi.time Weibull shape', failures)
      call check_close(time_fit%weibull_scale, 0.20674168905980708_dp, 1.0e-12_dp, &
         'diversi.time Weibull scale', failures)
      call check_close(time_fit%log_likelihood_b, -10.160512081140215_dp, 1.0e-11_dp, &
         'diversi.time model B likelihood', failures)
      call check_close(time_fit%early_rate, 0.18181818181818182_dp, 1.0e-14_dp, &
         'diversi.time early rate', failures)
      call check_close(time_fit%late_rate, 0.2222222222222222_dp, 1.0e-14_dp, &
         'diversi.time late rate', failures)
      call check_close(time_fit%log_likelihood_c, -10.61832167349155_dp, 1.0e-12_dp, &
         'diversi.time model C likelihood', failures)
      call check_close(time_fit%model_a_vs_b%p_value, 0.33104633874170913_dp, 1.0e-10_dp, &
         'diversi.time model A/B p-value', failures)
      call check_close(time_fit%model_a_vs_c%p_value, 0.8643611747061055_dp, 1.0e-10_dp, &
         'diversi.time model A/C p-value', failures)

      call slowinski_guyer_test(reshape([2.0_dp, 3.0_dp], [1, 2]), chi, individual_p, info)
      call check(info == 0 .and. chi%degrees_of_freedom == 2, 'Slowinski-Guyer status', failures)
      call check_close(individual_p(1), 0.75_dp, 1.0e-12_dp, 'Slowinski-Guyer individual p', failures)
      call check_close(chi%p_value, 0.75_dp, 1.0e-12_dp, 'Slowinski-Guyer combined p', failures)
      call mcconway_sims_test(reshape([4.0_dp, 7.0_dp, 8.0_dp, 3.0_dp], [2, 2]), chi, info)
      call check(info == 0 .and. chi%p_value >= 0.0_dp .and. chi%p_value <= 1.0_dp, &
         'McConway-Sims status', failures)

      call diversity_contrasts(reshape([2.0_dp, 8.0_dp, 4.0_dp, 8.0_dp, 2.0_dp, 4.0_dp], [3, 2]), &
         'ratiolog', contrasts, info)
      call check(info == 0, 'diversity contrasts status', failures)
      call check_close(contrasts(1), -3.0_dp, 1.0e-12_dp, 'diversity contrast negative', failures)
      call check_close(contrasts(2), 3.0_dp, 1.0e-12_dp, 'diversity contrast second row', failures)
      call check_close(contrasts(3), 0.0_dp, 1.0e-12_dp, 'diversity contrast tie', failures)
   end subroutine test_misc_statistics

   subroutine test_ordination(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed PCoA check.
      type(pcoa_result) :: result
      real(dp) :: distance(4, 4)
      integer :: info

      distance = reshape([ &
         0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
         1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, &
         1.0_dp, 1.0_dp, 0.0_dp, 3.0_dp, &
         1.0_dp, 3.0_dp, 3.0_dp, 0.0_dp], [4, 4])
      call pcoa(distance, result, info)
      call check(info == 0, 'PCoA uncorrected status', failures)
      call check(result%positive_rank == 2, 'PCoA positive rank', failures)
      call check_close(result%trace, 5.5_dp, 1.0e-12_dp, 'PCoA trace', failures)
      call check_close(result%eigenvalues(1), 5.964101615137754_dp, 1.0e-12_dp, &
         'PCoA leading eigenvalue', failures)
      call check_close(result%eigenvalues(2), 0.5_dp, 1.0e-12_dp, 'PCoA second eigenvalue', failures)
      call check_close(result%eigenvalues(4), -0.964101615137754_dp, 1.0e-12_dp, &
         'PCoA negative eigenvalue', failures)

      call pcoa(distance, result, info, 'lingoes')
      call check(info == 0 .and. result%correction == pcoa_lingoes, 'PCoA Lingoes status', failures)
      call check_close(result%correction_constant, 0.964101615137754_dp, 1.0e-12_dp, &
         'PCoA Lingoes constant', failures)
      call check_close(result%corrected_eigenvalues(1), 6.928203230275509_dp, 1.0e-12_dp, &
         'PCoA Lingoes leading eigenvalue', failures)
      call check_close(result%corrected_eigenvalues(2), 1.464101615137754_dp, 1.0e-12_dp, &
         'PCoA Lingoes second eigenvalue', failures)

      call pcoa(distance, result, info, 'cailliez')
      call check(info == 0 .and. result%correction == pcoa_cailliez, 'PCoA Cailliez status', failures)
      call check_close(result%correction_constant, 1.146264369941973_dp, 1.0e-11_dp, &
         'PCoA Cailliez constant', failures)
      call check_close(result%corrected_eigenvalues(1), 10.89897948556636_dp, 1.0e-10_dp, &
         'PCoA Cailliez leading eigenvalue', failures)
      call check_close(result%corrected_eigenvalues(2), 2.303225372841205_dp, 1.0e-10_dp, &
         'PCoA Cailliez second eigenvalue', failures)
   end subroutine test_ordination

   subroutine test_continuous_ace(failures)
      integer, intent(inout) :: failures !! Running count incremented for every continuous-ACE check.
      type(phylo_tree) :: tree
      type(ace_continuous_result) :: fit
      integer :: edge(6, 2)
      real(dp) :: edge_length(6)
      integer :: info

      edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
      edge_length = 1.0_dp
      tree = make_phylo_tree(4, edge, edge_length)
      call ace_continuous_ml(tree, [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp], fit, info)
      call check(info == 0, 'ace continuous ML status', failures)
      call check_close(fit%estimates(1), 4.0_dp, 1.0e-12_dp, 'ace ML root estimate', failures)
      call check_close(fit%estimates(2), 8.0_dp / 3.0_dp, 1.0e-12_dp, 'ace ML left estimate', failures)
      call check_close(fit%estimates(3), 16.0_dp / 3.0_dp, 1.0e-12_dp, 'ace ML right estimate', failures)
      call check_close(fit%sigma2, 14.0_dp / 9.0_dp, 1.0e-12_dp, 'ace ML sigma2', failures)
      call check_close(fit%sigma2_standard_error, (14.0_dp / 9.0_dp) / sqrt(6.0_dp), 1.0e-12_dp, &
         'ace ML sigma2 SE', failures)
      call check_close(fit%standard_error(1), 0.763762615825973_dp, 1.0e-12_dp, &
         'ace ML root SE', failures)
      call check_close(fit%log_likelihood, -4.325498256837117_dp, 1.0e-12_dp, &
         'ace ML log likelihood', failures)

      call ace_continuous_reml(tree, [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp], fit, info)
      call check(info == 0, 'ace continuous REML status', failures)
      call check_close(fit%estimates(1), 4.0_dp, 1.0e-12_dp, 'ace REML root estimate', failures)
      call check_close(fit%estimates(2), 8.0_dp / 3.0_dp, 1.0e-12_dp, 'ace REML left estimate', failures)
      call check_close(fit%estimates(3), 16.0_dp / 3.0_dp, 1.0e-12_dp, 'ace REML right estimate', failures)
      call check_close(fit%sigma2, 7.0_dp / 3.0_dp, 1.0e-12_dp, 'ace REML sigma2', failures)
      call check_close(fit%sigma2_standard_error, (7.0_dp / 3.0_dp) * sqrt(0.5_dp), 1.0e-12_dp, &
         'ace REML sigma2 SE', failures)
      call check_close(fit%standard_error(1), 1.322875655532295_dp, 1.0e-12_dp, &
         'ace REML root SE', failures)
      call check_close(fit%residual_log_likelihood, -4.541893581161611_dp, 1.0e-12_dp, &
         'ace REML residual likelihood', failures)

      call ace_continuous_gls(tree, [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp], 'brownian', result=fit, info=info)
      call check(info == 0, 'continuous ACE GLS Brownian status', failures)
      call check_close(fit%estimates(1), 4.0_dp, 2.0e-12_dp, 'continuous ACE GLS root', failures)
      call check_close(fit%estimates(2), 8.0_dp / 3.0_dp, 2.0e-12_dp, 'continuous ACE GLS left', failures)
      call check_close(fit%estimates(3), 16.0_dp / 3.0_dp, 2.0e-12_dp, 'continuous ACE GLS right', failures)
      call check_close(fit%standard_error(1), 0.0_dp, 2.0e-12_dp, 'continuous ACE GLS root SE', failures)
      call check_close(fit%standard_error(2), 1.0_dp / sqrt(3.0_dp), 2.0e-12_dp, &
         'continuous ACE GLS left SE', failures)
      call check_close(fit%sigma2, 14.0_dp / 3.0_dp, 2.0e-12_dp, 'continuous ACE GLS sigma2', failures)
   end subroutine test_continuous_ace


   subroutine test_discrete_ace(failures)
      integer, intent(inout) :: failures !! Running count incremented for every discrete ancestral-state likelihood check.
      type(phylo_tree) :: tree
      type(ace_discrete_result) :: fit
      integer, allocatable :: index_matrix(:, :)
      real(dp), allocatable :: ancestral(:, :)
      integer :: edge(6, 2)
      integer :: info
      integer :: states(4)
      real(dp) :: deviance
      real(dp) :: rates(1)

      edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
      tree = make_phylo_tree(4, edge, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
      states = [1, 1, 2, 2]
      call ace_rate_index_matrix(2, 'ER', index_matrix, info)
      call check(info == 0, 'discrete ACE ER index status', failures)
      call check(all(index_matrix == reshape([0, 1, 1, 0], [2, 2])), 'discrete ACE ER index', failures)

      rates = 0.2746530688361161_dp
      call ace_discrete_likelihood(tree, states, index_matrix, rates, deviance, ancestral, info)
      call check(info == 0, 'discrete ACE fixed likelihood status', failures)
      call check_close(deviance, 3.5835189384561095_dp, 2.0e-11_dp, 'discrete ACE fixed deviance', failures)
      call check_close(ancestral(1, 1), 0.5_dp, 2.0e-12_dp, 'discrete ACE root likelihood', failures)
      call check_close(ancestral(2, 1), 0.884900180742_dp, 2.0e-11_dp, 'discrete ACE smoothed child', failures)
      call check_close(ancestral(3, 1), 0.115099819258_dp, 2.0e-11_dp, 'discrete ACE smoothed sibling', failures)

      call ace_discrete_fit(tree, states, 2, 'ER', fit, info, tolerance=1.0e-10_dp, max_iter=1000)
      call check(info == 0, 'discrete ACE fit status', failures)
      call check(fit%converged, 'discrete ACE optimizer convergence', failures)
      call check_close(fit%rates(1), 0.2746530688361161_dp, 1.0e-7_dp, 'discrete ACE fitted rate', failures)
      call check_close(fit%log_likelihood, -1.7917594692280547_dp, 2.0e-11_dp, 'discrete ACE log likelihood', failures)
      call check_close(fit%standard_error(1), 0.25_dp, 2.0e-6_dp, 'discrete ACE rate SE', failures)

      call ace_rate_index_matrix(3, 'ARD', index_matrix, info)
      call check(info == 0, 'discrete ACE ARD index status', failures)
      call check(all(index_matrix == reshape([0, 1, 2, 3, 0, 4, 5, 6, 0], [3, 3])), &
         'discrete ACE ARD column-major index', failures)
      call ace_rate_index_matrix(3, 'SYM', index_matrix, info)
      call check(info == 0, 'discrete ACE SYM index status', failures)
      call check(all(index_matrix == reshape([0, 1, 2, 1, 0, 3, 2, 3, 0], [3, 3])), &
         'discrete ACE SYM index', failures)
   end subroutine test_discrete_ace


   subroutine test_pgls(failures)
      integer, intent(inout) :: failures !! Running count incremented for every phylogenetic GLS/correlation check.
      type(phylo_tree) :: tree
      type(pgls_result) :: fit
      real(dp), allocatable :: matrix(:, :)
      real(dp) :: design(4, 2)
      real(dp) :: response(4)
      integer :: edge(6, 2)
      integer :: info

      edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
      tree = make_phylo_tree(4, edge, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
      call cor_brownian(tree, matrix, info)
      call check(info == 0, 'PGLS Brownian status', failures)
      call check_close(matrix(1, 2), 0.5_dp, 2.0e-12_dp, 'PGLS Brownian sister correlation', failures)
      call cor_martins(tree, 0.5_dp, matrix, info)
      call check(info == 0, 'PGLS Martins status', failures)
      call check_close(matrix(1, 2), exp(-1.0_dp), 2.0e-12_dp, 'PGLS Martins sister correlation', failures)
      call check_close(matrix(1, 3), exp(-2.0_dp), 2.0e-12_dp, 'PGLS Martins cross correlation', failures)
      call cor_grafen(tree, 2.0_dp, matrix, info)
      call check(info == 0, 'PGLS Grafen status', failures)
      call check_close(matrix(1, 2), 8.0_dp / 9.0_dp, 2.0e-12_dp, 'PGLS Grafen sister correlation', failures)
      call cor_pagel(tree, 0.3_dp, matrix, info)
      call check(info == 0, 'PGLS Pagel status', failures)
      call check_close(matrix(1, 2), 0.15_dp, 2.0e-12_dp, 'PGLS Pagel sister correlation', failures)
      call cor_blomberg(tree, 2.0_dp, matrix, info)
      call check(info == 0, 'PGLS Blomberg status', failures)
      call check_close(matrix(1, 2), 1.0_dp / sqrt(2.0_dp), 2.0e-12_dp, 'PGLS Blomberg sister correlation', failures)

      design(:, 1) = 1.0_dp
      design(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
      response = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
      call cor_brownian(tree, matrix, info)
      call pgls_fit(response, design, matrix, fit, info)
      call check(info == 0, 'PGLS fixed covariance fit status', failures)
      call check_close(fit%coefficients(1), 3.0_dp / 14.0_dp, 2.0e-12_dp, 'PGLS intercept', failures)
      call check_close(fit%coefficients(2), 33.0_dp / 14.0_dp, 2.0e-12_dp, 'PGLS slope', failures)
      call check_close(fit%sigma2, 8.0_dp / 7.0_dp, 2.0e-12_dp, 'PGLS residual variance', failures)
      call check_close(fit%log_likelihood, -5.6551348456159545_dp, 2.0e-12_dp, 'PGLS ML log likelihood', failures)

      call pgls_fit_model(tree, response, design, 'pagel', fit, info, max_iter=500, tolerance=1.0e-10_dp)
      call check(info == 0, 'PGLS Pagel profile status', failures)
      call check(fit%converged, 'PGLS Pagel profile convergence', failures)
      call check_close(fit%correlation_parameter, 0.0_dp, 2.0e-8_dp, 'PGLS Pagel profiled lambda', failures)
      call check_close(fit%coefficients(1), 0.3_dp, 2.0e-10_dp, 'PGLS Pagel profiled intercept', failures)
      call check_close(fit%coefficients(2), 2.3_dp, 2.0e-10_dp, 'PGLS Pagel profiled slope', failures)
      call check_close(fit%log_likelihood, -4.568983656449124_dp, 2.0e-10_dp, &
         'PGLS Pagel profiled log likelihood', failures)
   end subroutine test_pgls


   subroutine test_chronopl(failures)
      integer, intent(inout) :: failures !! Running count incremented for every penalized-likelihood dating check.
      type(phylo_tree) :: tree
      type(chronopl_result) :: fit
      real(dp) :: all_age(7)
      real(dp) :: objective
      integer :: edge(6, 2)
      integer :: info

      edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
      tree = make_phylo_tree(4, edge, [0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp])
      all_age = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.5_dp, 0.5_dp]
      objective = chronopl_objective(tree, tree%edge_length, &
         [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp], all_age, 1.0_dp, info)
      call check(info == 0, 'chronopl fixed objective status', failures)
      call check_close(objective, 2.9898105231489573_dp, 2.0e-12_dp, 'chronopl fixed objective', failures)

      call chronopl_fit(tree, 1.0_dp, fit, info, tolerance=1.0e-10_dp, max_iter=1000)
      call check(info == 0, 'chronopl fit status', failures)
      call check(fit%converged, 'chronopl fit convergence', failures)
      call check_close(fit%node_age(1), 1.0_dp, 2.0e-12_dp, 'chronopl calibrated root age', failures)
      call check_close(fit%node_age(2), 0.5_dp, 4.0e-6_dp, 'chronopl left internal age', failures)
      call check_close(fit%node_age(3), 0.5_dp, 4.0e-6_dp, 'chronopl right internal age', failures)
      call check(maxval(abs(fit%rates - 0.5_dp)) < 5.0e-6_dp, 'chronopl equal fitted rates', failures)
      call check(maxval(abs(fit%tree%edge_length - 0.5_dp)) < 4.0e-6_dp, 'chronopl dated branch lengths', failures)
      call check_close(fit%penalized_log_likelihood, -2.9898105231489573_dp, 1.0e-9_dp, &
         'chronopl penalized log likelihood', failures)
   end subroutine test_chronopl



   subroutine test_chronos_clock(failures)
      integer, intent(inout) :: failures !! Running count incremented for every chronos clock-model check.
      type(phylo_tree) :: tree
      type(chronos_clock_result) :: fit
      integer :: edge(6, 2)
      integer :: info

      edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
      tree = make_phylo_tree(4, edge, [0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp])
      call chronos_clock_fit(tree, fit, info, tolerance=1.0e-10_dp, max_iter=1000)
      call check(info == 0, 'chronos clock status', failures)
      call check(fit%converged, 'chronos clock convergence', failures)
      call check_close(fit%rate, 0.5_dp, 2.0e-12_dp, 'chronos clock rate', failures)
      call check_close(fit%node_age(1), 1.0_dp, 2.0e-12_dp, 'chronos clock root age', failures)
      call check_close(fit%node_age(2), 0.5_dp, 2.0e-12_dp, 'chronos clock left age', failures)
      call check_close(fit%node_age(3), 0.5_dp, 2.0e-12_dp, 'chronos clock right age', failures)
      call check_close(fit%log_likelihood, -2.9898105231489573_dp, 2.0e-12_dp, &
         'chronos clock log likelihood', failures)
      call check(fit%n_parameters == 3, 'chronos clock parameter count', failures)
      call check_close(fit%phiic, 11.979621046297915_dp, 2.0e-12_dp, 'chronos clock PHIIC', failures)
   end subroutine test_chronos_clock

   subroutine test_chronos_models(failures)
      integer, intent(inout) :: failures !! Running count incremented for every failed multi-model chronos check.
      type(phylo_tree) :: tree
      type(chronos_result) :: fit
      integer :: edge(6, 2)
      real(dp) :: ages(7)
      real(dp) :: edge_length(6)
      real(dp) :: frequencies(2)
      real(dp) :: rates(6)
      real(dp) :: value
      integer :: info

      edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
      edge_length = 0.25_dp
      tree = make_phylo_tree(4, edge, edge_length)
      ages = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.5_dp, 0.5_dp]
      rates = [0.4_dp, 0.6_dp, 0.5_dp, 0.7_dp, 0.3_dp, 0.8_dp]

      value = chronos_objective(tree, edge_length, rates, ages, 1.2_dp, 'correlated', info)
      call check(info == 0, 'chronos correlated objective status', failures)
      call check_close(value, 3.328103461253784_dp, 2.0e-11_dp, 'chronos correlated objective', failures)

      value = chronos_objective(tree, edge_length, rates, ages, 1.2_dp, 'relaxed', info)
      call check(info == 0, 'chronos relaxed objective status', failures)
      call check_close(value, 3.413580249477147_dp, 2.0e-10_dp, 'chronos relaxed objective', failures)

      frequencies = [0.25_dp, 0.75_dp]
      value = chronos_objective(tree, edge_length, [0.4_dp, 0.8_dp], ages, 1.0_dp, 'discrete', info, frequencies)
      call check(info == 0, 'chronos discrete objective status', failures)
      call check_close(value, 3.144332090588235_dp, 2.0e-11_dp, 'chronos discrete objective', failures)

      call chronos_fit(tree, fit, info, lambda=1.0_dp, model='correlated', age_min=[1.0_dp], age_max=[1.0_dp], &
         calibration_nodes=[5], max_iter=1000, tolerance=1.0e-10_dp)
      call check(info == 0, 'chronos correlated fit status', failures)
      call check(fit%converged, 'chronos correlated convergence', failures)
      call check_close(fit%log_likelihood, -2.989810523148957_dp, 2.0e-10_dp, 'chronos correlated log likelihood', failures)
      call check_close(maxval(abs(fit%rates - 0.5_dp)), 0.0_dp, 2.0e-7_dp, 'chronos correlated equal rates', failures)
      call check_close(fit%phiic, 21.97962104629791_dp, 3.0e-9_dp, 'chronos correlated PHIIC', failures)
   end subroutine test_chronos_models

   subroutine test_compar_ou(failures)
      integer, intent(inout) :: failures !! Running count incremented for every OU comparative-model check.
      type(phylo_tree) :: tree
      type(compar_ou_result) :: fit
      real(dp) :: phenotype(4)
      real(dp) :: deviance
      integer :: edge(6, 2)
      integer :: info

      edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
      tree = make_phylo_tree(4, edge, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
      phenotype = [-0.989121350348_dp, -0.367786651468_dp, 1.287925261289_dp, 0.193974419133_dp]
      call compar_ou_fit(tree, phenotype, fit, info, max_iter=1000, tolerance=1.0e-10_dp)
      call check(info == 0, 'compar.ou fit status', failures)
      call check(fit%converged, 'compar.ou optimizer convergence', failures)
      call check_close(fit%alpha, 0.128810945177_dp, 8.0e-7_dp, 'compar.ou alpha', failures)
      call check_close(fit%sigma2, 1.74229135723_dp, 8.0e-6_dp, 'compar.ou sigma2', failures)
      call check_close(fit%theta(1), 0.137587793223_dp, 8.0e-7_dp, 'compar.ou theta', failures)
      call check_close(fit%log_likelihood, -4.7559726111697_dp, 2.0e-11_dp, 'compar.ou log likelihood', failures)
      deviance = compar_ou_likelihood(tree, phenotype, fit%alpha, fit%sigma2, fit%theta, info)
      call check(info == 0, 'compar.ou fixed likelihood status', failures)
      call check_close(deviance, 9.5119452223394_dp, 3.0e-11_dp, 'compar.ou fixed deviance', failures)
   end subroutine test_compar_ou

   subroutine test_compar_lynch(failures)
      integer, intent(inout) :: failures !! Running failure count incremented by Lynch comparative-method checks.
      type(compar_lynch_result) :: fit
      real(dp) :: g(4, 4)
      real(dp) :: x(4, 2)
      integer :: info

      x(:, 1) = [1.0_dp, 2.0_dp, 4.0_dp, 6.0_dp]
      x(:, 2) = [2.0_dp, 1.0_dp, 5.0_dp, 4.0_dp]
      g = reshape([ &
         2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, &
         0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp], [4, 4])

      call compar_lynch_fit(x, g, fit, info)
      call check(info == 0, 'compar.lynch status', failures)
      call check(fit%converged, 'compar.lynch convergence', failures)
      call check(fit%iterations == 174, 'compar.lynch iterations', failures)
      call check_close(fit%log_likelihood, -13.1409137297564_dp, 2.0e-11_dp, &
         'compar.lynch log likelihood', failures)
      call check_close(fit%mean(1), 3.25_dp, 1.0e-12_dp, 'compar.lynch mean trait 1', failures)
      call check_close(fit%mean(2), 3.0_dp, 1.0e-12_dp, 'compar.lynch mean trait 2', failures)
      call check_close(fit%phylogenetic_covariance(1, 1), 1.035734361732426_dp, 2.0e-12_dp, &
         'compar.lynch phylogenetic variance', failures)
      call check_close(fit%phylogenetic_covariance(1, 2), 0.882280788893202_dp, 2.0e-12_dp, &
         'compar.lynch phylogenetic covariance', failures)
      call check_close(fit%environmental_covariance(1, 1), 0.632347155545047_dp, 2.0e-12_dp, &
         'compar.lynch environmental variance', failures)
      call check_close(fit%environmental_covariance(1, 2), -0.368636598680333_dp, 2.0e-12_dp, &
         'compar.lynch environmental covariance', failures)
   end subroutine test_compar_lynch

   subroutine test_corphylo(failures)
      integer, intent(inout) :: failures !! Running failure count incremented by multivariate corphylo checks.
      type(phylo_tree) :: tree
      type(corphylo_result) :: fit
      integer :: edge(6, 2)
      real(dp) :: x(4, 2)
      real(dp) :: objective
      integer :: info

      edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
      tree = make_phylo_tree(4, edge, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
      x(:, 1) = [1.0_dp, 2.0_dp, 4.0_dp, 6.0_dp]
      x(:, 2) = [2.0_dp, 1.0_dp, 5.0_dp, 4.0_dp]

      objective = corphylo_objective(x, tree, [1.0_dp, 0.2_dp, 0.8_dp], [0.5_dp, 0.7_dp], info)
      call check(info == 0, 'corphylo fixed objective status', failures)
      call check_close(objective, 3.101923885322540_dp, 2.0e-12_dp, 'corphylo fixed objective', failures)

      call corphylo_fit(x, tree, fit, info, constrain_d=.true.)
      call check(info == 0, 'corphylo fit status', failures)
      call check(fit%converged, 'corphylo convergence', failures)
      call check_close(fit%log_likelihood, -7.271053234382986_dp, 3.0e-8_dp, 'corphylo log likelihood', failures)
      call check_close(fit%d(1), 0.9999546021312976_dp, 4.0e-7_dp, 'corphylo d trait 1', failures)
      call check_close(fit%d(2), 0.9999546021312976_dp, 4.0e-7_dp, 'corphylo d trait 2', failures)
      call check_close(fit%correlation(1, 2), 0.3897672580970114_dp, 2.0e-7_dp, &
         'corphylo trait correlation', failures)
   end subroutine test_corphylo

   subroutine test_binary_pglmm(failures)
      integer, intent(inout) :: failures !! Running failure count incremented by binary phylogenetic GLMM checks.
      type(phylo_tree) :: tree
      type(binary_pglmm_result) :: fit
      integer :: edge(14, 2)
      integer :: y(8)
      real(dp) :: design(8, 2)
      real(dp) :: fixed_design(4, 2)
      real(dp) :: h(4)
      real(dp) :: inv_weight(4)
      real(dp) :: objective
      real(dp) :: vphy(4, 4)
      integer :: info

      fixed_design(:, 1) = 1.0_dp
      fixed_design(:, 2) = [-1.0_dp, -0.2_dp, 0.4_dp, 1.2_dp]
      h = [0.2_dp, -0.1_dp, 0.4_dp, -0.3_dp]
      inv_weight = [4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp]
      vphy = reshape([ &
         1.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
         0.5_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 1.0_dp, 0.5_dp, &
         0.0_dp, 0.0_dp, 0.5_dp, 1.0_dp], [4, 4])
      objective = binary_pglmm_reml_objective(0.7_dp, inv_weight, h, vphy, fixed_design, info)
      call check(info == 0, 'binaryPGLMM fixed REML status', failures)
      call check_close(objective, 5.94086367491974_dp, 2.0e-12_dp, 'binaryPGLMM fixed REML objective', failures)

      edge = reshape([ &
         9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, &
         10, 11, 12, 13, 14, 15, 1, 2, 3, 4, 5, 6, 7, 8], [14, 2])
      tree = make_phylo_tree(8, edge, [ &
         1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
         1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
      design(:, 1) = 1.0_dp
      design(:, 2) = [-1.5_dp, -1.0_dp, -0.5_dp, 0.0_dp, 0.2_dp, 0.7_dp, 1.1_dp, 1.5_dp]
      y = [0, 1, 0, 0, 1, 1, 0, 1]
      call binary_pglmm_fit(y, design, tree, fit, info)
      call check(info == 0, 'binaryPGLMM fit status', failures)
      call check(fit%converged, 'binaryPGLMM convergence', failures)
      call check(fit%iterations == 13, 'binaryPGLMM PQL iterations', failures)
      call check_close(fit%s2, 3.81156655_dp, 3.0e-6_dp, 'binaryPGLMM phylogenetic variance', failures)
      call check_close(fit%coefficients(1), -0.035752787_dp, 3.0e-7_dp, 'binaryPGLMM intercept', failures)
      call check_close(fit%coefficients(2), 0.58172124_dp, 3.0e-7_dp, 'binaryPGLMM slope', failures)
      call check_close(fit%conditional_reml_log_likelihood, -14.25504169_dp, 2.0e-6_dp, &
         'binaryPGLMM conditional REML likelihood', failures)
      call check_close(fit%s2_p_value, 0.21011554_dp, 2.0e-6_dp, 'binaryPGLMM variance p-value', failures)
   end subroutine test_binary_pglmm

   subroutine test_reconstruct(failures)
      integer, intent(inout) :: failures !! Running failure count incremented by reconstruct ancestral-model checks.
      type(phylo_tree) :: tree
      type(reconstruct_result) :: fit
      integer :: edge(6, 2)
      real(dp) :: x(4)
      integer :: info

      edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
      tree = make_phylo_tree(4, edge, [1.0_dp, 0.8_dp, 0.7_dp, 1.2_dp, 1.0_dp, 1.4_dp])
      x = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp]

      call reconstruct_fit(x, tree, 'GLS', fit, info)
      call check(info == 0, 'reconstruct GLS status', failures)
      call check_close(fit%ancestral(1), 3.827693262962_dp, 2.0e-12_dp, 'reconstruct GLS root', failures)
      call check_close(fit%ancestral(2), 2.377832971127_dp, 2.0e-12_dp, 'reconstruct GLS node 6', failures)
      call check_close(fit%ancestral(3), 4.987581496430_dp, 2.0e-12_dp, 'reconstruct GLS node 7', failures)
      call check_close(fit%sigma2, 3.237089930663_dp, 2.0e-12_dp, 'reconstruct GLS variance', failures)

      call reconstruct_fit(x, tree, 'GLS_ABM', fit, info)
      call check(info == 0, 'reconstruct GLS ABM status', failures)
      call check_close(fit%trend, 4.880086647068_dp, 2.0e-12_dp, 'reconstruct ABM trend', failures)
      call check_close(fit%ancestral(1), -5.572799009748_dp, 2.0e-12_dp, 'reconstruct ABM root', failures)
      call check_close(fit%sigma2, 2.466346897726_dp, 2.0e-12_dp, 'reconstruct ABM variance', failures)

      call reconstruct_fit(x, tree, 'GLS_OUS', fit, info, alpha=0.4_dp)
      call check(info == 0, 'reconstruct GLS OUS status', failures)
      call check_close(fit%ancestral(1), 3.998509547246_dp, 2.0e-12_dp, 'reconstruct OUS root', failures)
      call check_close(fit%log_likelihood, -8.706025236551_dp, 2.0e-12_dp, 'reconstruct OUS likelihood', failures)

      call reconstruct_fit(x, tree, 'GLS_OU', fit, info, alpha=0.4_dp)
      call check(info == 0, 'reconstruct GLS OU status', failures)
      call check_close(fit%theta, 16.55123482837_dp, 2.0e-11_dp, 'reconstruct OU theta', failures)
      call check_close(fit%ancestral(1), -10.83649298203_dp, 2.0e-11_dp, 'reconstruct OU root', failures)
      call check_close(fit%log_likelihood, -7.996032862231_dp, 2.0e-12_dp, 'reconstruct OU likelihood', failures)

      call reconstruct_fit(x, tree, 'ML', fit, info)
      call check(info == 0, 'reconstruct ML status', failures)
      call check_close(fit%sigma2, 1.618544965332_dp, 2.0e-12_dp, 'reconstruct ML variance', failures)
      call check_close(fit%standard_error(1), 1.069008971018_dp, 2.0e-12_dp, 'reconstruct ML root SE', failures)

      call reconstruct_fit(x, tree, 'REML', fit, info)
      call check(info == 0, 'reconstruct REML status', failures)
      call check_close(fit%sigma2, 2.427817447998_dp, 2.0e-12_dp, 'reconstruct REML variance', failures)
      call check_close(fit%standard_error(1), 1.309263254726_dp, 2.0e-12_dp, 'reconstruct REML root SE', failures)
   end subroutine test_reconstruct

   subroutine test_skyline(failures)
      integer, intent(inout) :: failures !! Running count incremented for every deterministic skyline check.
      type(skyline_result) :: result
      integer, allocatable :: group(:)
      integer :: groups
      integer :: info

      call collapsed_intervals([0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp], 0.35_dp, group, groups, info)
      call check(info == 0 .and. groups == 2, 'collapsed intervals status', failures)
      call check(all(group == [1, 1, 1, 2]), 'collapsed intervals grouping', failures)

      call skyline_from_intervals([5, 4, 3, 2], [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp], 0.35_dp, result, info)
      call check(info == 0 .and. result%parameter_count == 2, 'skyline modern status', failures)
      call check_close(result%time(1), 0.6_dp, 1.0e-12_dp, 'skyline first time', failures)
      call check_close(result%time(2), 1.0_dp, 1.0e-12_dp, 'skyline total time', failures)
      call check_close(result%population_size(1), 31.0_dp / 30.0_dp, 1.0e-12_dp, &
         'skyline first population size', failures)
      call check_close(result%population_size(2), 0.4_dp, 1.0e-12_dp, &
         'skyline second population size', failures)
      call check_close(result%log_likelihood, 2.010878114295393_dp, 1.0e-12_dp, &
         'skyline log likelihood', failures)
      call check_close(result%log_likelihood_aicc, -5.989121885704607_dp, 1.0e-12_dp, &
         'skyline AICc likelihood', failures)

      call skyline_from_intervals([5, 4, 3, 2], [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp], &
         0.35_dp, result, info, old_style=.true.)
      call check(info == 0, 'skyline old-style status', failures)
      call check_close(result%population_size(1), 1.0_dp, 1.0e-12_dp, &
         'skyline old-style first size', failures)
      call check_close(result%log_likelihood, 2.009247582764366_dp, 1.0e-12_dp, &
         'skyline old-style likelihood', failures)
   end subroutine test_skyline

   subroutine test_birthdeath(failures)
      integer, intent(inout) :: failures !! Running count incremented for every standard birth-death check.
      type(birthdeath_result) :: fit
      real(dp) :: times(4)
      integer :: info

      times = [3.613551448368_dp, 3.246318433096_dp, 1.406539523638_dp, 0.085104551586_dp]
      call birthdeath_from_times(times, fit, info)
      call check(info == 0, 'birthdeath interior status', failures)
      call check_close(fit%death_birth_ratio, 0.54988562_dp, 2.0e-6_dp, &
         'birthdeath d/b estimate', failures)
      call check_close(fit%net_diversification, 0.13890576_dp, 2.0e-7_dp, &
         'birthdeath b-d estimate', failures)
      call check_close(fit%deviance, 7.90608755910693_dp, 1.0e-11_dp, &
         'birthdeath deviance', failures)
      call check_close(fit%death_birth_ratio_se, 1.65335_dp, 2.0e-5_dp, &
         'birthdeath d/b SE', failures)
      call check_close(fit%net_diversification_se, 0.419189_dp, 2.0e-6_dp, &
         'birthdeath b-d SE', failures)
      call check(fit%confidence_interval(1, 1) < fit%death_birth_ratio &
         .and. fit%confidence_interval(1, 2) > fit%death_birth_ratio, &
         'birthdeath ratio interval contains estimate', failures)

      call birthdeath_from_times([2.0_dp, 1.0_dp, 1.0_dp], fit, info)
      call check(info == 0, 'birthdeath boundary status', failures)
      call check_close(fit%death_birth_ratio, 0.0_dp, 1.0e-14_dp, &
         'birthdeath boundary d/b', failures)
      call check_close(fit%net_diversification, 1.0_dp / 3.0_dp, 3.0e-8_dp, &
         'birthdeath boundary b-d', failures)
      call check_close(fit%deviance, 4.81093021621633_dp, 1.0e-11_dp, &
         'birthdeath boundary deviance', failures)
   end subroutine test_birthdeath


   subroutine test_birthdeath_extended(failures)
      integer, intent(inout) :: failures !! Running count incremented for every extended birth-death check.
      type(birthdeath_extended_result) :: fit
      type(phylo_tree) :: tree
      real(dp) :: terminal_length(4)
      real(dp) :: times(3)
      integer :: edge(6, 2)
      integer :: species_count(4)
      integer :: info

      times = [4.85113079542_dp, 1.420611460025_dp, 0.845227443437_dp]
      terminal_length = [2.590021556028_dp, 0.597010750104_dp, 3.124978982577_dp, 3.887881909425_dp]
      species_count = [11, 18, 23, 27]
      call check_close(birthdeath_extended_deviance(times, terminal_length, species_count, 0.2_dp, 0.7_dp, .true.), &
         56.176564750257995_dp, 2.0e-12_dp, 'bd.ext conditional fixed deviance', failures)
      call check_close(birthdeath_extended_deviance(times, terminal_length, species_count, 0.2_dp, 0.7_dp, .false.), &
         57.54108422675587_dp, 2.0e-12_dp, 'bd.ext unconditional fixed deviance', failures)

      call birthdeath_extended_from_data(times, terminal_length, species_count, fit, info, .true.)
      call check(info == 0, 'bd.ext conditional status', failures)
      call check_close(fit%death_birth_ratio, 0.38669621_dp, 8.0e-7_dp, 'bd.ext conditional d/b', failures)
      call check_close(fit%net_diversification, 1.22567227_dp, 8.0e-7_dp, 'bd.ext conditional b-d', failures)
      call check_close(fit%deviance, 43.78698130255824_dp, 2.0e-11_dp, 'bd.ext conditional deviance', failures)

      call birthdeath_extended_from_data(times, terminal_length, species_count, fit, info, .false.)
      call check(info == 0, 'bd.ext unconditional status', failures)
      call check_close(fit%death_birth_ratio, 0.0_dp, 1.0e-14_dp, 'bd.ext unconditional boundary d/b', failures)
      call check_close(fit%net_diversification, 1.47064037_dp, 2.0e-6_dp, 'bd.ext unconditional b-d', failures)
      call check_close(fit%deviance, 44.20532536249572_dp, 2.0e-11_dp, 'bd.ext unconditional deviance', failures)

      edge = reshape([5, 5, 6, 6, 7, 7, 1, 6, 2, 7, 3, 4], [6, 2])
      tree = make_phylo_tree(4, edge, [3.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
      call birthdeath_extended_fit(tree, [5, 2, 3, 4], fit, info)
      call check(info == 0, 'bd.ext tree entry status', failures)
      call check_close(fit%death_birth_ratio, 0.0_dp, 1.0e-14_dp, 'bd.ext tree entry d/b', failures)
      call check_close(fit%net_diversification, 0.89242813_dp, 2.0e-6_dp, 'bd.ext tree entry b-d', failures)
   end subroutine test_birthdeath_extended

   function edge_length_to_child(tree, child, found) result(value)
      type(phylo_tree), intent(in) :: tree !! Tree whose incoming edge length is requested.
      integer, intent(in) :: child !! Child node number whose unique incoming edge is searched.
      logical, intent(out) :: found !! True when the requested child occurs in the edge matrix.
      real(dp) :: value
      integer :: i

      value = 0.0_dp
      found = .false.
      if (.not. allocated(tree%edge_length)) return
      do i = 1, tree%nedge()
         if (tree%edge(i, 2) /= child) cycle
         value = tree%edge_length(i)
         found = .true.
         return
      end do
   end function edge_length_to_child

   subroutine check(condition, label, failures)
      logical, intent(in) :: condition !! Boolean test condition that must be true for the check to pass.
      character(len=*), intent(in) :: label !! Short diagnostic label printed if the check fails.
      integer, intent(inout) :: failures !! Running failure count incremented when `condition` is false.

      if (.not. condition) then
         failures = failures + 1
         print '(a,1x,a)', 'FAIL:', trim(label)
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual !! Computed scalar value under test.
      real(dp), intent(in) :: expected !! Reference scalar value for the deterministic check.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute error.
      character(len=*), intent(in) :: label !! Short diagnostic label printed if the check fails.
      integer, intent(inout) :: failures !! Running failure count incremented when the tolerance is exceeded.

      call check(abs(actual - expected) <= tolerance, label, failures)
   end subroutine check_close

end program test_ape
