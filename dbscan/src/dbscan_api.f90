module dbscan_api
    use dbscan_kinds, only : dp
    use dbscan_types, only : knn_result, frnn_result, clustering_result, optics_result, dendrogram_result
    use dbscan_types, only : hdbscan_result, dbcv_result
    use dbscan_nn, only : knn, knn_query, frnn, frnn_query, knn_from_dist, frnn_from_dist
    use dbscan_nn, only : adjacencylist_knn, knn_dist
    use dbscan_graph, only : connected_components_knn, connected_components_frnn
    use dbscan_cluster, only : dbscan, dbscan_from_frnn, is_corepoint
    use dbscan_snn, only : snn, snnclust, jpclust
    use dbscan_outlier, only : lof, pointdensity
    use dbscan_optics, only : optics, extract_dbscan, extract_xi
    use dbscan_hdbscan, only : coredist, mrdist, mst_dense, hdbscan, glosh, extract_fosc
    use dbscan_dbcv, only : dbcv
    use dbscan_conversion, only : optics_to_dendrogram, dendrogram_to_reachability
    use dbscan_utils, only : unique_positive_count, count_zero
    implicit none
    private

    public :: dp
    public :: knn_result, frnn_result, clustering_result, optics_result, dendrogram_result
    public :: hdbscan_result, dbcv_result
    public :: knn, knn_query, frnn, frnn_query, knn_from_dist, frnn_from_dist
    public :: adjacencylist_knn, adjacencylist_frnn, knn_dist
    public :: comps_knn, comps_frnn
    public :: dbscan, dbscan_from_frnn, is_corepoint
    public :: snn, snnclust, jpclust
    public :: lof, pointdensity
    public :: optics, extract_dbscan, extract_xi
    public :: coredist, mrdist, mst_dense, hdbscan, glosh, extract_fosc
    public :: dbcv
    public :: optics_to_dendrogram, dendrogram_to_reachability
    public :: ncluster, nnoise

    interface ncluster
        module procedure ncluster_clustering
        module procedure ncluster_optics
        module procedure ncluster_hdbscan
    end interface ncluster

    interface nnoise
        module procedure nnoise_clustering
        module procedure nnoise_optics
        module procedure nnoise_hdbscan
    end interface nnoise

contains

    function adjacencylist_frnn(nn) result(adj)
        type(frnn_result), intent(in) :: nn !! Fixed-radius neighbor object returned unchanged as an adjacency list.
        type(frnn_result) :: adj
        adj = nn
    end function adjacencylist_frnn

    subroutine comps_knn(nn, labels, mutual)
        type(knn_result), intent(in) :: nn !! k-nearest-neighbor graph.
        integer, allocatable, intent(out) :: labels(:) !! Connected-component labels.
        logical, intent(in), optional :: mutual !! Require reciprocal edges when true.
        call connected_components_knn(nn, labels, mutual)
    end subroutine comps_knn

    subroutine comps_frnn(nn, labels, mutual)
        type(frnn_result), intent(in) :: nn !! Fixed-radius adjacency graph.
        integer, allocatable, intent(out) :: labels(:) !! Connected-component labels.
        logical, intent(in), optional :: mutual !! Require reciprocal edges when true.
        call connected_components_frnn(nn, labels, mutual)
    end subroutine comps_frnn

    pure integer function ncluster_clustering(x) result(n)
        type(clustering_result), intent(in) :: x !! Flat clustering whose positive labels are counted.
        n = unique_positive_count(x%cluster)
    end function ncluster_clustering

    pure integer function ncluster_optics(x) result(n)
        type(optics_result), intent(in) :: x !! OPTICS result with an extracted cluster vector.
        if (allocated(x%cluster)) then
            n = unique_positive_count(x%cluster)
        else
            n = 0
        end if
    end function ncluster_optics

    pure integer function ncluster_hdbscan(x) result(n)
        type(hdbscan_result), intent(in) :: x !! HDBSCAN result whose positive labels are counted.
        n = unique_positive_count(x%cluster)
    end function ncluster_hdbscan

    pure integer function nnoise_clustering(x) result(n)
        type(clustering_result), intent(in) :: x !! Flat clustering where zero denotes noise.
        n = count_zero(x%cluster)
    end function nnoise_clustering

    pure integer function nnoise_optics(x) result(n)
        type(optics_result), intent(in) :: x !! OPTICS result with an extracted cluster vector.
        if (allocated(x%cluster)) then
            n = count_zero(x%cluster)
        else
            n = 0
        end if
    end function nnoise_optics

    pure integer function nnoise_hdbscan(x) result(n)
        type(hdbscan_result), intent(in) :: x !! HDBSCAN result where zero denotes noise.
        n = count_zero(x%cluster)
    end function nnoise_hdbscan

end module dbscan_api
