module rpart_types
   use rpart_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: RPART_ANOVA = 1
   integer, parameter, public :: RPART_POISSON = 2
   integer, parameter, public :: RPART_CLASS = 3
   integer, parameter, public :: RPART_EXP = 4
   integer, parameter, public :: RPART_GINI = 1
   integer, parameter, public :: RPART_INFORMATION = 2
   integer, parameter, public :: RPART_DEVIANCE = 1
   integer, parameter, public :: RPART_SQRT = 2
   integer, parameter, public :: RPART_LEFT = -1
   integer, parameter, public :: RPART_MISSING = 0
   integer, parameter, public :: RPART_RIGHT = 1

   type, public :: rpart_control
      integer :: minsplit = 20
      integer :: minbucket = 7
      real(dp) :: cp = 0.01_dp
      integer :: maxcompete = 4
      integer :: maxsurrogate = 5
      integer :: usesurrogate = 2
      integer :: surrogatestyle = 0
      integer :: maxdepth = 30
      integer :: xval = 10
   end type rpart_control

   type, public :: rpart_split
      integer :: var = 0
      integer :: ncat = 0
      integer :: count = 0
      real(dp) :: improve = 0.0_dp
      real(dp) :: adj = 0.0_dp
      real(dp) :: spoint = 0.0_dp
      integer :: direction = RPART_LEFT
      integer, allocatable :: csplit(:)
   end type rpart_split

   type, public :: rpart_node
      integer :: id = 1
      integer :: depth = 0
      integer :: nobs = 0
      integer :: nfinal = 0
      real(dp) :: sum_wt = 0.0_dp
      real(dp) :: risk = 0.0_dp
      real(dp) :: complexity = 0.0_dp
      real(dp), allocatable :: response(:)
      type(rpart_split), allocatable :: primary(:)
      type(rpart_split), allocatable :: surrogate(:)
      integer :: lastsurrogate = RPART_MISSING
      type(rpart_node), allocatable :: left
      type(rpart_node), allocatable :: right
   end type rpart_node

   type, public :: rpart_cp_row
      real(dp) :: cp = 0.0_dp
      integer :: nsplit = 0
      real(dp) :: rel_error = 0.0_dp
      real(dp) :: xerror = 0.0_dp
      real(dp) :: xstd = 0.0_dp
   end type rpart_cp_row

   type, public :: rpart_model
      integer :: method = RPART_ANOVA
      integer :: nvar = 0
      integer :: nclass = 0
      integer :: nresp = 1
      integer :: poisson_method = RPART_DEVIANCE
      integer :: split_rule = RPART_GINI
      real(dp) :: poisson_shrink = 1.0_dp
      real(dp) :: poisson_alpha = 0.0_dp
      real(dp) :: poisson_beta = 0.0_dp
      real(dp) :: root_risk = 0.0_dp
      real(dp) :: total_weight = 0.0_dp
      type(rpart_control) :: control
      integer, allocatable :: ncat(:)
      real(dp), allocatable :: vcost(:)
      real(dp), allocatable :: prior(:)
      real(dp), allocatable :: altered_prior(:)
      real(dp), allocatable :: class_freq(:)
      real(dp), allocatable :: loss(:,:)
      type(rpart_node), allocatable :: root
      type(rpart_cp_row), allocatable :: cptable(:)
      real(dp), allocatable :: variable_importance(:)
      integer, allocatable :: where(:)
      real(dp), allocatable :: fitted(:)
   end type rpart_model

end module rpart_types
