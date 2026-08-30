module e1071_svm_io
    use e1071_kinds, only: dp
    use e1071_svm, only: svm_model, svm_options, svm_pair_model, svm_c_classification, &
                         svm_nu_classification, svm_one_classification, svm_eps_regression, &
                         svm_nu_regression, svm_linear, svm_polynomial, svm_radial, svm_sigmoid
    implicit none
    private

    public :: svm_write_libsvm, svm_read_libsvm

contains

    subroutine svm_write_libsvm(model, svm_file, scale_file, yscale_file)
        type(svm_model), intent(in) :: model !! Fitted SVM whose prediction representation is written in LIBSVM text format.
        character(len=*), intent(in) :: svm_file !! Output LIBSVM model pathname, replaced when it already exists.
        character(len=*), intent(in), optional :: scale_file !! Optional predictor center/scale sidecar, one row per feature.
        character(len=*), intent(in), optional :: yscale_file !! Optional response scaling sidecar for epsilon- or nu-SVR models.
        real(dp), allocatable :: support_vectors(:, :)
        real(dp), allocatable :: coefficients(:, :)
        integer, allocatable :: nsv(:)
        integer :: nr_class
        integer :: total_sv
        integer :: unit
        integer :: p
        integer :: i
        integer :: j
        integer :: k

        call validate_serializable_model(model)
        call expand_libsvm_layout(model, support_vectors, coefficients, nsv, nr_class)
        total_sv = size(support_vectors, 1)

        open(newunit=unit, file=svm_file, status="replace", action="write")
        write(unit, '(a,1x,a)') "svm_type", trim(svm_type_name(model%options%svm_type))
        write(unit, '(a,1x,a)') "kernel_type", trim(kernel_name(model%options%kernel))
        if (model%options%kernel == svm_polynomial) write(unit, '(a,1x,i0)') "degree", model%options%degree
        if (model%options%kernel == svm_polynomial .or. model%options%kernel == svm_radial .or. &
            model%options%kernel == svm_sigmoid) then
            write(unit, '(a,1x,g0.17)') "gamma", model%options%gamma
        end if
        if (model%options%kernel == svm_polynomial .or. model%options%kernel == svm_sigmoid) then
            write(unit, '(a,1x,g0.17)') "coef0", model%options%coef0
        end if
        write(unit, '(a,1x,i0)') "nr_class", nr_class
        write(unit, '(a,1x,i0)') "total_sv", total_sv
        write(unit, '(a)', advance='no') "rho"
        do p = 1, size(model%pair)
            write(unit, '(1x,g0.17)', advance='no') model%pair(p)%rho
        end do
        write(unit, '()')
        if (is_classification(model)) then
            write(unit, '(a)', advance='no') "label"
            do i = 1, size(model%class_labels)
                write(unit, '(1x,i0)', advance='no') model%class_labels(i)
            end do
            write(unit, '()')
            if (model%probability_fitted) then
                write(unit, '(a)', advance='no') "probA"
                do p = 1, size(model%pair)
                    write(unit, '(1x,g0.17)', advance='no') model%pair(p)%prob_a
                end do
                write(unit, '()')
                write(unit, '(a)', advance='no') "probB"
                do p = 1, size(model%pair)
                    write(unit, '(1x,g0.17)', advance='no') model%pair(p)%prob_b
                end do
                write(unit, '()')
            end if
            write(unit, '(a)', advance='no') "nr_sv"
            do i = 1, size(nsv)
                write(unit, '(1x,i0)', advance='no') nsv(i)
            end do
            write(unit, '()')
        else if (is_regression(model) .and. model%probability_fitted) then
            write(unit, '(a,1x,g0.17)') "probA", model%svr_probability_sigma
        end if
        write(unit, '(a)') "SV"
        do i = 1, total_sv
            do j = 1, size(coefficients, 1)
                write(unit, '(g0.17,1x)', advance='no') coefficients(j, i)
            end do
            do k = 1, model%n_features
                if (abs(support_vectors(i, k)) <= 0.0_dp) cycle
                write(unit, '(i0,a,g0.8,1x)', advance='no') k, ':', support_vectors(i, k)
            end do
            write(unit, '()')
        end do
        close(unit)

        if (present(scale_file)) call write_predictor_scale(model, scale_file)
        if (present(yscale_file)) call write_response_scale(model, yscale_file)
    end subroutine svm_write_libsvm

    subroutine svm_read_libsvm(svm_file, model, scale_file, yscale_file, n_features)
        character(len=*), intent(in) :: svm_file !! Existing LIBSVM model pathname to parse into the native prediction model.
        type(svm_model), intent(out) :: model !! Native SVM reconstructed from the serialized support-vector decision functions.
        character(len=*), intent(in), optional :: scale_file !! Optional predictor center/scale sidecar applied on prediction.
        character(len=*), intent(in), optional :: yscale_file !! Optional response center/scale sidecar for loaded SVR predictions.
        integer, intent(in), optional :: n_features !! Optional feature count when the model omits trailing all-zero variables.
        character(len=65536) :: line
        character(len=128) :: key
        character(len=128) :: value
        real(dp), allocatable :: rho(:)
        real(dp), allocatable :: prob_a(:)
        real(dp), allocatable :: prob_b(:)
        real(dp), allocatable :: coefficients(:, :)
        real(dp), allocatable :: support_vectors(:, :)
        integer, allocatable :: labels(:)
        integer, allocatable :: nsv(:)
        integer :: nr_class
        integer :: total_sv
        integer :: max_feature
        integer :: unit
        integer :: ios
        integer :: ncoef
        integer :: row
        logical :: in_sv
        logical :: have_prob_a
        logical :: have_prob_b

        model%options = svm_options()
        nr_class = 0
        total_sv = 0
        max_feature = 0
        have_prob_a = .false.
        have_prob_b = .false.
        in_sv = .false.
        open(newunit=unit, file=svm_file, status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "svm_read_libsvm: unable to open model file"
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            if (.not. in_sv) then
                call split_key_value(line, key, value)
                select case (trim(key))
                case ("svm_type")
                    model%options%svm_type = parse_svm_type(trim(value))
                case ("kernel_type")
                    model%options%kernel = parse_kernel(trim(value))
                case ("degree")
                    read(value, *, iostat=ios) model%options%degree
                    if (ios /= 0) error stop "svm_read_libsvm: invalid degree"
                case ("gamma")
                    read(value, *, iostat=ios) model%options%gamma
                    if (ios /= 0) error stop "svm_read_libsvm: invalid gamma"
                case ("coef0")
                    read(value, *, iostat=ios) model%options%coef0
                    if (ios /= 0) error stop "svm_read_libsvm: invalid coef0"
                case ("nr_class")
                    read(value, *, iostat=ios) nr_class
                    if (ios /= 0 .or. nr_class < 1) error stop "svm_read_libsvm: invalid nr_class"
                case ("total_sv")
                    read(value, *, iostat=ios) total_sv
                    if (ios /= 0 .or. total_sv < 0) error stop "svm_read_libsvm: invalid total_sv"
                case ("rho")
                    call parse_real_vector(value, rho)
                case ("label")
                    call parse_integer_vector(value, labels)
                case ("probA")
                    call parse_real_vector(value, prob_a)
                    have_prob_a = .true.
                case ("probB")
                    call parse_real_vector(value, prob_b)
                    have_prob_b = .true.
                case ("nr_sv")
                    call parse_integer_vector(value, nsv)
                case ("SV")
                    in_sv = .true.
                case default
                    error stop "svm_read_libsvm: unsupported model header"
                end select
            else
                ncoef = max(1, nr_class - 1)
                call scan_sv_line(line, ncoef, max_feature)
            end if
        end do
        close(unit)

        if (nr_class < 1 .or. total_sv < 0) error stop "svm_read_libsvm: incomplete header"
        if (present(scale_file)) max_feature = max(max_feature, count_scale_rows(scale_file))
        if (present(n_features)) then
            if (n_features < max_feature) error stop "svm_read_libsvm: n_features is smaller than serialized feature index"
            max_feature = n_features
        end if
        if (max_feature < 1) error stop "svm_read_libsvm: model feature count is unknown"
        ncoef = max(1, nr_class - 1)
        allocate(coefficients(ncoef, total_sv), support_vectors(total_sv, max_feature))
        coefficients = 0.0_dp
        support_vectors = 0.0_dp
        open(newunit=unit, file=svm_file, status="old", action="read")
        in_sv = .false.
        row = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (.not. in_sv) then
                call split_key_value(line, key, value)
                if (trim(key) == "SV") in_sv = .true.
                cycle
            end if
            if (len_trim(line) == 0) cycle
            row = row + 1
            if (row > total_sv) error stop "svm_read_libsvm: too many support-vector rows"
            call parse_sv_line(line, ncoef, coefficients(:, row), support_vectors(row, :))
        end do
        close(unit)
        if (row /= total_sv) error stop "svm_read_libsvm: support-vector count mismatch"

        call rebuild_native_model(model, nr_class, rho, labels, nsv, prob_a, prob_b, have_prob_a, have_prob_b, &
                                  coefficients, support_vectors)
        model%n_features = max_feature
        call set_default_scaling(model)
        if (present(scale_file)) call read_predictor_scale(model, scale_file)
        if (present(yscale_file)) call read_response_scale(model, yscale_file)
    end subroutine svm_read_libsvm

    subroutine expand_libsvm_layout(model, support_vectors, coefficients, nsv, nr_class)
        type(svm_model), intent(in) :: model !! Native fitted model whose pair expansions are converted to LIBSVM layout.
        real(dp), allocatable, intent(out) :: support_vectors(:, :) !! Global support-vector rows in class-block order.
        real(dp), allocatable, intent(out) :: coefficients(:, :) !! LIBSVM coefficient rows by global support-vector column.
        integer, allocatable, intent(out) :: nsv(:) !! Number of emitted support-vector rows associated with each class.
        integer, intent(out) :: nr_class !! LIBSVM class count; two for one-class and regression models.
        integer :: total_sv
        integer :: ncoef
        integer :: class_index
        integer :: pair_index
        integer :: a
        integer :: b
        integer :: i
        integer :: row

        if (.not. is_classification(model)) then
            nr_class = 2
            total_sv = size(model%pair(1)%coefficients)
            allocate(nsv(0), support_vectors(total_sv, model%n_features), coefficients(1, total_sv))
            support_vectors = model%pair(1)%support_vectors
            coefficients(1, :) = model%pair(1)%coefficients
            return
        end if

        nr_class = size(model%class_labels)
        ncoef = nr_class - 1
        allocate(nsv(nr_class))
        nsv = 0
        pair_index = 0
        do a = 1, nr_class - 1
            do b = a + 1, nr_class
                pair_index = pair_index + 1
                do i = 1, size(model%pair(pair_index)%coefficients)
                    if (model%pair(pair_index)%coefficients(i) > 0.0_dp) then
                        nsv(a) = nsv(a) + 1
                    else
                        nsv(b) = nsv(b) + 1
                    end if
                end do
            end do
        end do
        total_sv = sum(nsv)
        allocate(support_vectors(total_sv, model%n_features), coefficients(ncoef, total_sv))
        support_vectors = 0.0_dp
        coefficients = 0.0_dp
        row = 0
        do class_index = 1, nr_class
            pair_index = 0
            do a = 1, nr_class - 1
                do b = a + 1, nr_class
                    pair_index = pair_index + 1
                    do i = 1, size(model%pair(pair_index)%coefficients)
                        if (class_index == a .and. model%pair(pair_index)%coefficients(i) > 0.0_dp) then
                            row = row + 1
                            support_vectors(row, :) = model%pair(pair_index)%support_vectors(i, :)
                            coefficients(b - 1, row) = model%pair(pair_index)%coefficients(i)
                        else if (class_index == b .and. model%pair(pair_index)%coefficients(i) < 0.0_dp) then
                            row = row + 1
                            support_vectors(row, :) = model%pair(pair_index)%support_vectors(i, :)
                            coefficients(a, row) = model%pair(pair_index)%coefficients(i)
                        end if
                    end do
                end do
            end do
        end do
    end subroutine expand_libsvm_layout

    subroutine rebuild_native_model(model, nr_class, rho, labels, nsv, prob_a, prob_b, have_prob_a, have_prob_b, &
                                    coefficients, support_vectors)
        type(svm_model), intent(inout) :: model !! Model header already parsed and completed here with native pair expansions.
        integer, intent(in) :: nr_class !! Class count declared by the LIBSVM file.
        real(dp), intent(in) :: rho(:) !! Pairwise or scalar LIBSVM rho values in standard pair order.
        integer, allocatable, intent(in) :: labels(:) !! Classification labels when present in the serialized model header.
        integer, allocatable, intent(in) :: nsv(:) !! Per-class support-vector counts when present in the model header.
        real(dp), allocatable, intent(in) :: prob_a(:) !! Pairwise sigmoid A values or scalar SVR probability scale when present.
        real(dp), allocatable, intent(in) :: prob_b(:) !! Pairwise sigmoid B values for probabilistic classification when present.
        logical, intent(in) :: have_prob_a !! True when the model header contained a probA record.
        logical, intent(in) :: have_prob_b !! True when the model header contained a probB record.
        real(dp), intent(in) :: coefficients(:, :) !! Parsed LIBSVM coefficient rows by support-vector column.
        real(dp), intent(in) :: support_vectors(:, :) !! Parsed global support-vector matrix in scaled predictor units.
        integer :: npair
        integer :: pair_index
        integer :: a
        integer :: b
        integer :: start_a
        integer :: start_b
        integer :: end_a
        integer :: end_b
        integer :: i
        integer :: count_pair
        integer :: pos

        if (.not. is_classification(model)) then
            if (size(rho) /= 1) error stop "svm_read_libsvm: scalar model must contain one rho"
            allocate(model%pair(1))
            allocate(model%pair(1)%support_vectors(size(support_vectors, 1), size(support_vectors, 2)))
            allocate(model%pair(1)%coefficients(size(support_vectors, 1)))
            model%pair(1)%support_vectors = support_vectors
            model%pair(1)%coefficients = coefficients(1, :)
            model%pair(1)%rho = rho(1)
            if (is_regression(model) .and. have_prob_a) then
                if (.not. allocated(prob_a)) error stop "svm_read_libsvm: missing probA values"
                if (size(prob_a) /= 1) error stop "svm_read_libsvm: invalid SVR probA length"
                model%svr_probability_sigma = prob_a(1)
                model%probability_fitted = .true.
            end if
            return
        end if

        if (.not. allocated(labels) .or. .not. allocated(nsv)) then
            error stop "svm_read_libsvm: classification model requires label and nr_sv records"
        end if
        if (size(labels) /= nr_class .or. size(nsv) /= nr_class) then
            error stop "svm_read_libsvm: classification header length mismatch"
        end if
        if (sum(nsv) /= size(support_vectors, 1)) error stop "svm_read_libsvm: nr_sv does not sum to total_sv"
        npair = nr_class * (nr_class - 1) / 2
        if (size(rho) /= npair) error stop "svm_read_libsvm: rho count mismatch"
        allocate(model%class_labels(nr_class), model%pair(npair))
        model%class_labels = labels
        pair_index = 0
        start_a = 1
        do a = 1, nr_class - 1
            end_a = start_a + nsv(a) - 1
            start_b = end_a + 1
            do b = a + 1, nr_class
                if (b > a + 1) start_b = start_b + nsv(b - 1)
                end_b = start_b + nsv(b) - 1
                pair_index = pair_index + 1
                count_pair = 0
                do i = start_a, end_a
                    if (abs(coefficients(b - 1, i)) > 0.0_dp) count_pair = count_pair + 1
                end do
                do i = start_b, end_b
                    if (abs(coefficients(a, i)) > 0.0_dp) count_pair = count_pair + 1
                end do
                allocate(model%pair(pair_index)%support_vectors(count_pair, size(support_vectors, 2)))
                allocate(model%pair(pair_index)%coefficients(count_pair))
                pos = 0
                do i = start_a, end_a
                    if (abs(coefficients(b - 1, i)) <= 0.0_dp) cycle
                    pos = pos + 1
                    model%pair(pair_index)%support_vectors(pos, :) = support_vectors(i, :)
                    model%pair(pair_index)%coefficients(pos) = coefficients(b - 1, i)
                end do
                do i = start_b, end_b
                    if (abs(coefficients(a, i)) <= 0.0_dp) cycle
                    pos = pos + 1
                    model%pair(pair_index)%support_vectors(pos, :) = support_vectors(i, :)
                    model%pair(pair_index)%coefficients(pos) = coefficients(a, i)
                end do
                model%pair(pair_index)%rho = rho(pair_index)
                model%pair(pair_index)%positive_label = labels(a)
                model%pair(pair_index)%negative_label = labels(b)
                if (have_prob_a .and. have_prob_b) then
                    if (.not. allocated(prob_a) .or. .not. allocated(prob_b)) then
                        error stop "svm_read_libsvm: incomplete probability arrays"
                    end if
                    model%pair(pair_index)%prob_a = prob_a(pair_index)
                    model%pair(pair_index)%prob_b = prob_b(pair_index)
                end if
            end do
            start_a = end_a + 1
        end do
        model%probability_fitted = have_prob_a .and. have_prob_b
    end subroutine rebuild_native_model

    subroutine validate_serializable_model(model)
        type(svm_model), intent(in) :: model !! Native SVM checked for fields required by LIBSVM prediction serialization.

        if (.not. allocated(model%pair)) error stop "svm_write_libsvm: model has no support-vector pairs"
        if (model%n_features < 1) error stop "svm_write_libsvm: invalid feature count"
        if (model%options%kernel < svm_linear .or. model%options%kernel > svm_sigmoid) then
            error stop "svm_write_libsvm: unsupported kernel"
        end if
        if (is_classification(model)) then
            if (.not. allocated(model%class_labels)) error stop "svm_write_libsvm: classification labels are missing"
            if (size(model%pair) /= size(model%class_labels) * (size(model%class_labels) - 1) / 2) then
                error stop "svm_write_libsvm: inconsistent pair count"
            end if
        else if (size(model%pair) /= 1) then
            error stop "svm_write_libsvm: scalar SVM must contain one support-vector expansion"
        end if
    end subroutine validate_serializable_model

    subroutine set_default_scaling(model)
        type(svm_model), intent(inout) :: model !! Loaded model receiving identity predictor and response scaling defaults.

        allocate(model%center(model%n_features), model%scale_value(model%n_features), model%scale_mask(model%n_features))
        model%center = 0.0_dp
        model%scale_value = 1.0_dp
        model%scale_mask = .false.
        model%y_center = 0.0_dp
        model%y_scale = 1.0_dp
    end subroutine set_default_scaling

    subroutine write_predictor_scale(model, filename)
        type(svm_model), intent(in) :: model !! Model whose predictor centers and scales are written one row per feature.
        character(len=*), intent(in) :: filename !! Output sidecar pathname containing two numeric columns: center and scale.
        integer :: unit
        integer :: j

        open(newunit=unit, file=filename, status="replace", action="write")
        do j = 1, model%n_features
            if (allocated(model%scale_mask) .and. model%scale_mask(j)) then
                write(unit, '(g0.17,1x,g0.17)') model%center(j), model%scale_value(j)
            else
                write(unit, '(g0.17,1x,g0.17)') 0.0_dp, 1.0_dp
            end if
        end do
        close(unit)
    end subroutine write_predictor_scale

    function count_scale_rows(filename) result(nrow)
        character(len=*), intent(in) :: filename !! Predictor scaling sidecar whose nonempty rows are counted.
        integer :: nrow
        character(len=1024) :: line
        integer :: unit
        integer :: ios

        nrow = 0
        open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "svm_read_libsvm: unable to open predictor scale file"
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0) nrow = nrow + 1
        end do
        close(unit)
    end function count_scale_rows

    subroutine read_predictor_scale(model, filename)
        type(svm_model), intent(inout) :: model !! Loaded model whose predictor preprocessing is replaced from the sidecar.
        character(len=*), intent(in) :: filename !! Sidecar pathname with exactly one center/scale row per model feature.
        integer :: unit
        integer :: ios
        integer :: j
        real(dp) :: center
        real(dp) :: scale_value

        open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "svm_read_libsvm: unable to open predictor scale file"
        do j = 1, model%n_features
            read(unit, *, iostat=ios) center, scale_value
            if (ios /= 0) error stop "svm_read_libsvm: predictor scale row count mismatch"
            if (scale_value <= 0.0_dp) error stop "svm_read_libsvm: predictor scale must be positive"
            model%center(j) = center
            model%scale_value(j) = scale_value
            model%scale_mask(j) = abs(center) > 0.0_dp .or. abs(scale_value - 1.0_dp) > 0.0_dp
        end do
        read(unit, *, iostat=ios) center, scale_value
        if (ios == 0) error stop "svm_read_libsvm: predictor scale contains extra rows"
        close(unit)
    end subroutine read_predictor_scale

    subroutine write_response_scale(model, filename)
        type(svm_model), intent(in) :: model !! Regression model whose response center and scale are written for prediction.
        character(len=*), intent(in) :: filename !! Output sidecar pathname containing one center/scale row.
        integer :: unit

        if (.not. is_regression(model)) error stop "svm_write_libsvm: yscale_file is valid only for regression"
        open(newunit=unit, file=filename, status="replace", action="write")
        write(unit, '(g0.17,1x,g0.17)') model%y_center, model%y_scale
        close(unit)
    end subroutine write_response_scale

    subroutine read_response_scale(model, filename)
        type(svm_model), intent(inout) :: model !! Loaded regression model receiving response center and scale from the sidecar.
        character(len=*), intent(in) :: filename !! Sidecar pathname containing one numeric center/scale row.
        integer :: unit
        integer :: ios

        if (.not. is_regression(model)) error stop "svm_read_libsvm: yscale_file is valid only for regression"
        open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "svm_read_libsvm: unable to open response scale file"
        read(unit, *, iostat=ios) model%y_center, model%y_scale
        if (ios /= 0 .or. model%y_scale <= 0.0_dp) error stop "svm_read_libsvm: invalid response scale file"
        close(unit)
    end subroutine read_response_scale

    subroutine scan_sv_line(line, ncoef, max_feature)
        character(len=*), intent(in) :: line !! Serialized support-vector row scanned only to discover its largest feature index.
        integer, intent(in) :: ncoef !! Number of leading coefficient tokens before sparse `index:value` features.
        integer, intent(inout) :: max_feature !! Largest positive feature index seen across all scanned support-vector rows.
        character(len=256) :: token
        integer :: position
        integer :: i
        integer :: colon
        integer :: index_value
        integer :: ios

        position = 1
        do i = 1, ncoef
            call next_token(line, position, token)
            if (len_trim(token) == 0) error stop "svm_read_libsvm: missing support-vector coefficient"
        end do
        do
            call next_token(line, position, token)
            if (len_trim(token) == 0) exit
            colon = index(token, ':')
            if (colon <= 1) error stop "svm_read_libsvm: malformed support-vector feature"
            read(token(:colon - 1), *, iostat=ios) index_value
            if (ios /= 0 .or. index_value < 1) error stop "svm_read_libsvm: invalid support-vector feature index"
            max_feature = max(max_feature, index_value)
        end do
    end subroutine scan_sv_line

    subroutine parse_sv_line(line, ncoef, coefficients, support_vector)
        character(len=*), intent(in) :: line !! Serialized support-vector row containing coefficients followed by sparse features.
        integer, intent(in) :: ncoef !! Number of leading LIBSVM coefficient fields to parse.
        real(dp), intent(out) :: coefficients(:) !! Parsed coefficient fields; shape must equal `ncoef`.
        real(dp), intent(out) :: support_vector(:) !! Dense support-vector row reconstructed from one-based sparse features.
        character(len=256) :: token
        integer :: position
        integer :: i
        integer :: colon
        integer :: index_value
        integer :: ios
        real(dp) :: feature_value

        if (size(coefficients) /= ncoef) error stop "svm_read_libsvm: coefficient shape mismatch"
        coefficients = 0.0_dp
        support_vector = 0.0_dp
        position = 1
        do i = 1, ncoef
            call next_token(line, position, token)
            read(token, *, iostat=ios) coefficients(i)
            if (ios /= 0) error stop "svm_read_libsvm: invalid support-vector coefficient"
        end do
        do
            call next_token(line, position, token)
            if (len_trim(token) == 0) exit
            colon = index(token, ':')
            if (colon <= 1) error stop "svm_read_libsvm: malformed support-vector feature"
            read(token(:colon - 1), *, iostat=ios) index_value
            if (ios /= 0 .or. index_value < 1 .or. index_value > size(support_vector)) then
                error stop "svm_read_libsvm: support-vector feature index out of range"
            end if
            read(token(colon + 1:), *, iostat=ios) feature_value
            if (ios /= 0) error stop "svm_read_libsvm: invalid support-vector feature value"
            support_vector(index_value) = feature_value
        end do
    end subroutine parse_sv_line

    subroutine split_key_value(line, key, value)
        character(len=*), intent(in) :: line !! Header line split at its first whitespace separator.
        character(len=*), intent(out) :: key !! Header keyword such as `svm_type`, `rho`, or `SV`.
        character(len=*), intent(out) :: value !! Remaining trimmed header text after the keyword; empty for `SV`.
        integer :: space

        key = ""
        value = ""
        if (len_trim(line) == 0) return
        space = scan(trim(line), ' ' // achar(9))
        if (space == 0) then
            key = trim(line)
        else
            key = trim(line(:space - 1))
            value = adjustl(line(space + 1:))
        end if
    end subroutine split_key_value

    subroutine parse_real_vector(text, values)
        character(len=*), intent(in) :: text !! Whitespace-separated real-valued header payload to decode.
        real(dp), allocatable, intent(out) :: values(:) !! Parsed real values in their serialized order.
        character(len=256) :: token
        integer :: position
        integer :: count_value
        integer :: i
        integer :: ios

        position = 1
        count_value = 0
        do
            call next_token(text, position, token)
            if (len_trim(token) == 0) exit
            count_value = count_value + 1
        end do
        allocate(values(count_value))
        position = 1
        do i = 1, count_value
            call next_token(text, position, token)
            read(token, *, iostat=ios) values(i)
            if (ios /= 0) error stop "svm_read_libsvm: invalid real header value"
        end do
    end subroutine parse_real_vector

    subroutine parse_integer_vector(text, values)
        character(len=*), intent(in) :: text !! Whitespace-separated integer-valued header payload to decode.
        integer, allocatable, intent(out) :: values(:) !! Parsed integers in their serialized order.
        character(len=256) :: token
        integer :: position
        integer :: count_value
        integer :: i
        integer :: ios

        position = 1
        count_value = 0
        do
            call next_token(text, position, token)
            if (len_trim(token) == 0) exit
            count_value = count_value + 1
        end do
        allocate(values(count_value))
        position = 1
        do i = 1, count_value
            call next_token(text, position, token)
            read(token, *, iostat=ios) values(i)
            if (ios /= 0) error stop "svm_read_libsvm: invalid integer header value"
        end do
    end subroutine parse_integer_vector

    subroutine next_token(text, position, token)
        character(len=*), intent(in) :: text !! Source text containing whitespace-delimited tokens.
        integer, intent(inout) :: position !! One-based scan position advanced past the returned token.
        character(len=*), intent(out) :: token !! Next token, or an empty string when no token remains.
        integer :: start
        integer :: finish
        integer :: n

        token = ""
        n = len_trim(text)
        do while (position <= n)
            if (text(position:position) /= ' ' .and. text(position:position) /= achar(9)) exit
            position = position + 1
        end do
        if (position > n) return
        start = position
        finish = start
        do while (finish <= n)
            if (text(finish:finish) == ' ' .or. text(finish:finish) == achar(9)) exit
            finish = finish + 1
        end do
        token = text(start:finish - 1)
        position = finish
    end subroutine next_token

    pure function is_classification(model) result(answer)
        type(svm_model), intent(in) :: model !! SVM inspected to determine whether it is C-SVC or nu-SVC.
        logical :: answer

        answer = model%options%svm_type == svm_c_classification .or. &
                 model%options%svm_type == svm_nu_classification
    end function is_classification

    pure function is_regression(model) result(answer)
        type(svm_model), intent(in) :: model !! SVM inspected to determine whether it is epsilon-SVR or nu-SVR.
        logical :: answer

        answer = model%options%svm_type == svm_eps_regression .or. model%options%svm_type == svm_nu_regression
    end function is_regression

    pure function svm_type_name(svm_type) result(name)
        integer, intent(in) :: svm_type !! Native SVM type constant converted to the standard LIBSVM header spelling.
        character(len=16) :: name

        select case (svm_type)
        case (svm_c_classification)
            name = "c_svc"
        case (svm_nu_classification)
            name = "nu_svc"
        case (svm_one_classification)
            name = "one_class"
        case (svm_eps_regression)
            name = "epsilon_svr"
        case (svm_nu_regression)
            name = "nu_svr"
        case default
            name = ""
        end select
    end function svm_type_name

    pure function kernel_name(kernel) result(name)
        integer, intent(in) :: kernel !! Native kernel constant converted to the standard LIBSVM header spelling.
        character(len=16) :: name

        select case (kernel)
        case (svm_linear)
            name = "linear"
        case (svm_polynomial)
            name = "polynomial"
        case (svm_radial)
            name = "rbf"
        case (svm_sigmoid)
            name = "sigmoid"
        case default
            name = ""
        end select
    end function kernel_name

    pure function parse_svm_type(name) result(svm_type)
        character(len=*), intent(in) :: name !! LIBSVM svm_type header token converted to a native type constant.
        integer :: svm_type

        select case (trim(name))
        case ("c_svc")
            svm_type = svm_c_classification
        case ("nu_svc")
            svm_type = svm_nu_classification
        case ("one_class")
            svm_type = svm_one_classification
        case ("epsilon_svr")
            svm_type = svm_eps_regression
        case ("nu_svr")
            svm_type = svm_nu_regression
        case default
            error stop "svm_read_libsvm: unsupported svm_type"
        end select
    end function parse_svm_type

    pure function parse_kernel(name) result(kernel)
        character(len=*), intent(in) :: name !! LIBSVM kernel_type header token converted to a native kernel constant.
        integer :: kernel

        select case (trim(name))
        case ("linear")
            kernel = svm_linear
        case ("polynomial")
            kernel = svm_polynomial
        case ("rbf")
            kernel = svm_radial
        case ("sigmoid")
            kernel = svm_sigmoid
        case default
            error stop "svm_read_libsvm: unsupported kernel_type"
        end select
    end function parse_kernel

end module e1071_svm_io
