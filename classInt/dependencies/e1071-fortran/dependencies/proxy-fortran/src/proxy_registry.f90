module proxy_registry
    use proxy_kinds, only: dp
    use proxy_utils, only: proxy_convert_default, proxy_normalize_name
    implicit none
    private

    integer, parameter :: max_registry_entries = 64
    integer, parameter :: registry_name_length = 64

    abstract interface
        function proxy_numeric_callback(x, y, p) result(value)
            import dp
            real(dp), intent(in) :: x(:) !! First observation vector supplied by the proximity engine.
            real(dp), intent(in) :: y(:) !! Second observation vector supplied by the proximity engine.
            real(dp), intent(in), optional :: p(:) !! Optional user parameter vector forwarded unchanged by the engine.
            real(dp) :: value
        end function proxy_numeric_callback

        function proxy_binary_callback(a, b, c, d, n) result(value)
            import dp
            integer, intent(in) :: a !! Count of TRUE/TRUE component pairs.
            integer, intent(in) :: b !! Count of TRUE/FALSE component pairs.
            integer, intent(in) :: c !! Count of FALSE/TRUE component pairs.
            integer, intent(in) :: d !! Count of FALSE/FALSE component pairs.
            integer, intent(in) :: n !! Total number of valid binary component pairs.
            real(dp) :: value
        end function proxy_binary_callback
    end interface

    type :: numeric_registry_entry
        character(len=registry_name_length) :: name = ''
        logical :: is_distance = .true.
        integer :: conversion = proxy_convert_default
        procedure(proxy_numeric_callback), pointer, nopass :: fun => null()
    end type numeric_registry_entry

    type :: binary_registry_entry
        character(len=registry_name_length) :: name = ''
        logical :: is_distance = .false.
        integer :: conversion = proxy_convert_default
        procedure(proxy_binary_callback), pointer, nopass :: fun => null()
    end type binary_registry_entry

    type(numeric_registry_entry), save :: numeric_entries(max_registry_entries)
    type(binary_registry_entry), save :: binary_entries(max_registry_entries)
    integer, save :: numeric_count = 0
    integer, save :: binary_count = 0

    public :: proxy_numeric_callback, proxy_binary_callback
    public :: proxy_register_numeric, proxy_register_binary, proxy_registry_clear
    public :: proxy_registry_numeric_names, proxy_registry_binary_names
    public :: proxy_lookup_numeric, proxy_lookup_binary

contains

    subroutine proxy_register_numeric(name, function_pointer, is_distance, conversion, status)
        character(len=*), intent(in) :: name !! Unique method name or alias used for case-insensitive lookup after punctuation
        !! normalization.
        procedure(proxy_numeric_callback), pointer, intent(in) :: function_pointer !! Numeric callback accepting two observation
        !! vectors and an optional parameter vector.
        logical, intent(in) :: is_distance !! True when the callback returns a distance/dissimilarity; false when it returns a
        !! similarity.
        integer, intent(in), optional :: conversion !! Optional conversion code used when callers request the opposite proximity
        !! type; defaults to proxy's standard conversion.
        integer, intent(out), optional :: status !! Zero on success; nonzero if the registry is full or the supplied normalized
        !! name is empty.
        character(len=:), allocatable :: key
        integer :: i
        integer :: slot

        if (present(status)) status = 0
        key = proxy_normalize_name(name)
        if (len(key) == 0) then
            if (present(status)) status = 2
            return
        end if
        slot = 0
        do i = 1, numeric_count
            if (trim(numeric_entries(i)%name) == key) then
                slot = i
                exit
            end if
        end do
        if (slot == 0) then
            if (numeric_count >= max_registry_entries) then
                if (present(status)) status = 1
                return
            end if
            numeric_count = numeric_count + 1
            slot = numeric_count
        end if
        numeric_entries(slot)%name = key
        numeric_entries(slot)%is_distance = is_distance
        numeric_entries(slot)%conversion = proxy_convert_default
        if (present(conversion)) numeric_entries(slot)%conversion = conversion
        numeric_entries(slot)%fun => function_pointer
    end subroutine proxy_register_numeric

    subroutine proxy_register_binary(name, function_pointer, is_distance, conversion, status)
        character(len=*), intent(in) :: name !! Unique binary-method name used for normalized case-insensitive lookup.
        procedure(proxy_binary_callback), pointer, intent(in) :: function_pointer !! Binary callback receiving the precomputed
        !! `(a,b,c,d,n)` contingency counts.
        logical, intent(in) :: is_distance !! True when the binary callback yields a distance; false when it yields a similarity.
        integer, intent(in), optional :: conversion !! Optional conversion code for requests of the opposite proximity type;
        !! defaults to standard proxy conversion.
        integer, intent(out), optional :: status !! Zero on success; nonzero when the registry is full or the normalized name is
        !! empty.
        character(len=:), allocatable :: key
        integer :: i
        integer :: slot

        if (present(status)) status = 0
        key = proxy_normalize_name(name)
        if (len(key) == 0) then
            if (present(status)) status = 2
            return
        end if
        slot = 0
        do i = 1, binary_count
            if (trim(binary_entries(i)%name) == key) then
                slot = i
                exit
            end if
        end do
        if (slot == 0) then
            if (binary_count >= max_registry_entries) then
                if (present(status)) status = 1
                return
            end if
            binary_count = binary_count + 1
            slot = binary_count
        end if
        binary_entries(slot)%name = key
        binary_entries(slot)%is_distance = is_distance
        binary_entries(slot)%conversion = proxy_convert_default
        if (present(conversion)) binary_entries(slot)%conversion = conversion
        binary_entries(slot)%fun => function_pointer
    end subroutine proxy_register_binary

    subroutine proxy_lookup_numeric(name, function_pointer, found, is_distance, conversion)
        character(len=*), intent(in) :: name !! Numeric registry name or alias to locate after normalization.
        procedure(proxy_numeric_callback), pointer, intent(out) :: function_pointer !! Callback pointer on success; null on a
        !! missing method.
        logical, intent(out) :: found !! True when a matching registered numeric method was found.
        logical, intent(out) :: is_distance !! Distance/similarity classification stored with the found method; false when not
        !! found.
        integer, intent(out) :: conversion !! Stored conversion code; defaults to standard conversion when not found.
        character(len=:), allocatable :: key
        integer :: i

        nullify(function_pointer)
        found = .false.
        is_distance = .false.
        conversion = proxy_convert_default
        key = proxy_normalize_name(name)
        do i = 1, numeric_count
            if (trim(numeric_entries(i)%name) == key) then
                function_pointer => numeric_entries(i)%fun
                found = associated(function_pointer)
                is_distance = numeric_entries(i)%is_distance
                conversion = numeric_entries(i)%conversion
                return
            end if
        end do
    end subroutine proxy_lookup_numeric

    subroutine proxy_lookup_binary(name, function_pointer, found, is_distance, conversion)
        character(len=*), intent(in) :: name !! Binary registry name or alias to locate after normalization.
        procedure(proxy_binary_callback), pointer, intent(out) :: function_pointer !! Binary callback pointer on success; null
        !! when no entry matches.
        logical, intent(out) :: found !! True when the requested registered binary method exists.
        logical, intent(out) :: is_distance !! Stored classification of the registered method as distance or similarity.
        integer, intent(out) :: conversion !! Stored conversion code, or standard conversion when no method is found.
        character(len=:), allocatable :: key
        integer :: i

        nullify(function_pointer)
        found = .false.
        is_distance = .false.
        conversion = proxy_convert_default
        key = proxy_normalize_name(name)
        do i = 1, binary_count
            if (trim(binary_entries(i)%name) == key) then
                function_pointer => binary_entries(i)%fun
                found = associated(function_pointer)
                is_distance = binary_entries(i)%is_distance
                conversion = binary_entries(i)%conversion
                return
            end if
        end do
    end subroutine proxy_lookup_binary

    subroutine proxy_registry_clear()
        integer :: i

        do i = 1, numeric_count
            numeric_entries(i)%name = ''
            nullify(numeric_entries(i)%fun)
        end do
        do i = 1, binary_count
            binary_entries(i)%name = ''
            nullify(binary_entries(i)%fun)
        end do
        numeric_count = 0
        binary_count = 0
    end subroutine proxy_registry_clear

    subroutine proxy_registry_numeric_names(names)
        character(len=registry_name_length), allocatable, intent(out) :: names(:) !! Allocated normalized names of all currently
        !! registered custom numeric methods.
        integer :: i

        allocate(names(numeric_count))
        do i = 1, numeric_count
            names(i) = numeric_entries(i)%name
        end do
    end subroutine proxy_registry_numeric_names

    subroutine proxy_registry_binary_names(names)
        character(len=registry_name_length), allocatable, intent(out) :: names(:) !! Allocated normalized names of all currently
        !! registered custom binary methods.
        integer :: i

        allocate(names(binary_count))
        do i = 1, binary_count
            names(i) = binary_entries(i)%name
        end do
    end subroutine proxy_registry_binary_names

end module proxy_registry
