module proxy_utils
    use proxy_kinds, only: dp
    implicit none
    private

    integer, parameter, public :: proxy_convert_default = 1
    integer, parameter, public :: proxy_convert_one_minus = 2
    integer, parameter, public :: proxy_convert_none = 3

    public :: proxy_normalize_name, proxy_simil_to_dist, proxy_dist_to_simil

contains

    pure function proxy_normalize_name(text) result(key)
        character(len=*), intent(in) :: text !! User-facing method name; ASCII letters are lowercased and nonalphanumeric
        !! separators are removed for alias matching.
        character(len=:), allocatable :: key
        character(len=len_trim(text)) :: work
        integer :: c
        integer :: i
        integer :: j

        work = ''
        j = 0
        do i = 1, len_trim(text)
            c = iachar(text(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) c = c + 32
            if ((c >= iachar('a') .and. c <= iachar('z')) .or. (c >= iachar('0') .and. c <= iachar('9'))) then
                j = j + 1
                work(j:j) = achar(c)
            end if
        end do
        key = work(:j)
    end function proxy_normalize_name

    pure function proxy_simil_to_dist(value, conversion) result(converted)
        real(dp), intent(in) :: value !! Similarity value to convert to a dissimilarity/distance representation.
        integer, intent(in) :: conversion !! Conversion code: default uses `1-abs(s)`, one-minus uses `1-s`, and none leaves the
        !! value unchanged.
        real(dp) :: converted

        select case (conversion)
        case (proxy_convert_default)
            converted = 1.0_dp - abs(value)
        case (proxy_convert_one_minus)
            converted = 1.0_dp - value
        case default
            converted = value
        end select
    end function proxy_simil_to_dist

    pure function proxy_dist_to_simil(value, conversion) result(converted)
        real(dp), intent(in) :: value !! Distance/dissimilarity value to convert to a similarity representation.
        integer, intent(in) :: conversion !! Conversion code: default uses `1/(1+d)`, one-minus uses `1-d`, and none leaves the
        !! value unchanged.
        real(dp) :: converted

        select case (conversion)
        case (proxy_convert_default)
            converted = 1.0_dp / (1.0_dp + value)
        case (proxy_convert_one_minus)
            converted = 1.0_dp - value
        case default
            converted = value
        end select
    end function proxy_dist_to_simil

end module proxy_utils
