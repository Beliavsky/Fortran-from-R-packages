module proxy_binary_measures
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
    use proxy_kinds, only: dp
    use proxy_ieee, only: proxy_nan, proxy_is_missing
    implicit none
    private

    public :: binary_counts, binary_counts_numeric, binary_similarity_from_counts
    public :: binary_jaccard_similarity, binary_matching_distance

contains

    pure subroutine binary_counts(x, y, a, b, c, d, n)
        logical, intent(in) :: x(:) !! First Boolean feature vector; every component is treated as observed because Fortran
        !! LOGICAL has no NA value.
        logical, intent(in) :: y(:) !! Second Boolean feature vector; components beyond the shorter input are ignored.
        integer, intent(out) :: a !! Number of components where both vectors are TRUE.
        integer, intent(out) :: b !! Number of components where `x` is TRUE and `y` is FALSE.
        integer, intent(out) :: c !! Number of components where `x` is FALSE and `y` is TRUE.
        integer, intent(out) :: d !! Number of components where both vectors are FALSE.
        integer, intent(out) :: n !! Number of compared components, equal to `min(size(x),size(y))`.
        integer :: i

        a = 0
        b = 0
        c = 0
        d = 0
        n = min(size(x), size(y))
        do i = 1, n
            if (x(i) .and. y(i)) then
                a = a + 1
            else if (x(i)) then
                b = b + 1
            else if (y(i)) then
                c = c + 1
            else
                d = d + 1
            end if
        end do
    end subroutine binary_counts

    pure subroutine binary_counts_numeric(x, y, a, b, c, d, n)
        real(dp), intent(in) :: x(:) !! First numeric vector, coerced to Boolean using zero=FALSE and nonzero=TRUE; NaN
        !! components are omitted.
        real(dp), intent(in) :: y(:) !! Second numeric vector with matching Boolean coercion and pairwise NaN omission.
        integer, intent(out) :: a !! Number of valid components where both values coerce to TRUE.
        integer, intent(out) :: b !! Number of valid components where only `x` coerces to TRUE.
        integer, intent(out) :: c !! Number of valid components where only `y` coerces to TRUE.
        integer, intent(out) :: d !! Number of valid components where both values coerce to FALSE.
        integer, intent(out) :: n !! Number of valid nonmissing component pairs included in the counts.
        integer :: i
        logical :: tx
        logical :: ty

        a = 0
        b = 0
        c = 0
        d = 0
        n = 0
        do i = 1, min(size(x), size(y))
            if (proxy_is_missing(x(i)) .or. proxy_is_missing(y(i))) cycle
            tx = abs(x(i)) > tiny(1.0_dp)
            ty = abs(y(i)) > tiny(1.0_dp)
            if (tx .and. ty) then
                a = a + 1
            else if (tx) then
                b = b + 1
            else if (ty) then
                c = c + 1
            else
                d = d + 1
            end if
            n = n + 1
        end do
    end subroutine binary_counts_numeric

    pure function binary_jaccard_similarity(x, y) result(similarity)
        logical, intent(in) :: x(:) !! First Boolean feature vector in the Jaccard coefficient.
        logical, intent(in) :: y(:) !! Second Boolean feature vector with the same conceptual feature set as `x`.
        real(dp) :: similarity
        integer :: a
        integer :: b
        integer :: c
        integer :: d
        integer :: n

        call binary_counts(x, y, a, b, c, d, n)
        similarity = binary_similarity_from_counts('jaccard', a, b, c, d, n)
    end function binary_jaccard_similarity

    pure function binary_matching_distance(x, y) result(distance)
        logical, intent(in) :: x(:) !! First Boolean vector for simple Hamming/matching dissimilarity.
        logical, intent(in) :: y(:) !! Second Boolean vector; the result is the proportion of discordant valid components.
        real(dp) :: distance
        integer :: a
        integer :: b
        integer :: c
        integer :: d
        integer :: n

        call binary_counts(x, y, a, b, c, d, n)
        if (n == 0) then
            distance = proxy_nan()
        else
            distance = real(b + c, dp) / real(n, dp)
        end if
    end function binary_matching_distance

    pure function binary_similarity_from_counts(method, a, b, c, d, n) result(similarity)
        character(len=*), intent(in) :: method !! Binary similarity name or alias, matched case-insensitively after normalization
        !! by the caller.
        integer, intent(in) :: a !! Count of TRUE/TRUE pairs in the binary contingency table.
        integer, intent(in) :: b !! Count of TRUE/FALSE pairs in the binary contingency table.
        integer, intent(in) :: c !! Count of FALSE/TRUE pairs in the binary contingency table.
        integer, intent(in) :: d !! Count of FALSE/FALSE pairs in the binary contingency table.
        integer, intent(in) :: n !! Number of valid component pairs; normally `a+b+c+d` and possibly zero when all data are missing.
        real(dp) :: similarity
        real(dp) :: ad
        real(dp) :: bc
        real(dp) :: den
        character(len=:), allocatable :: key

        key = normalize(method)
        if (n == 0) then
            similarity = proxy_nan()
            return
        end if

        select case (key)
        case ('jaccard', 'binary', 'reyssac', 'roux')
            den = real(a + b + c, dp)
            if (a + b + c == 0) then
                similarity = 1.0_dp
            else
                similarity = real(a, dp) / den
            end if
        case ('kulczynski1')
            similarity = ratio(real(a, dp), real(b + c, dp))
        case ('kulczynski2')
            similarity = 0.5_dp * (ratio(real(a, dp), real(a + b, dp)) + &
                                   ratio(real(a, dp), real(a + c, dp)))
        case ('mountford')
            similarity = ratio(2.0_dp * real(a, dp), &
                               real(a * (b + c) + 2 * b * c, dp))
        case ('fager', 'mcgowan', 'fagermcgowan')
            den = sqrt(real((a + b) * (a + c), dp))
            similarity = ratio(real(a, dp), den) - sqrt(real(a + c, dp)) / 2.0_dp
        case ('russel', 'rao', 'russelrao')
            similarity = real(a, dp) / real(n, dp)
        case ('simplematching', 'sokalmichener')
            similarity = real(a + d, dp) / real(n, dp)
        case ('hamman')
            similarity = real(a + d - b - c, dp) / real(n, dp)
        case ('faith')
            similarity = (real(a, dp) + real(d, dp) / 2.0_dp) / real(n, dp)
        case ('tanimoto', 'rogers', 'rogerstanimoto')
            similarity = ratio(real(a + d, dp), real(a + 2 * (b + c) + d, dp))
        case ('dice', 'czekanowski', 'sorensen')
            similarity = ratio(2.0_dp * real(a, dp), real(2 * a + b + c, dp))
        case ('phi')
            den = sqrt(real(a + b, dp)) * sqrt(real(c + d, dp)) * &
                  sqrt(real(a + c, dp)) * sqrt(real(b + d, dp))
            similarity = ratio(real(a * d - b * c, dp), den)
        case ('stiles')
            den = abs(real(a * d - b * c, dp)) - 0.5_dp * real(n, dp)
            if (den < 0.0_dp) then
                similarity = proxy_nan()
            else if (den <= tiny(1.0_dp)) then
                similarity = ieee_value(0.0_dp, ieee_negative_inf)
            else if (a + b <= 0 .or. c + d <= 0 .or. a + c <= 0 .or. b + d <= 0) then
                similarity = proxy_nan()
            else
                similarity = log(real(n, dp)) + 2.0_dp * log(den) - &
                             log(real(a + b, dp)) - log(real(c + d, dp)) - &
                             log(real(a + c, dp)) - log(real(b + d, dp))
            end if
        case ('michael')
            similarity = ratio(4.0_dp * real(a * d - b * c, dp), &
                               real((a + d)**2 + (b + c)**2, dp))
        case ('mozley', 'margalef', 'mozleymargalef')
            similarity = ratio(real(a * n, dp), real((a + b) * (a + c), dp))
        case ('yule')
            ad = real(a * d, dp)
            bc = real(b * c, dp)
            similarity = ratio(ad - bc, ad + bc)
        case ('yule2')
            ad = sqrt(real(a * d, dp))
            bc = sqrt(real(b * c, dp))
            similarity = ratio(ad - bc, ad + bc)
        case ('ochiai')
            similarity = ratio(real(a, dp), sqrt(real((a + b) * (a + c), dp)))
        case ('simpson')
            similarity = ratio(real(a, dp), real(min(a + b, a + c), dp))
        case ('braunblanquet')
            similarity = ratio(real(a, dp), real(max(a + b, a + c), dp))
        case default
            similarity = proxy_nan()
        end select
    end function binary_similarity_from_counts

    pure function ratio(numerator, denominator) result(value)
        real(dp), intent(in) :: numerator !! Numerator of a binary coefficient ratio.
        real(dp), intent(in) :: denominator !! Denominator of the ratio; values near zero return NaN rather than trapping.
        real(dp) :: value

        if (abs(denominator) <= tiny(1.0_dp)) then
            if (abs(numerator) <= tiny(1.0_dp)) then
                value = proxy_nan()
            else if (numerator < 0.0_dp) then
                value = ieee_value(0.0_dp, ieee_negative_inf)
            else
                value = ieee_value(0.0_dp, ieee_positive_inf)
            end if
        else
            value = numerator / denominator
        end if
    end function ratio

    pure function normalize(text) result(key)
        character(len=*), intent(in) :: text !! Method name or alias; punctuation and whitespace are removed and ASCII letters
        !! are lowercased.
        character(len=:), allocatable :: key
        character(len=len_trim(text)) :: work
        integer :: c
        integer :: i
        integer :: j

        work = ''
        j = 0
        do i = 1, len_trim(text)
            c = iachar(text(i:i))
            if ((c >= iachar('A') .and. c <= iachar('Z'))) c = c + 32
            if ((c >= iachar('a') .and. c <= iachar('z')) .or. (c >= iachar('0') .and. c <= iachar('9'))) then
                j = j + 1
                work(j:j) = achar(c)
            end if
        end do
        key = work(:j)
    end function normalize

end module proxy_binary_measures
