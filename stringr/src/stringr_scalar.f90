! SPDX-License-Identifier: MIT
! SPDX-FileComment: Scalar fixed-string operations for the Fortran stringr translation.
module stringr_scalar
   !! Implements elemental ASCII and byte-oriented string transformations.
   use stringr_types, only: location_type
   implicit none
   private

   public :: str_count, str_detect, str_dup, str_ends, str_equal, str_escape
   public :: str_extract, str_length, str_like, str_locate, str_pad, str_remove
   public :: str_remove_all, str_replace, str_replace_all, str_squish, str_starts
   public :: str_sub, str_to_camel, str_to_kebab, str_to_lower, str_to_sentence
   public :: str_to_snake, str_to_title, str_to_upper, str_trim, str_trunc, str_width, word

contains

   pure elemental integer function str_length(text) result(length)
      !! Returns the trimmed byte length of a string.
      character(len=*), intent(in) :: text !! Input text.
      length = len_trim(text)
   end function str_length

   pure elemental integer function str_width(text) result(width)
      !! Returns display width under the portable one-byte-one-column model.
      character(len=*), intent(in) :: text !! Input text.
      width = len_trim(text)
   end function str_width

   pure elemental logical function str_detect(text, pattern) result(found)
      !! Detects a literal fixed substring.
      character(len=*), intent(in) :: text    !! Input text.
      character(len=*), intent(in) :: pattern !! Literal pattern.
      found = index(text, pattern) > 0
   end function str_detect

   pure elemental logical function str_starts(text, pattern) result(found)
      !! Tests whether text starts with a literal pattern.
      character(len=*), intent(in) :: text    !! Input text.
      character(len=*), intent(in) :: pattern !! Literal prefix.
      found = len(pattern) <= len_trim(text)
      if (found) found = text(:len(pattern)) == pattern
   end function str_starts

   pure elemental logical function str_ends(text, pattern) result(found)
      !! Tests whether trimmed text ends with a literal pattern.
      character(len=*), intent(in) :: text    !! Input text.
      character(len=*), intent(in) :: pattern !! Literal suffix.
      integer :: n
      n = len_trim(text)
      found = len(pattern) <= n
      if (found) found = text(n - len(pattern) + 1:n) == pattern
   end function str_ends

   pure elemental logical function str_equal(left, right, ignore_case) result(equal)
      !! Compares strings with optional portable ASCII case folding.
      character(len=*), intent(in) :: left  !! First text.
      character(len=*), intent(in) :: right !! Second text.
      logical, intent(in), optional :: ignore_case !! Ignore ASCII case when true.
      logical :: folded
      folded = .false.
      if (present(ignore_case)) folded = ignore_case
      if (folded) then
         equal = trim(str_to_lower(left)) == trim(str_to_lower(right))
      else
         equal = trim(left) == trim(right)
      end if
   end function str_equal

   pure elemental integer function str_count(text, pattern) result(count)
      !! Counts nonoverlapping occurrences of a literal pattern.
      character(len=*), intent(in) :: text    !! Input text.
      character(len=*), intent(in) :: pattern !! Nonempty literal pattern.
      integer :: offset, position
      count = 0
      if (len(pattern) == 0) return
      offset = 1
      do while (offset <= len_trim(text))
         position = index(text(offset:len_trim(text)), pattern)
         if (position == 0) exit
         count = count + 1
         offset = offset + position - 1 + len(pattern)
      end do
   end function str_count

   pure elemental function str_locate(text, pattern) result(location)
      !! Returns inclusive bounds of the first literal match.
      character(len=*), intent(in) :: text    !! Input text.
      character(len=*), intent(in) :: pattern !! Literal pattern.
      type(location_type) :: location
      location%start = index(text, pattern)
      if (location%start > 0) location%end = location%start + len(pattern) - 1
   end function str_locate

   pure elemental function str_extract(text, pattern) result(extracted)
      !! Extracts the first literal match or a blank value when absent.
      character(len=*), intent(in) :: text    !! Input text.
      character(len=*), intent(in) :: pattern !! Literal pattern.
      character(len=len(pattern)) :: extracted
      integer :: position
      extracted = ''
      position = index(text, pattern)
      if (position > 0) extracted = text(position:position + len(pattern) - 1)
   end function str_extract

   pure elemental function str_sub(text, start, finish) result(substring)
      !! Extracts inclusive bounds, with negative positions counted from the end.
      character(len=*), intent(in) :: text !! Input text.
      integer, intent(in) :: start         !! First position.
      integer, intent(in), optional :: finish !! Final position; defaults to end.
      character(len=len(text)) :: substring
      integer :: first, last, n
      n = len_trim(text)
      first = merge(start, n + start + 1, start >= 0)
      last = n
      if (present(finish)) last = merge(finish, n + finish + 1, finish >= 0)
      first = max(1, first)
      last = min(n, last)
      substring = ''
      if (first <= last) substring = text(first:last)
   end function str_sub

   pure elemental function str_to_lower(text) result(output)
      !! Converts ASCII uppercase letters to lowercase.
      character(len=*), intent(in) :: text !! Input text.
      character(len=len(text)) :: output
      integer :: code, i
      output = text
      do i = 1, len(text)
         code = iachar(output(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) output(i:i) = achar(code + 32)
      end do
   end function str_to_lower

   pure elemental function str_to_upper(text) result(output)
      !! Converts ASCII lowercase letters to uppercase.
      character(len=*), intent(in) :: text !! Input text.
      character(len=len(text)) :: output
      integer :: code, i
      output = text
      do i = 1, len(text)
         code = iachar(output(i:i))
         if (code >= iachar('a') .and. code <= iachar('z')) output(i:i) = achar(code - 32)
      end do
   end function str_to_upper

   pure elemental function str_to_title(text) result(output)
      !! Converts the first ASCII letter of each whitespace-delimited word to uppercase.
      character(len=*), intent(in) :: text !! Input text.
      character(len=len(text)) :: output
      logical :: new_word
      integer :: i
      output = str_to_lower(text)
      new_word = .true.
      do i = 1, len_trim(output)
         if (is_space(output(i:i))) then
            new_word = .true.
         else if (new_word) then
            output(i:i) = str_to_upper(output(i:i))
            new_word = .false.
         end if
      end do
   end function str_to_title

   pure elemental function str_to_sentence(text) result(output)
      !! Converts text to lowercase and capitalizes its first ASCII character.
      character(len=*), intent(in) :: text !! Input text.
      character(len=len(text)) :: output
      integer :: i
      output = str_to_lower(text)
      do i = 1, len_trim(output)
         if (.not. is_space(output(i:i))) then
            output(i:i) = str_to_upper(output(i:i))
            exit
         end if
      end do
   end function str_to_sentence

   pure elemental function str_trim(text, side) result(output)
      !! Removes ASCII whitespace from both sides or one selected side.
      character(len=*), intent(in) :: text !! Input text.
      character(len=*), intent(in), optional :: side !! `both`, `left`, or `right`.
      character(len=len(text)) :: output
      character(len=5) :: selected
      integer :: first, last
      selected = 'both'
      if (present(side)) selected = side
      first = 1
      last = len_trim_whitespace(text)
      if (selected == 'both' .or. selected == 'left') then
         do while (first <= last .and. is_space(text(first:first)))
            first = first + 1
         end do
      end if
      if (selected == 'left') last = len(text)
      if (selected == 'right') first = 1
      output = ''
      if (first <= last) output = text(first:last)
   end function str_trim

   pure elemental function str_squish(text) result(output)
      !! Trims text and collapses ASCII whitespace runs to one space.
      character(len=*), intent(in) :: text !! Input text.
      character(len=len(text)) :: output
      logical :: pending_space
      integer :: i, n
      output = ''
      n = 0
      pending_space = .false.
      do i = 1, len_trim_whitespace(text)
         if (is_space(text(i:i))) then
            if (n > 0) pending_space = .true.
         else
            if (pending_space) then
               n = n + 1
               output(n:n) = ' '
            end if
            n = n + 1
            output(n:n) = text(i:i)
            pending_space = .false.
         end if
      end do
   end function str_squish

   pure elemental function str_dup(text, times) result(output)
      !! Repeats text a nonnegative number of times.
      character(len=*), intent(in) :: text !! Input text.
      integer, intent(in) :: times         !! Repetition count.
      character(len=max(0, len_trim(text) * times)) :: output
      integer :: i, n
      output = ''
      if (times < 0) return
      n = len_trim(text)
      do i = 1, times
         output((i - 1) * n + 1:i * n) = text(:n)
      end do
   end function str_dup

   pure elemental function str_pad(text, width, side, pad) result(output)
      !! Pads text to a minimum byte width.
      character(len=*), intent(in) :: text !! Input text.
      integer, intent(in) :: width         !! Minimum output width.
      character(len=*), intent(in), optional :: side !! `left`, `right`, or `both`.
      character(len=1), intent(in), optional :: pad !! Padding byte; space by default.
      character(len=max(len_trim(text), width)) :: output
      character(len=5) :: selected
      character(len=1) :: fill
      integer :: left_count, n
      selected = 'left'
      if (present(side)) selected = side
      fill = ' '
      if (present(pad)) fill = pad
      n = len_trim(text)
      left_count = max(0, width - n)
      if (selected == 'right') left_count = 0
      if (selected == 'both') left_count = max(0, width - n) / 2
      output = repeat(fill, len(output))
      if (n > 0) output(left_count + 1:left_count + n) = text(:n)
   end function str_pad

   pure elemental function str_trunc(text, width, side, ellipsis) result(output)
      !! Truncates text to a requested width and inserts an ellipsis marker.
      character(len=*), intent(in) :: text !! Input text.
      integer, intent(in) :: width         !! Maximum output width.
      character(len=*), intent(in), optional :: side !! `right` or `left`.
      character(len=*), intent(in), optional :: ellipsis !! Marker, default `...`.
      character(len=max(0, width)) :: output
      character(len=:), allocatable :: marker
      character(len=5) :: selected
      integer :: keep, n
      marker = '...'
      if (present(ellipsis)) marker = ellipsis
      selected = 'right'
      if (present(side)) selected = side
      output = ''
      n = len_trim(text)
      if (n <= width) then
         if (n > 0) output(:n) = text(:n)
         return
      end if
      keep = max(0, width - len(marker))
      if (selected == 'left') then
         output = marker // text(n - keep + 1:n)
      else
         output = text(:keep) // marker
      end if
   end function str_trunc

   pure function str_replace(text, pattern, replacement) result(output)
      !! Replaces the first literal occurrence.
      character(len=*), intent(in) :: text        !! Input text.
      character(len=*), intent(in) :: pattern     !! Literal pattern.
      character(len=*), intent(in) :: replacement !! Replacement text.
      character(len=:), allocatable :: output
      integer :: position
      position = index(text(:len_trim(text)), pattern)
      if (position == 0 .or. len(pattern) == 0) then
         output = text(:len_trim(text))
      else
         output = text(:position - 1) // replacement // &
                  text(position + len(pattern):len_trim(text))
      end if
   end function str_replace

   pure function str_replace_all(text, pattern, replacement) result(output)
      !! Replaces every nonoverlapping literal occurrence.
      character(len=*), intent(in) :: text        !! Input text.
      character(len=*), intent(in) :: pattern     !! Literal pattern.
      character(len=*), intent(in) :: replacement !! Replacement text.
      character(len=:), allocatable :: output
      integer :: position, start
      output = ''
      if (len(pattern) == 0) then
         output = text(:len_trim(text))
         return
      end if
      start = 1
      do
         position = index(text(start:len_trim(text)), pattern)
         if (position == 0) then
            output = output // text(start:len_trim(text))
            exit
         end if
         position = start + position - 1
         output = output // text(start:position - 1) // replacement
         start = position + len(pattern)
         if (start > len_trim(text)) exit
      end do
   end function str_replace_all

   pure function str_remove(text, pattern) result(output)
      !! Removes the first literal occurrence.
      character(len=*), intent(in) :: text    !! Input text.
      character(len=*), intent(in) :: pattern !! Literal pattern.
      character(len=:), allocatable :: output
      output = str_replace(text, pattern, '')
   end function str_remove

   pure function str_remove_all(text, pattern) result(output)
      !! Removes every nonoverlapping literal occurrence.
      character(len=*), intent(in) :: text    !! Input text.
      character(len=*), intent(in) :: pattern !! Literal pattern.
      character(len=:), allocatable :: output
      output = str_replace_all(text, pattern, '')
   end function str_remove_all

   pure elemental logical function str_like(text, pattern, ignore_case) result(matches)
      !! Matches SQL-style `%` and `_` wildcards using portable fixed bytes.
      character(len=*), intent(in) :: text    !! Input text.
      character(len=*), intent(in) :: pattern !! SQL-like pattern.
      logical, intent(in), optional :: ignore_case !! Ignore ASCII case when true.
      character(len=len(text)) :: subject
      character(len=len(pattern)) :: model
      logical :: folded
      folded = .false.
      if (present(ignore_case)) folded = ignore_case
      subject = text
      model = pattern
      if (folded) then
         subject = str_to_lower(text)
         model = str_to_lower(pattern)
      end if
      matches = wildcard_match(trim(subject), trim(model), 1, 1)
   end function str_like

   pure elemental function str_escape(text) result(output)
      !! Prefixes common regular-expression metacharacters with a backslash.
      character(len=*), intent(in) :: text !! Input text.
      character(len=2 * len(text)) :: output
      integer :: i, n
      output = ''
      n = 0
      do i = 1, len_trim(text)
         if (index('.^$|?*+()[]{}\', text(i:i)) > 0) then
            n = n + 1
            output(n:n) = '\'
         end if
         n = n + 1
         output(n:n) = text(i:i)
      end do
   end function str_escape

   pure elemental function str_to_snake(text) result(output)
      !! Converts ASCII words to lower snake case.
      character(len=*), intent(in) :: text !! Input text.
      character(len=len(text)) :: output
      output = separated_case(text, '_')
   end function str_to_snake

   pure elemental function str_to_kebab(text) result(output)
      !! Converts ASCII words to lower kebab case.
      character(len=*), intent(in) :: text !! Input text.
      character(len=len(text)) :: output
      output = separated_case(text, '-')
   end function str_to_kebab

   pure elemental function str_to_camel(text) result(output)
      !! Converts ASCII words to lower camel case.
      character(len=*), intent(in) :: text !! Input text.
      character(len=len(text)) :: output
      logical :: capitalize
      integer :: i, n
      output = ''
      capitalize = .false.
      n = 0
      do i = 1, len_trim(text)
         if (is_separator(text(i:i))) then
            capitalize = n > 0
         else
            n = n + 1
            output(n:n) = str_to_lower(text(i:i))
            if (capitalize) output(n:n) = str_to_upper(output(n:n))
            capitalize = .false.
         end if
      end do
   end function str_to_camel

   pure elemental function word(text, start, finish, separator) result(output)
      !! Extracts an inclusive range of whitespace- or delimiter-separated words.
      character(len=*), intent(in) :: text !! Input text.
      integer, intent(in), optional :: start !! First word, default one.
      integer, intent(in), optional :: finish !! Last word, default first.
      character(len=*), intent(in), optional :: separator !! Literal separator.
      character(len=len(text)) :: output
      character(len=:), allocatable :: delimiter
      integer :: first_word, last_word, i, first, last, current
      delimiter = ' '
      if (present(separator)) delimiter = separator
      first_word = 1
      if (present(start)) first_word = start
      last_word = first_word
      if (present(finish)) last_word = finish
      output = ''
      first = 1
      current = 1
      do while (first <= len_trim(text))
         i = index(text(first:len_trim(text)), delimiter)
         if (i == 0) then
            last = len_trim(text)
         else
            last = first + i - 2
         end if
         if (current >= first_word .and. current <= last_word) then
            if (len_trim(output) > 0) output = trim(output) // delimiter
            output = trim(output) // text(first:last)
         end if
         if (i == 0 .or. current >= last_word) exit
         first = last + 1 + len(delimiter)
         current = current + 1
      end do
   end function word

   pure recursive logical function wildcard_match(text, pattern, i, j) result(matches)
      !! Recursively matches one SQL-like wildcard pattern.
      character(len=*), intent(in) :: text    !! Input text.
      character(len=*), intent(in) :: pattern !! Pattern text.
      integer, intent(in) :: i                !! Current text position.
      integer, intent(in) :: j                !! Current pattern position.
      if (j > len(pattern)) then
         matches = i > len(text)
      else if (pattern(j:j) == '%') then
         matches = wildcard_match(text, pattern, i, j + 1)
         if (.not. matches .and. i <= len(text)) matches = wildcard_match(text, pattern, i + 1, j)
      else if (i > len(text)) then
         matches = .false.
      else if (pattern(j:j) == '_' .or. pattern(j:j) == text(i:i)) then
         matches = wildcard_match(text, pattern, i + 1, j + 1)
      else
         matches = .false.
      end if
   end function wildcard_match

   pure elemental function separated_case(text, separator) result(output)
      !! Converts ASCII words to a lower separated case.
      character(len=*), intent(in) :: text !! Input text.
      character(len=1), intent(in) :: separator !! Output separator.
      character(len=len(text)) :: output
      logical :: pending
      integer :: i, n
      output = ''
      pending = .false.
      n = 0
      do i = 1, len_trim(text)
         if (is_separator(text(i:i))) then
            pending = n > 0
         else
            if (pending .and. n < len(output)) then
               n = n + 1
               output(n:n) = separator
            end if
            n = n + 1
            output(n:n) = str_to_lower(text(i:i))
            pending = .false.
         end if
      end do
   end function separated_case

   pure elemental logical function is_space(character) result(space)
      !! Recognizes portable ASCII whitespace bytes.
      character(len=1), intent(in) :: character !! Character to classify.
      space = index(' ' // achar(9) // achar(10) // achar(13), character) > 0
   end function is_space

   pure elemental logical function is_separator(character) result(separator)
      !! Recognizes common ASCII word separators.
      character(len=1), intent(in) :: character !! Character to classify.
      separator = is_space(character) .or. index('_-.', character) > 0
   end function is_separator

   pure elemental integer function len_trim_whitespace(text) result(last)
      !! Finds the final non-whitespace byte.
      character(len=*), intent(in) :: text !! Input text.
      last = len(text)
      do while (last > 0 .and. is_space(text(last:last)))
         last = last - 1
      end do
   end function len_trim_whitespace

end module stringr_scalar
