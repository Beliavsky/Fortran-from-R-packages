! SPDX-License-Identifier: MIT
! SPDX-FileComment: Name selectors for the Fortran tidyselect translation.
module tidyselect
   !! Implements deterministic selectors returning one-based name positions.
   implicit none
   private
   public :: all_of, any_of, contains, ends_with, eval_relocate, eval_rename
   public :: eval_select, everything, last_col, matches, num_range, one_of
   public :: starts_with, vars_pull, vars_rename, vars_select

contains
   pure function everything(names) result(indices)
      !! Selects every name.
      character(len=*), intent(in) :: names(:) !! Available names.
      integer, allocatable :: indices(:)
      integer :: i
      indices = [(i, i=1,size(names))]
   end function everything

   pure function all_of(names, requested) result(indices)
      !! Resolves requested names and rejects absent names.
      character(len=*), intent(in) :: names(:) !! Available names.
      character(len=*), intent(in) :: requested(:) !! Required names.
      integer, allocatable :: indices(:)
      integer :: i
      allocate(indices(size(requested)))
      do i=1,size(requested)
         indices(i)=find_name(names,requested(i))
         if(indices(i)==0) error stop 'all_of: requested name is absent'
      end do
   end function all_of

   pure function any_of(names, requested) result(indices)
      !! Resolves requested names while ignoring absent names.
      character(len=*), intent(in) :: names(:) !! Available names.
      character(len=*), intent(in) :: requested(:) !! Candidate names.
      integer, allocatable :: indices(:), found(:)
      integer :: i
      allocate(found(size(requested)))
      do i=1,size(requested)
         found(i)=find_name(names,requested(i))
      end do
      indices=pack(found,found>0)
   end function any_of

   pure function one_of(names, requested) result(indices)
      !! Compatibility alias for `any_of`.
      character(len=*), intent(in) :: names(:) !! Available names.
      character(len=*), intent(in) :: requested(:) !! Candidate names.
      integer, allocatable :: indices(:)
      indices=any_of(names,requested)
   end function one_of

   pure function starts_with(names, prefix, ignore_case) result(indices)
      !! Selects names beginning with a fixed prefix.
      character(len=*), intent(in) :: names(:) !! Available names.
      character(len=*), intent(in) :: prefix !! Prefix.
      logical, intent(in), optional :: ignore_case !! Ignore ASCII case.
      integer, allocatable :: indices(:)
      indices=select_pattern(names,prefix,1,ignore_case)
   end function starts_with

   pure function ends_with(names, suffix, ignore_case) result(indices)
      !! Selects names ending with a fixed suffix.
      character(len=*), intent(in) :: names(:) !! Available names.
      character(len=*), intent(in) :: suffix !! Suffix.
      logical, intent(in), optional :: ignore_case !! Ignore ASCII case.
      integer, allocatable :: indices(:)
      indices=select_pattern(names,suffix,2,ignore_case)
   end function ends_with

   pure function contains(names, pattern, ignore_case) result(indices)
      !! Selects names containing a fixed substring.
      character(len=*), intent(in) :: names(:) !! Available names.
      character(len=*), intent(in) :: pattern !! Substring.
      logical, intent(in), optional :: ignore_case !! Ignore ASCII case.
      integer, allocatable :: indices(:)
      indices=select_pattern(names,pattern,3,ignore_case)
   end function contains

   pure function matches(names, pattern, ignore_case) result(indices)
      !! Selects fixed substring matches; regular expressions are not interpreted.
      character(len=*), intent(in) :: names(:) !! Available names.
      character(len=*), intent(in) :: pattern !! Fixed substring.
      logical, intent(in), optional :: ignore_case !! Ignore ASCII case.
      integer, allocatable :: indices(:)
      indices=contains(names,pattern,ignore_case)
   end function matches

   pure function last_col(names, offset) result(indices)
      !! Selects one column by offset from the final column.
      character(len=*), intent(in) :: names(:) !! Available names.
      integer, intent(in), optional :: offset !! Nonnegative offset.
      integer, allocatable :: indices(:)
      integer :: skip
      skip=0
      if(present(offset)) skip=offset
      if(skip<0 .or. skip>=size(names)) error stop 'last_col: invalid offset'
      indices=[size(names)-skip]
   end function last_col

   pure function num_range(names, prefix, range, width) result(indices)
      !! Selects names formed from a prefix and integer range.
      character(len=*), intent(in) :: names(:) !! Available names.
      character(len=*), intent(in) :: prefix !! Name prefix.
      integer, intent(in) :: range(:) !! Integer suffixes.
      integer, intent(in), optional :: width !! Zero-padding width.
      integer, allocatable :: indices(:)
      character(len=64), allocatable :: requested(:)
      character(len=32) :: number
      integer :: i,w
      w=0
      if(present(width)) w=width
      allocate(requested(size(range)))
      do i=1,size(range)
         if(w>0) then
            write(number,'(i0)') range(i)
            requested(i)=prefix//repeat('0',max(0,w-len_trim(number)))//trim(number)
         else
            write(requested(i),'(a,i0)') trim(prefix),range(i)
         end if
      end do
      indices=any_of(names,requested)
   end function num_range

   pure function eval_select(names, selections) result(indices)
      !! Combines signed selections, with negative positions excluding names.
      character(len=*), intent(in) :: names(:) !! Available names.
      integer, intent(in) :: selections(:) !! Signed one-based positions.
      integer, allocatable :: indices(:)
      logical, allocatable :: keep(:)
      integer :: i
      allocate(keep(size(names)),source=.false.)
      if(any(selections<0)) keep=.true.
      do i=1,size(selections)
         if(abs(selections(i))<1 .or. abs(selections(i))>size(names)) error stop 'eval_select: index out of range'
         keep(abs(selections(i)))=selections(i)>0
      end do
      indices=pack([(i,i=1,size(names))],keep)
   end function eval_select

   pure function vars_select(names, requested) result(indices)
      !! Resolves exact required names.
      character(len=*), intent(in) :: names(:) !! Available names.
      character(len=*), intent(in) :: requested(:) !! Required names.
      integer, allocatable :: indices(:)
      indices=all_of(names,requested)
   end function vars_select

   pure integer function vars_pull(names, requested) result(index_value)
      !! Resolves exactly one required name.
      character(len=*), intent(in) :: names(:) !! Available names.
      character(len=*), intent(in) :: requested !! Required name.
      index_value=find_name(names,requested)
      if(index_value==0) error stop 'vars_pull: name is absent'
   end function vars_pull

   pure function vars_rename(names, positions, replacements) result(output)
      !! Renames selected positions and rejects duplicate output names.
      character(len=*), intent(in) :: names(:) !! Existing names.
      integer, intent(in) :: positions(:) !! Positions to rename.
      character(len=*), intent(in) :: replacements(:) !! Replacement names.
      character(len=:), allocatable :: output(:)
      integer :: i,j
      if(size(positions)/=size(replacements)) error stop 'vars_rename: size mismatch'
      allocate(character(len=len(names)) :: output(size(names)))
      output=names
      do i=1,size(positions)
         if(positions(i)<1 .or. positions(i)>size(names)) error stop 'vars_rename: index out of range'
         output(positions(i))=replacements(i)
      end do
      do i=1,size(output)
         do j=1,i-1
            if(output(i)==output(j)) error stop 'vars_rename: duplicate output name'
         end do
      end do
   end function vars_rename

   pure function eval_rename(names, positions, replacements) result(output)
      !! Renames selected positions.
      character(len=*), intent(in) :: names(:) !! Existing names.
      integer, intent(in) :: positions(:) !! Positions.
      character(len=*), intent(in) :: replacements(:) !! New names.
      character(len=:), allocatable :: output(:)
      allocate(character(len=len(names)) :: output(size(names)))
      output=vars_rename(names,positions,replacements)
   end function eval_rename

   pure function eval_relocate(names, selected, before) result(indices)
      !! Returns a permutation moving selected positions before an anchor.
      character(len=*), intent(in) :: names(:) !! Available names.
      integer, intent(in) :: selected(:) !! Positions to move.
      integer, intent(in), optional :: before !! Anchor; one by default.
      integer, allocatable :: indices(:)
      logical, allocatable :: moved(:)
      integer :: anchor,i,j,k
      allocate(indices(size(names)))
      allocate(moved(size(names)),source=.false.)
      moved(selected)=.true.
      anchor=1
      if(present(before)) anchor=before
      k=0
      do j=1,size(names)
         if(j==anchor) then
            do i=1,size(selected)
            k=k+1
            indices(k)=selected(i)
            end do
         end if
         if(.not.moved(j)) then
         k=k+1
         indices(k)=j
         end if
      end do
   end function eval_relocate

   pure function select_pattern(names,pattern,mode,ignore_case) result(indices)
      !! Implements prefix, suffix, and containment selection modes.
      character(len=*),intent(in)::names(:),pattern
      integer,intent(in)::mode
      logical,intent(in),optional::ignore_case
      integer,allocatable::indices(:)
      logical,allocatable::keep(:)
      character(len=len(names))::candidate
      character(len=len(pattern))::target
      logical::fold
      integer::i,n
      fold=.false.
      if(present(ignore_case)) fold=ignore_case
      target=pattern
      if(fold) target=lower(pattern)
      allocate(keep(size(names)))
      do i=1,size(names)
         candidate=names(i)
         if(fold) candidate=lower(names(i))
         n=len_trim(candidate)
         select case(mode)
         case(1)
            keep(i)=len_trim(target)<=n
            if(keep(i)) keep(i)=candidate(:len_trim(target))==trim(target)
         case(2)
            keep(i)=len_trim(target)<=n
            if(keep(i)) keep(i)=candidate(n-len_trim(target)+1:n)==trim(target)
         case default
         keep(i)=index(candidate,trim(target))>0
         end select
      end do
      indices=pack([(i,i=1,size(names))],keep)
   end function select_pattern

   pure elemental function lower(text) result(output)
      !! Converts ASCII uppercase bytes to lowercase.
      character(len=*),intent(in)::text
      character(len=len(text))::output
      integer::i,c
      output=text
      do i=1,len(text)
      c=iachar(output(i:i))
      if(c>=65.and.c<=90)output(i:i)=achar(c+32)
      end do
   end function lower

   pure integer function find_name(names,name) result(position)
      !! Finds an exact name and returns zero when absent.
      character(len=*),intent(in)::names(:),name
      integer::i
      position=0
      do i=1,size(names)
      if(names(i)==name)then
      position=i
      return
      end if
      end do
   end function find_name
end module tidyselect
