#!/usr/bin/env python3
"""Tests for conservative Fortran syntax modernization."""

from __future__ import annotations

import unittest

from modernize_fortran_syntax import transform_text


class ModernizeFortranSyntaxTests(unittest.TestCase):
    def test_removes_module_procedure_implicit_none(self) -> None:
        source = """module example
  implicit none
contains
  subroutine work(x)
    implicit none
    real, intent(inout) :: x
    x = x + 1.0
  end subroutine work
end module example
"""
        transformed, implicit_count, _, _, _ = transform_text(source)
        self.assertEqual(implicit_count, 1)
        self.assertEqual(transformed.count("implicit none"), 1)

    def test_keeps_standalone_and_interface_implicit_none(self) -> None:
        source = """subroutine standalone(x)
  implicit none
  real, intent(inout) :: x
end subroutine standalone

module interfaces
  implicit none
  interface
    subroutine callback(x)
      implicit none
      real, intent(inout) :: x
    end subroutine callback
  end interface
end module interfaces
"""
        transformed, implicit_count, _, _, _ = transform_text(source)
        self.assertEqual(implicit_count, 0)
        self.assertEqual(transformed, source)

    def test_removes_unreferenced_do_name(self) -> None:
        source = """loop_10: do i = 1, n
  x(i) = i
end do loop_10
"""
        transformed, _, do_count, _, _ = transform_text(source)
        self.assertEqual(do_count, 1)
        self.assertEqual(transformed, "do i = 1, n\n  x(i) = i\nend do\n")

    def test_keeps_do_name_targeted_by_cycle(self) -> None:
        source = """search: do i = 1, n
  if (x(i) < 0.0) cycle search
end do search
"""
        transformed, _, do_count, _, _ = transform_text(source)
        self.assertEqual(do_count, 0)
        self.assertEqual(transformed, source)

    def test_keeps_legacy_labelled_do(self) -> None:
        source = """do 100 i = 1, n
  x(i) = i
100 continue
"""
        transformed, _, do_count, _, _ = transform_text(source)
        self.assertEqual(do_count, 0)
        self.assertEqual(transformed, source)

    def test_ignores_names_in_comments_and_strings(self) -> None:
        source = """unused: do i = 1, n
  print *, "cycle unused" ! exit unused
end do unused
"""
        transformed, _, do_count, _, _ = transform_text(source)
        self.assertEqual(do_count, 1)
        self.assertNotIn("unused:", transformed)

    def test_splits_an_overlong_old_style_declaration_without_continuations(self) -> None:
        source = "integer first_argument, second_argument, third_argument, fourth_argument\n"
        transformed, _, _, _, declaration_count = transform_text(source, max_line_length=45)
        self.assertEqual(declaration_count, 1)
        self.assertGreater(len(transformed.splitlines()), 1)
        self.assertTrue(all("::" in line for line in transformed.splitlines()))
        self.assertNotIn("&", transformed)

    def test_adds_explicit_save_to_initialized_procedure_local(self) -> None:
        source = """module example
  implicit none
contains
  subroutine work(x)
    real, intent(out) :: x
    real :: initial_value = 1.0
    x = initial_value
  end subroutine work
end module example
"""
        transformed, _, _, save_count, _ = transform_text(source)
        self.assertEqual(save_count, 1)
        self.assertIn("real, save :: initial_value = 1.0", transformed)

    def test_does_not_add_save_to_parameters_module_data_or_components(self) -> None:
        source = """module example
  real :: module_value = 1.0
  type :: item
    real :: component = 2.0
  end type item
contains
  subroutine work(x)
    real, intent(out) :: x
    real, parameter :: constant = 3.0
    type :: local_item
      real :: component = 4.0
    end type local_item
    x = constant
  end subroutine work
end module example
"""
        transformed, _, _, save_count, _ = transform_text(source)
        self.assertEqual(save_count, 0)
        self.assertEqual(transformed, source)

    def test_preserves_existing_save_attribute(self) -> None:
        source = """subroutine work(x)
  real, intent(out) :: x
  real, save :: initial_value = 1.0
  x = initial_value
end subroutine work
"""
        transformed, _, _, save_count, _ = transform_text(source)
        self.assertEqual(save_count, 0)
        self.assertEqual(transformed, source)

    def test_splits_mixed_declaration_without_saving_automatic_locals(self) -> None:
        source = """subroutine work(x)
  real, intent(out) :: x
  real :: zero = 0.0, temporary, two = 2.0, workspace(4)
  x = zero + two
end subroutine work
"""
        transformed, _, _, save_count, _ = transform_text(source)
        self.assertEqual(save_count, 2)
        self.assertIn("real, save :: zero = 0.0\n", transformed)
        self.assertIn("real :: temporary\n", transformed)
        self.assertIn("real, save :: two = 2.0\n", transformed)
        self.assertIn("real :: workspace(4)\n", transformed)

    def test_handles_simple_continued_initialized_declaration(self) -> None:
        source = """subroutine work(x)
  real, intent(out) :: x
  real :: large = huge(1.0), zero = 0.0, temporary, workspace(4), &
          small = tiny(1.0)
  x = large + zero + small
end subroutine work
"""
        transformed, _, _, save_count, _ = transform_text(source)
        self.assertEqual(save_count, 2)
        self.assertIn("real, save :: large = huge(1.0), zero = 0.0\n", transformed)
        self.assertIn("real :: temporary, workspace(4)\n", transformed)
        self.assertIn("real, save :: small = tiny(1.0)\n", transformed)
        self.assertNotIn("&", transformed)

    def test_preserves_continued_array_constructor(self) -> None:
        source = """module generated_data
  real, parameter :: values(8) = [ &
    1.0, 2.0, 3.0, 4.0, &
    5.0, 6.0, 7.0, 8.0 ]
end module generated_data
"""
        transformed, _, _, _, _ = transform_text(source, max_line_length=40)
        self.assertEqual(transformed, source)
        self.assertEqual(transformed.count("&"), 2)


if __name__ == "__main__":
    unittest.main()
