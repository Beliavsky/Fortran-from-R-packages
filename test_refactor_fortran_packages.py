#!/usr/bin/env python3
"""Unit tests for the package refactoring driver."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from refactor_fortran_packages import (
    dependency_order,
    discover_packages,
    has_numeric_statement_labels,
    maintained_sources,
    xparam_command,
)


class RefactorDriverTests(unittest.TestCase):
    def test_xparam_pipeline_uses_initialized_only_without_backups(self) -> None:
        command = xparam_command(Path("xparam.py"), [Path("source.f90")])
        self.assertIn("--initialized-only", command)
        self.assertIn("--no-backup", command)
        self.assertEqual(command[-1], "source.f90")

    def test_dependency_order_places_dependencies_first(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "base").mkdir()
            (root / "base" / "fpm.toml").write_text('name = "base"\nversion = "0.1.0"\n', encoding="utf-8")
            (root / "user").mkdir()
            (root / "user" / "fpm.toml").write_text(
                'name = "user"\nversion = "0.1.0"\n[dependencies]\nbase = { path = "../base" }\n',
                encoding="utf-8",
            )
            packages = discover_packages(root)
            names = [package.name for package in dependency_order(packages)]
            self.assertLess(names.index("base"), names.index("user"))

    def test_dependency_order_prioritizes_widely_used_foundations(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in ("foundation", "unused", "first_user", "second_user"):
                (root / name).mkdir()
            (root / "foundation" / "fpm.toml").write_text(
                'name = "foundation"\nversion = "0.1.0"\n', encoding="utf-8"
            )
            (root / "unused" / "fpm.toml").write_text(
                'name = "unused"\nversion = "0.1.0"\n', encoding="utf-8"
            )
            for name in ("first_user", "second_user"):
                (root / name / "fpm.toml").write_text(
                    f'name = "{name}"\nversion = "0.1.0"\n'
                    '[dependencies]\nfoundation = { path = "../foundation" }\n',
                    encoding="utf-8",
                )
            names = [package.name for package in dependency_order(discover_packages(root))]
            self.assertEqual(names[0], "foundation")
            self.assertLess(names.index("unused"), names.index("first_user"))
            self.assertLess(names.index("unused"), names.index("second_user"))

    def test_sources_exclude_vendor_upstream_build_and_fixed_form(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            package = Path(temporary)
            expected = package / "src" / "kept.f90"
            expected.parent.mkdir()
            expected.write_text("module kept\nend module kept\n", encoding="utf-8")
            for directory in ("vendor", "upstream", "build"):
                path = package / directory / "ignored.f90"
                path.parent.mkdir()
                path.write_text("ignored\n", encoding="utf-8")
            (package / "src" / "fixed.f").write_text("      end\n", encoding="utf-8")
            self.assertEqual(maintained_sources(package), [expected])

    def test_packages_without_any_dependencies_precede_external_dependencies(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "independent").mkdir()
            (root / "independent" / "fpm.toml").write_text(
                'name = "independent"\nversion = "0.1.0"\n', encoding="utf-8"
            )
            (root / "external_user").mkdir()
            (root / "external_user" / "fpm.toml").write_text(
                'name = "external_user"\nversion = "0.1.0"\n'
                '[dependencies]\nexternal = { git = "https://example.com/external.git" }\n',
                encoding="utf-8",
            )
            packages = discover_packages(root)
            ordered = dependency_order(packages)
            self.assertEqual([package.name for package in ordered], ["independent", "external_user"])
            self.assertEqual(ordered[1].dependencies, ("external",))

    def test_detects_numeric_statement_labels_that_fprettify_can_remove(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "legacy_control_flow.f90"
            source.write_text(
                "subroutine step()\n"
                "  integer :: iteration = 60\n"
                "  goto 60\n"
                "60 call advance()\n"
                "end subroutine step\n",
                encoding="utf-8",
            )
            self.assertTrue(has_numeric_statement_labels(source))

    def test_does_not_treat_numbers_in_code_or_comments_as_labels(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "modern.f90"
            source.write_text(
                "subroutine step()\n"
                "  integer :: iteration = 60\n"
                "  ! 60 is only mentioned in a comment\n"
                "end subroutine step\n",
                encoding="utf-8",
            )
            self.assertFalse(has_numeric_statement_labels(source))


if __name__ == "__main__":
    unittest.main()
