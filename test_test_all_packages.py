import tempfile
import unittest
from pathlib import Path

from test_all_packages import (
    discover_packages,
    package_has_tests,
    package_test_prerequisites,
    parse_args,
    safe_log_name,
    select_packages,
)


class TestAllPackagesTests(unittest.TestCase):
    def test_discovers_only_top_level_fpm_packages_case_insensitively_sorted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ("zeta", "Alpha"):
                package = root / name
                package.mkdir()
                (package / "fpm.toml").write_text(f'name = "{name}"\n', encoding="utf-8")
            nested = root / "container" / "nested"
            nested.mkdir(parents=True)
            (nested / "fpm.toml").write_text('name = "nested"\n', encoding="utf-8")
            self.assertEqual([path.name for path in discover_packages(root)], ["Alpha", "zeta"])

    def test_detects_conventional_test_source(self):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory)
            (package / "fpm.toml").write_text('name = "demo"\n', encoding="utf-8")
            (package / "test").mkdir()
            (package / "test" / "main.f90").write_text("program main\nend program\n", encoding="utf-8")
            self.assertTrue(package_has_tests(package))

    def test_respects_disabled_auto_tests(self):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory)
            (package / "fpm.toml").write_text(
                'name = "demo"\n[build]\nauto-tests = false\n', encoding="utf-8"
            )
            (package / "test").mkdir()
            (package / "test" / "main.f90").write_text("program main\nend program\n", encoding="utf-8")
            self.assertFalse(package_has_tests(package))

    def test_reports_missing_declared_test_prerequisite(self):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory)
            (package / "fpm.toml").write_text(
                'name = "demo"\n[extra.r-fpm-test]\n'
                'required-files = ["vendor/lib/backend.a"]\n'
                'preparation = "run the backend builder"\n',
                encoding="utf-8",
            )
            missing, preparation = package_test_prerequisites(package)
            self.assertEqual(missing, [package / "vendor/lib/backend.a"])
            self.assertEqual(preparation, "run the backend builder")

    def test_selects_requested_packages_without_duplicates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            package = root / "Mixtools"
            package.mkdir()
            (package / "fpm.toml").write_text('name = "mixtools"\n', encoding="utf-8")
            selected, missing = select_packages(root, ["mixtools", "MIXTOOLS", "absent"])
            self.assertEqual([path.name for path in selected], ["Mixtools"])
            self.assertEqual(missing, ["absent"])

    def test_makes_safe_log_names(self):
        self.assertEqual(safe_log_name("package/name:one"), "package_name_one")

    def test_cleans_all_build_products_by_default(self):
        self.assertEqual(parse_args([]).clean, "all")

    def test_cleanup_mode_can_be_changed(self):
        self.assertEqual(parse_args(["--clean", "skip"]).clean, "skip")


if __name__ == "__main__":
    unittest.main()
