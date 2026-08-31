import unittest

from generate_task_views import clean_heading, group_packages, metadata


class GenerateTaskViewsTests(unittest.TestCase):
    def test_cleans_linked_heading_and_anchor(self):
        self.assertEqual(clean_heading("[Base functionality]{#base}"), "Base functionality")

    def test_reads_source_metadata(self):
        source = "name: Demo\ntopic: Demo Methods\nmaintainer: A. Author\nversion: 2026-01-02\n"
        self.assertEqual(metadata(source)["topic"], "Demo Methods")

    def test_groups_first_non_skipped_package_occurrence(self):
        source = """### Table of contents
`r pkg("alpha")`
# Discrete distributions
`r pkg("alpha")` and `r pkg("beta")`
# Continuous distributions
`r pkg("alpha")`
"""
        groups, missing = group_packages(
            source,
            {"alpha": "Alpha", "beta": "beta"},
            frozenset({"Introduction", "Table of contents"}),
        )
        self.assertEqual(groups["Discrete distributions"], ["Alpha", "beta"])
        self.assertEqual(missing, set())

    def test_groups_heading_free_view_under_custom_initial_section(self):
        groups, missing = group_packages(
            '`r pkg("jomo")` and `r pkg("pan")`',
            {"jomo": "jomo", "pan": "pan"},
            frozenset({"Links"}),
            "Methods and packages",
        )
        self.assertEqual(groups["Methods and packages"], ["jomo", "pan"])
        self.assertEqual(missing, set())


if __name__ == "__main__":
    unittest.main()
