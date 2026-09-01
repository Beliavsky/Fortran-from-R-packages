import unittest

from remove_fortran_semicolons import FortranScanError, transform_text


class RemoveFortranSemicolonsTests(unittest.TestCase):
    def transform(self, source: str, keep_simple_case: bool = False) -> str:
        return transform_text(source, keep_simple_case=keep_simple_case).text

    def test_splits_statements_and_preserves_indentation(self) -> None:
        source = "   a = 1; b = 2; call work()\n"
        expected = "   a = 1\n   b = 2\n   call work()\n"
        self.assertEqual(self.transform(source), expected)

    def test_does_not_touch_strings_or_comments(self) -> None:
        source = "s = 'a;b'; t = \"c;d\" ! leave; this; comment\n"
        expected = "s = 'a;b'\nt = \"c;d\" ! leave; this; comment\n"
        self.assertEqual(self.transform(source), expected)

    def test_handles_doubled_quote_escapes(self) -> None:
        source = "s = 'don''t;split'; x = \"a\"\";b\"; y = 3\n"
        expected = "s = 'don''t;split'\nx = \"a\"\";b\"\ny = 3\n"
        self.assertEqual(self.transform(source), expected)

    def test_preserves_semicolon_in_continued_character_literal(self) -> None:
        source = "s = 'first&\n&;second'; x = 1\n"
        expected = "s = 'first&\n&;second'\nx = 1\n"
        self.assertEqual(self.transform(source), expected)

    def test_moves_trailing_comment_to_last_statement(self) -> None:
        source = "  x = 1; y = 2 ! result; retained\n"
        expected = "  x = 1\n  y = 2 ! result; retained\n"
        self.assertEqual(self.transform(source), expected)

    def test_removes_empty_and_trailing_statements(self) -> None:
        source = "  ; x = 1;; y = 2;\n"
        expected = "  x = 1\n  y = 2\n"
        self.assertEqual(self.transform(source), expected)

    def test_can_keep_simple_case_branch(self) -> None:
        source = 'case ("sin"); y = sin(x)\n'
        self.assertEqual(self.transform(source, keep_simple_case=True), source)
        self.assertEqual(self.transform(source), 'case ("sin")\ny = sin(x)\n')

    def test_leaves_preprocessor_line_unchanged(self) -> None:
        source = "#define RUN(x) x; x\na = 1; b = 2\n"
        expected = "#define RUN(x) x; x\na = 1\nb = 2\n"
        self.assertEqual(self.transform(source), expected)

    def test_preserves_crlf_and_missing_final_newline(self) -> None:
        self.assertEqual(self.transform("a=1; b=2\r\n"), "a=1\r\nb=2\r\n")
        self.assertEqual(self.transform("a=1; b=2"), "a=1\nb=2")

    def test_is_idempotent(self) -> None:
        once = self.transform("a=1; b='x;y'; c=3 ! d;e\n")
        self.assertEqual(self.transform(once), once)

    def test_rejects_unterminated_character_literal(self) -> None:
        with self.assertRaises(FortranScanError):
            self.transform("s = 'unterminated\n")


if __name__ == "__main__":
    unittest.main()
