import unittest
from audit_chapter_prose import prose_errors

class ProseRegressionTests(unittest.TestCase):
    def test_valid_math_and_code_regex(self):
        self.assertEqual(prose_errors('Use $R^2$.\n```{r}\npattern <- "\\("\n```\n`\\[code\\]`'),[])
    def test_reject_broken_math(self):
        self.assertIn('legacy_math_delimiters',prose_errors(r'Use \(R^2\).'))
    def test_group_and_maintainer_regressions(self):
        self.assertIn('wrong_wetland_groups',prose_errors('conventional/intensive/temporal wetland'))
        self.assertIn('maintainer_prose',prose_errors('这是教程合同'))

if __name__=='__main__':unittest.main()
