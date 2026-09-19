import tempfile
import unittest
from pathlib import Path

from lxml import html

from wechat_math import replace_math, render_formula


class MathDeliveryTests(unittest.TestCase):
    def test_inline_and_display_keep_surrounding_text(self):
        doc = html.fromstring(r'<section><p>before <span class="math inline">\(R^2\)</span> after</p><span class="math display">\[x=\frac{a}{b}\]</span></section>')
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(replace_math(doc, Path(tmp)), 2)
            self.assertEqual(''.join(doc.itertext()), 'before  after')
            images = doc.xpath('.//img')
            self.assertIn('inline-block', images[0].get('data-math-style'))
            self.assertIn('display:block', images[1].get('data-math-style'))
            self.assertTrue(all(Path(img.get('src')).is_file() for img in images))

    def test_file_commands_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                render_formula(r'\input{private}', True, Path(tmp))

    def test_malformed_math_stops_delivery(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                replace_math(html.fromstring('<section><span class="math inline">broken</span></section>'), Path(tmp))

    def test_citeproc_native_greek_is_preserved(self):
        doc = html.fromstring('<section>before <span class="math inline"><em>β</em></span>-Catenin</section>')
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(replace_math(doc, Path(tmp)), 0)
            self.assertEqual(''.join(doc.itertext()), 'before β-Catenin')
            self.assertEqual(len(doc.xpath('.//em')), 1)


if __name__ == '__main__':
    unittest.main()
