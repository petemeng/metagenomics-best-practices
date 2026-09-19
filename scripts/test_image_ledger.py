import unittest

from lxml import html

from build_wechat_review_bundle import image_ledger_matches


class ImageLedgerTests(unittest.TestCase):
    def test_repeated_formula_keeps_both_occurrences_with_one_asset(self):
        doc = html.fromstring('<section><img src="images/formula.jpg"><p>Again:</p><img src="images/formula.jpg"></section>')
        records = [{"relative_src": "images/formula.jpg"}]
        self.assertTrue(image_ledger_matches(doc, records, 1))
        self.assertEqual(len(doc.xpath('.//img')), 2)

    def test_missing_unused_or_duplicate_asset_fails(self):
        doc = html.fromstring('<section><img src="images/a.jpg"></section>')
        records = [{"relative_src": "images/a.jpg"}]
        self.assertFalse(image_ledger_matches(doc, [], 0))
        self.assertFalse(image_ledger_matches(doc, [{"relative_src": "images/b.jpg"}], 1))
        self.assertFalse(image_ledger_matches(doc, records * 2, 2))
        self.assertFalse(image_ledger_matches(doc, records, 2))

    def test_remote_image_does_not_need_a_local_upload_record(self):
        doc = html.fromstring('<section><img src="https://example.org/figure.png"><img src="images/a.jpg"></section>')
        self.assertTrue(image_ledger_matches(doc, [{"relative_src": "images/a.jpg"}], 1))


if __name__ == '__main__':
    unittest.main()
