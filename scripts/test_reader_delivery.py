import unittest
from lxml import html
from build_wechat_review_bundle import remove_wechat_bootstrap, remove_explicit_wechat_omissions
from audit_reader_delivery import FORBIDDEN, PRIVATE, has_decision_recommendation


class ReaderDeliveryTest(unittest.TestCase):
    def test_topic_specific_recommendation_is_allowed(self):
        body = '::: {.callout-important title="先决定要估计哪一种差异"}\n先确定比较对象和归一化分母；只在同一批受试者中比较调整方案。\n:::'
        self.assertTrue(has_decision_recommendation(body))

    def test_empty_title_or_code_example_does_not_pass(self):
        self.assertFalse(has_decision_recommendation('先确定这一点'))
        body = '::: {.callout-tip title="建议"}\n\n:::'
        self.assertFalse(has_decision_recommendation(body))
        self.assertFalse(has_decision_recommendation('```text\n::: {.callout-tip}\n' + '示例' * 20 + '\n:::\n```'))

    def test_setup_id_does_not_delete_data(self):
        root=html.fromstring('<main><section id="sec-setup"><h2>数据</h2><p>90 samples</p><pre>download.file(url, path)</pre></section></main>')
        remove_wechat_bootstrap(root)
        self.assertIn('download.file',root.text_content())
        self.assertIn('90 samples',root.text_content())

    def test_explicit_bootstrap_only_is_removed(self):
        root=html.fromstring('<main><section id="sec-preparation"><h2>数据</h2><p>771 subjects</p><details class="wechat-omit"><summary>安装</summary><pre>install.packages("x")</pre></details></section></main>')
        self.assertEqual(remove_explicit_wechat_omissions(root),1)
        remove_wechat_bootstrap(root)
        self.assertIn('771 subjects',root.text_content())
        self.assertNotIn('install.packages',root.text_content())

    def test_empty_setup_heading_removed(self):
        root=html.fromstring('<main><section id="sec-setup"><h2>环境</h2><details class="wechat-omit"><p>helpers</p></details></section></main>')
        remove_explicit_wechat_omissions(root);remove_wechat_bootstrap(root)
        self.assertFalse(root.xpath('.//h2'))

    def test_legacy_author_narration_is_rejected(self):
        for text in ['旧草稿只取两个队列','真实论文为“PCoA”提供了什么证据起点','附录 A：Methods / Results 模板','这里不复制论文原图']:
            self.assertIsNotNone(FORBIDDEN.search(text))

    def test_public_source_and_scientific_limits_allowed(self):
        self.assertIsNone(FORBIDDEN.search('检出相关序列不等于表型耐药；结果只支持遗传潜力。'))
        self.assertIsNone(PRIVATE.search('https://github.com/petemeng/metagenomics-best-practices'))

if __name__=='__main__':unittest.main()
