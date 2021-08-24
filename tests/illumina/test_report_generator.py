import unittest
from illumina.report_generator import read_bcl2fastq_stats_data_from_pandas

class Report_generator1(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read_bcl(self):
        a = 1
        self.assertEqual(a, 1)

    def test_read_bcl2(self):
        self.assertTrue(False)

if __name__ == '__main__':
  unittest.main()