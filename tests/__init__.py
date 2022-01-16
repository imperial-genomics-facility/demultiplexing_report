import unittest

def get_tests():
    return full_suite()

def full_suite():
    from .illumina.test_report_generator import Report_generator1
    from .illumina.test_report_generator_for_db import Test_report_generator_for_db
    return unittest.TestSuite([
        unittest.TestLoader().loadTestsFromTestCase(t)
        for t in [
            Report_generator1,
            Test_report_generator_for_db
            ]
    ])