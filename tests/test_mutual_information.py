import os
import unittest
from ..bitools import mutual_information as mi


class TestMutualInformation(unittest.TestCase):
    def test_compute_mutual_information(self):
        mutual_information_file = os.path.join('tests', 'test_mutual_information.fa')
        base_pair_count, mean_mutual_information = mi.compute_mutual_information(mutual_information_file)
        self.assertEqual(3, base_pair_count)
        self.assertEqual(1, mean_mutual_information)


if __name__ == '__main__':
    unittest.main()
