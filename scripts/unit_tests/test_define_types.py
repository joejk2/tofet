from io import print_lists_to_columns
import os
from nose.tools import assert_equals


class Test_tofet_define_types(object):
    def setUp(self):
        x = [0] * 5
        y = [0] * 5
        z = [0, 1, 2, 3, 4]
        E = [3] * 5
        self.file_with_energies = "test_with_energies.tmp"
        self.file_without_energies = "test_without_energies.tmp"
        self.output = "test_out.tmp"
        print_lists_to_columns(x, y, z,
                               filename=self.file_without_energies)
        print_lists_to_columns(x, y, z, E,
                               filename=self.file_with_energies)
    def tearDown(self):
        os.remove(self.file_without_energies)
        os.remove(self.file_with_energies)
    def test_without_energies(self):
        """Test get correct number of collectors and generators from file with
        energies"""
        command = "../tft_define_types.py " + self.file_without_energies +\
                  " 0.11" 
        output_lines = os.popen(command).readlines()
        c, g = 0, 0 
        for line in output_lines:
            if "c" in line:
                c += 1
            if "g" in line:
                g += 1
        assert_equals(c, 1)
        assert_equals(g, 1)
    def test_with_energies(self):
        """Test get correct number of collectors and generators from file with
        energies"""
        command = "../tft_define_types.py -e " + self.file_with_energies +\
                  " 0.11" 
        print command
        output_lines = os.popen(command).readlines()
        c, g = 0, 0 
        for line in output_lines:
            if "c" in line:
                c += 1
            if "g" in line:
                g += 1
        assert_equals(c, 1)
        assert_equals(g, 1)

        
