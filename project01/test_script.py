from regex_script import read_file as rf_regex
from regex_script import parse_function
from project01 import read_file as rf_group
from project01 import parse_line as parse_line
import unittest


class test_io(unittest.TestCase):

    def test_equality(self):
        '''
        Test whether both functions find the same entries and values
        '''
        regex_result = rf_regex()
        group_result = rf_group("clinvar_20190923_short.vcf")
        self.assertDictEqual(regex_result, group_result)
    
    def test_threshold(self):
        '''
        Test whether functions find same entries/values at different threshold 
        '''
        t = 0.001
        regex_result = rf_regex(threshold=t)
        group_result = rf_group("clinvar_20190923_short.vcf", threshold=t)
        self.assertDictEqual(regex_result, group_result)
    
    def test_parse_line(self):
        '''
        Test that same correct disease found
        '''
        diag = "CRS|COOF"
        correct = diag.split("|")
        test_line = f"AF_EXAC=0;Dead_River_Company=BangorME;CLNDN={diag}"
        regex_result = parse_function(test_line, "AF_EXAC", "CLNDN")
        group_result = parse_line(test_line, 0.0001)
        self.assertEqual(regex_result, correct)
        self.assertEqual(group_result, correct)
    
    def test_parse_line2(self):
        '''
        Test that disease not found
        '''
        diag = "CRS|COOF"
        correct = diag.split("|")
        test_line = f"Dead_River_Company=BangorME;CLNDN={diag}"
        regex_result = parse_function(test_line, "AF_EXAC", "CLNDN")
        group_result = parse_line(test_line, 0.0001)
        self.assertNotEqual(regex_result, correct)
        self.assertNotEqual(group_result, correct)


if __name__ == "__main__":
    unittest.main()