################################################################################
### Unit tests for vcf2nc, implemented in Python


from __future__ import print_function

import os
import unittest


import utils_vcf_format

class ConvertTest(unittest.TestCase):
    def test_1(self):
        uf = utils_vcf_format.UtilsVCFFormat(1,1)

        test_vcf = os.path.join(os.environ['HOME'], "tmp/test.vcf")
        test_netcdf = os.path.join(os.environ['HOME'], "tmp/test.nc")


        uf.write_vcf(test_vcf)
        os.system('rm ' + test_netcdf)
        cmd = ' '.join([ "./vcf2nc",
                             "-o", test_netcdf,
                             "-i", test_vcf])
        #print(cmd)
        os.system(cmd)
        self.assertTrue(uf.compare_variables(test_netcdf))

    def test_2(self):
        uf = utils_vcf_format.UtilsVCFFormat(10,20)

        test_vcf = os.path.join(os.environ['HOME'], "tmp/test.vcf")
        test_netcdf = os.path.join(os.environ['HOME'], "tmp/test2.nc")

        uf.write_vcf(test_vcf)
        os.system('rm ' + test_netcdf)
        cmd = ' '.join([ "./vcf2nc",
                             "-o", test_netcdf,
                             "-i", test_vcf])
        #print(cmd)
        os.system(cmd)
        self.assertTrue(uf.compare_variables(test_netcdf))




