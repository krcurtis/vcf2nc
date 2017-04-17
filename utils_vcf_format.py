# Copyright 2017 Fred Hutchinson Cancer Research Center
################################################################################
### Functional tests for vcf2nc

from __future__ import print_function

import argparse
import gzip


from math import *


import h5py
from netCDF4 import Dataset


import numpy as np


import gecco


def convert_matrix_to_strings(M):
    """convert numpy character matrix to list of strings"""
    strings = [];
    (nrow, ncol) = M.shape;
    for i in range(nrow):
        v = M[i];
        s = v.tostring()
        s = s.strip("\x00")
        strings.append(s);
    return strings;


class UtilsVCFFormat:
    def __init__(self, n_samples, n_snps):
        self.n_samples = n_samples
        self.n_snps = n_snps
        
        self.snp_name = [ 'Mito' + str(i+1) for i in range(n_snps)]
        self.sample_id = ['HSample' + str(i+1) for i in range(n_samples)]
        self.chromosome = np.random.randint(1,26, n_snps)
        self.position = np.random.randint(1, 1000, n_snps)
        self.ref = np.random.choice(['A', 'C', 'G', 'T'], size=n_snps)
        self.alt1 = np.random.choice(['A', 'C', 'G', 'T'], size=n_snps)
        self.alt2 = np.random.choice(['A', 'C', 'G', 'T'], size=n_snps)
        self.alt3 = np.random.choice(['A', 'C', 'G', 'T'], size=n_snps)

        self.Qual = np.random.uniform(size=(n_snps,))

        self.Filter = [ 'A;PASS' + str(i+1) for i in range(n_snps)]
        self.infoSB = np.random.randint(1, 1000, n_snps)
        self.infoRD = np.random.randint(10, 50, n_snps)
        self.infoBQ = np.random.randint(100, 200, n_snps)


        self.allele1_category = np.random.randint(1,3, size=(n_samples, n_snps))
        self.allele2_category = np.random.randint(1,3, size=(n_samples, n_snps))
        self.genotype_phase = np.random.randint(1,4, size=(n_samples, n_snps))
        self.read_depth = np.random.randint(10,50, size=(n_samples, n_snps))
        
        self.likelihoodAA = np.random.uniform(size=(n_samples, n_snps))
        self.likelihoodAB = np.random.uniform(size=(n_samples, n_snps))
        self.likelihoodBB = np.random.uniform(size=(n_samples, n_snps))

    def write_header(self, f_out):
        f_out.write( "##fileformat=VCFv4.1\n");
        f_out.write( "##INFO=<ID=%s,Number=%i,Type=%s,Description=\"Some value\">\n" % ("SB", 1, "Integer"))
        f_out.write( "##INFO=<ID=%s,Number=%i,Type=%s,Description=\"Some value\">\n" % ("RD", 1, "Integer"))
        f_out.write( "##INFO=<ID=%s,Number=%i,Type=%s,Description=\"Some value\">\n"%  ("BQ", 1, "Float"))
        f_out.write( "##FILTER=<ID=q10,Description=\"Quality below some level\">\n");
        f_out.write( "##FORMAT=<ID=%s,Number=%i,Type=%s,Description=\"Some value\">\n" % ("GT", 1, "String"))
        f_out.write( "##FORMAT=<ID=%s,Number=%i,Type=%s,Description=\"Some value\">\n" % ("RD", 1, "Integer"))
        f_out.write( "##FORMAT=<ID=%s,Number=%i,Type=%s,Description=\"Some value\">\n" % ("PL", 3, "Float"))

        f_out.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");


        for i in range(self.n_samples):
            f_out.write("\t" + self.sample_id[i]);
        f_out.write( "\n");

    def write_vcf(self, output_file):
        f_out = open(output_file, "w")
        self.write_header(f_out)

        lookup_phase = "-|\\/";  # the '\' is defined in VCF v3.3 but dropped in v4.0 and v4.1

        chromosome = [ "-", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                               "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", 
                               "21", "22", "X", "Y", "XY", "MT", "U"]


        for k in range(self.n_snps):
            f_out.write("{chrom}\t{pos}\t{snpid}\t{ref}\t{alt1},{alt2},{alt3}\t{qual}\t{Filter}".format(
                        chrom=chromosome[self.chromosome[k]],
                        pos=self.position[k],
                        snpid=self.snp_name[k],
                        ref=self.ref[k],
                        alt1=self.alt1[k],
                        alt2=self.alt2[k],
                        alt3=self.alt3[k],
                        qual=self.Qual[k],
                        Filter=self.Filter[k]))
            # info collection
            f_out.write("\tSB=%i;RD=%i;BQ=%lf" % (self.infoSB[k], self.infoRD[k], self.infoBQ[k]))

            #format description
            f_out.write("\tGT:RD:PL");

            # data
            for i in range(self.n_samples):
                f_out.write("\t%i%c%i:%i:%lf,%lf,%lf" % 
                            (self.allele1_category[i,k]-1,             #apply mapping for file generation
                             lookup_phase[ self.genotype_phase[i,k] ],  #apply mapping for file generation
                             self.allele2_category[i,k]-1,             #apply mapping for file generation
                             self.read_depth[i,k],
                             self.likelihoodAA[i,k],
                             self.likelihoodAB[i,k],
                             self.likelihoodBB[i,k]))
            f_out.write("\n")
        f_out.close()


    def compare_vector(self, input_netcdf, varname, b, epsilon=None):
        input_ncvars = Dataset(input_netcdf, 'r', format='NETCDF4')
        a = input_ncvars.variables[varname][:]

        # auto convert character arrays to strings
        if '|S1' == a.dtype and 2 == len(a.shape):
            a = convert_matrix_to_strings(a)

        success = True
        try:
            if len(a) != len(b):
                success = False
                print("ERROR different sizes {a} {b}".format(a=len(a), b=len(b)))
                raise Exception("ERROR different sizes {a} {b}".format(a=len(a), b=len(b)))

            if None == epsilon:
                for i in range(len(a)):
                    if a[i] != b[i]:
                        success = False
                        print("ERROR compare failed at index {i}: {a} != {b}".format(i=i,a=a[i], b=b[i]))
                        raise Exception("ERROR compare failed at index {i}: {a} != {b}".format(i=i,a=a[i], b=b[i]))
            else:
                for i in range(len(a)):
                    if np.abs(a[i] - b[i]) > epsilon:
                        success = False
                        print("ERROR compare failed at index {i}: {a} != {b}".format(i=i,a=a[i], b=b[i]))
                        raise Exception("ERROR compare failed at index {i}: {a} != {b}".format(i=i,a=a[i], b=b[i]))

        finally:
            input_ncvars.close()

        return success




    def compare_matrix_values(self, a, b, epsilon=None, msg=None):
        msg="" if None == msg else " " + msg
        if None == epsilon:
            for i in range(a.shape[0]):
                for j in range(a.shape[1]):
                    if a[i,j] != b[i,j]:
                        print("ERROR compare failed at index {i},{j}: {a} vs {b}{msg}".format(i=i,j=j,a=a[i,j], b=b[i,j], msg=msg))
                        return False

        else:
            for i in range(a.shape[0]):
                for j in range(a.shape[1]):
                    if np.abs(a[i,j] -  b[i,j]) > epsilon:
                        print("ERROR compare failed at index {i},{j}: {a} vs {b}{msg}".format(i=i,j=j,a=a[i,j],b=b[i,j], msg=msg))
                        return False
        return True

    def compare_matrix(self, input_netcdf, varname, b, epsilon=None):
        input_ncvars = Dataset(input_netcdf, 'r', format='NETCDF4')
        a = input_ncvars.variables[varname][:]

        if a.shape != b.shape:
            print("ERROR different sizes {a} {b}".format(a=a.shape, b=b.shape))
            return False

        result = self.compare_matrix_values(a, b, epsilon)
        input_ncvars.close()
        return result



    def compare_three_matrix(self, input_netcdf, varname, A, B, C, epsilon=None):
        input_ncvars = Dataset(input_netcdf, 'r', format='NETCDF4')
        a = input_ncvars.variables[varname][:]

        shape = (a.shape[0], a.shape[1])
        if shape != A.shape or shape != B.shape or shape != C.shape or 3 != a.shape[2]:
            print("ERROR different sizes {a} {A} {B} {C}".format(a=a.shape, A=A.shape, B=B.shape, C=C.shape))
            return False

        result1 = self.compare_matrix_values(a[:,:,0], A, epsilon, msg="Slice1")
        result2 = self.compare_matrix_values(a[:,:,1], B, epsilon, msg="Slice2")
        result3 = self.compare_matrix_values(a[:,:,2], C, epsilon, msg="Slice3")
        input_ncvars.close()
        return result1 and result2 and result3






    def compare_variables(self, input_file):
        epsilon = 1e-6;

        # load and compare

        if not self.compare_vector(input_file, "Sample_ID", self.sample_id):
            return False

        if not self.compare_vector(input_file, "Chromosome", self.chromosome):
            return False

        if not self.compare_vector(input_file, "Position", self.position):
            return False

        if not self.compare_vector(input_file, "ID", self.snp_name):
            return False


        if not self.compare_vector(input_file, "Reference_Allele", self.ref):
            return False
        if not self.compare_vector(input_file, "Alternate1_Allele", self.alt1):
            return False
        if not self.compare_vector(input_file, "Alternate2_Allele", self.alt2):
            return False
        if not self.compare_vector(input_file, "Alternate3_Allele", self.alt3):
            return False


        if not self.compare_vector(input_file, "info_SB", self.infoSB):
            return False

        if not self.compare_vector(input_file, "info_RD", self.infoRD):
            return False
        if not self.compare_vector(input_file, "info_BQ", self.infoBQ, 1e-6):
            return False

        if not self.compare_three_matrix(input_file, "array_GT", self.allele1_category, self.genotype_phase, self.allele2_category):
            return False

        if not self.compare_matrix(input_file, "array_RD", self.read_depth):
            return False

        if not self.compare_three_matrix(input_file, "array_PL", self.likelihoodAA, self.likelihoodAB, self.likelihoodBB, epsilon):
            return False


        if not self.compare_vector(input_file, "FILTER", self.Filter):
            return False

        return True

