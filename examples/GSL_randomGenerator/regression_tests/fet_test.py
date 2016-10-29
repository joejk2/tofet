#!/usr/bin/python
"""
Test simple tof simulations
"""
from nose.tools import assert_equal, assert_true, assert_false
import os
from test_functions import *


class TestFeT(object):
    
    
    def setup(self):
        self.run_output = "testing.tmp"
        self.ref_output = "fet/regression.out"


    def teardown(self):
        if os.path.exists('testing.tmp'):
            os.remove('testing.tmp')
    

    def test_run(self):
        command = 'tftOccupation fet/fet.sim \
                 fet/scl_fet.edge\
        fet/scl_fet.xyz'
        assert(run_and_check_sim(command, self.run_output))
        

    def test_occupation(self):
        self.test_run()
        assert equal_output(self.ref_output, self.run_output, '--occ_prob')


    def test_energies(self):
        self.test_run()
        assert equal_output(self.ref_output, self.run_output,
                                  '--energies')


    def test_fail_occ_mol(self):
        self.test_run()
        assert equal_output(self.ref_output, self.run_output, '--occ_mol')


    def test_sat_mu(self):
        output = run_return_stdout('tft_calc_sat_mu.py -l 99 -w 1\
                               ../fet/fet*.out').split('\n')
        assert_equal(output[0], "0.449783903961  cm^2 / V.s")


    def test_ecp(self):
        run = run_return_stdout('tft_calc_ecp.py \
                        ../fet_electrochemicalpotential/output.out')
        ref = open('../fet_electrochemicalpotential/output.ecPot').read()
        assert_equal(len(run), len(ref))
        assert_equal(run, ref)

    def test_average_xy(self):
        run = run_return_stdout('tft_average_xy.py\
                        ../fet_electrochemicalpotential/output.ecPot 100')
        ref = open('../fet_electrochemicalpotential/output.ecPot_zAvg').read()
        assert_equal(run, ref)


    
