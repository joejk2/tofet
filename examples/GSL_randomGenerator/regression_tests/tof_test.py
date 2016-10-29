#!/usr/bin/python
"""
Test simple tof simulations
"""
import subprocess as sp
from nose.tools import assert_equal
import os
from test_functions import *


class TestToF(object):


    def setup(self):
        self.run_output = "testing.tmp"
        self.ref_output = "tof/tof_0.out"
        self.description = 'regenerate_occ:  '


    def teardown(self):
        if os.path.exists('testing.tmp'):
            print "get here"
            os.remove('testing.tmp')


    def test_simple_run(self):
        command = 'tft \
                   tof/tof.sim tof/scl.xyz tof/scl.edge'
        assert(run_and_check_sim(command, self.run_output))
 
    def test_tft_run_batch(self):
        command = 'tft_run_batch.py \
                   tof/tof.sim tof/scl.xyz tof/scl.edge'
        run(command)
        outputs = ['0.out', '1.out', '2.out', '3.out'] 
        for output in outputs: 
            assert equal_output(output, 'tof/tof_' + output, '--mob_dis')
            assert equal_output(output, 'tof/tof_' + output, '--mob_vel')
            assert equal_output(output, 'tof/tof_' + output, '--transient')
            assert equal_output(output, 'tof/tof_' + output, '--temp')
            assert equal_output(output, 'tof/tof_' + output, '--fieldZ')
            assert equal_output(output, 'tof/tof_' + output, '--hop_start')
            assert equal_output(output, 'tof/tof_' + output, '--mob_vel')
            assert equal_output(output, 'tof/tof_' + output, '--transient')
            os.remove(output)


    def test_mob_vel(self):
        assert equal_output('0.out', 'tof/tof_0.out', '--hop_start')
        


    def test_mob_vel(self):
        self.test_simple_run()
        assert equal_output(self.run_output, self.ref_output, '--mob_vel')


    def test_mob_dis(self):
        self.test_simple_run()
        assert(equal_output(self.run_output, self.ref_output, '--mob_dis'))

    
    def test_photocurrent(self):
        self.test_simple_run()
        assert equal_output(self.run_output, self.ref_output, '--transient')


    def test_reorg(self):
        self.test_simple_run()
        assert(equal_output(self.run_output, self.ref_output, '--reorg'))


    def test_hop_start(self):
        self.test_simple_run()
        assert(equal_output(self.run_output, self.ref_output, '--hop_start'))


    def test_hop_maxTime(self):
        self.test_simple_run()
        assert(equal_output(self.run_output, self.ref_output, '--maxTime'))


    def test_hop_deltaTime(self):
        self.test_simple_run()
        assert(equal_output(self.run_output, self.ref_output, '--deltaTime'))


    def test_hop_alpha(self):
        self.test_simple_run()
        assert(equal_output(self.run_output, self.ref_output, '--deltaTime'))


    def test_hop_tol(self):
        self.test_simple_run()
        assert(equal_output(self.run_output, self.ref_output, '--deltaTime'))



        

        

        



        

    
    
