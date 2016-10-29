#!/usr/bin/python
"""
Test track hoppers
"""
import subprocess as sp
from nose.tools import assert_equal
import os
from test_functions import *


class TestTrackHoppers(object):


    def setup(self):
        self.run_output = "testing.tmp"
        self.ref_output = "track_hoppers/track.out"


    def teardown(self):
        os.remove('testing.tmp')


    def test_run(self):
        command = 'tftOccupation \
                  track_hoppers/track.sim track track_hoppers/scl.xyz track_hoppers/scl.edge'
        assert(run_and_check_sim(command, self.run_output))

    def test_exact_equality_of_file(self):
        self.test_run()
        ref = open(self.ref_output).readlines()
        run = open(self.run_output).readlines()
        assert_equal(ref[30:], run[30:])

