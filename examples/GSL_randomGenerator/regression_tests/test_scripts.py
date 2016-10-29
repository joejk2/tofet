#!/usr/bin/python

import subprocess as sp

def is_script_on_path(script):
    try:
        # Catch stdout and stderr.
        proc = sp.Popen(script, stdout=sp.PIPE,
                       stderr=sp.PIPE)
    except OSError: 
        print "OSError raised by", script
        return False
    return True


def test_scripts_on_path():
    scripts = ['ls', 'tft', 'tftOccupation', 'tft_average_xy.py',
               'tft_calc_sat_mu.py', 'tft_get_series.py', 
               'tft_make_cubic_lattice.py', 'tft_rotate_coords.py',
               'tft_calc_ecp.py', 'tft_define_types.py',
                'tft_get_values.py', 'tft_plot_pf.py', 'tft_run_batch.py']
    for script in scripts:
        print "Testing", script, "is on your path..."
        assert(is_script_on_path(script))








        
