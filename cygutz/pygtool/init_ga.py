#!/usr/bin/env python

import sys, glob, warnings
from pygrisb.gutz.init import screen_init
from pygrisb.iface.ifwien import h4set_indmfl


s_init = screen_init()
s_init.run()

# further work if Wien2K+GRISB
if '--nowien' not in sys.argv:
    fstruct = glob.glob('*struct')
    if len(fstruct) != 1:
        warnings.warn(f"none or multiple struct files exist: {fstruct}")
    if '-no-indmfl' not in sys.argv and len(fstruct) == 1:
        h4set_indmfl()
