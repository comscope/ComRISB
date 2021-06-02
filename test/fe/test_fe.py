#!/usr/bin/env python

import unittest, tempfile, os, shutil, subprocess
from pygrisb.basic.data import compare_array as compare_data


class KnowValues(unittest.TestCase):
    def test_data_prepare(self):
        cmd = ['./run.sh']
        subprocess.run(cmd)

    def test_nphy_matrix(self):
        self.assertTrue(compare_data('/impurity_0/NC_PHY',
            fname1='GLog.h5',
            fname2='dftg/u0j0/lowh/GLog.h5',
            chk_sum=True,   # avoid random rotation effect.
            )
            )


if __name__ == "__main__":
    print("Test ComRISB calculation for Fe with U=J=0")
    cwd = os.getcwd()
    tempd = tempfile.mkdtemp(prefix='gtest_tmp_')
    for f in ["bcc.cif", "comrisb.ini", "ginit.json", "run.sh", "GLog.h5"]:
        shutil.copy(f, tempd)
    os.chdir(tempd)
    unittest.main()
