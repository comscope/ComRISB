#!/usr/bin/env python
import argparse, h5py, sys
from pygrisb.gsolver.gs_driver import driver as gs_driver
from pygrisb.run.cygutz import get_cygtail, get_groot

def solve_hembed(imp=0,
        edinit=0,
        mpi="mpirun -np 1",
        nval=5,
        analysis=False,
        cmdlargs=False,
        iembeddiag=None,
        ):
    # imp is needed first
    if "-i" in sys.argv:
        imp_pos = sys.argv.index("-i") + 1
    elif "--imp" in sys.argv:
        imp_pos = sys.argv.index("--imp") + 1
    else:
        imp_pos = None
    if imp_pos is not None:
        imp = int(sys.argv[imp_pos])

    with h5py.File("GParam.h5", "r") as f:
        if iembeddiag is None:
            iembeddiag = f["./"].attrs["giembeddiag"]
        nval_bot = f["./"].attrs["nval_bot_list"][imp]
        nval_top = f["./"].attrs["nval_top_list"][imp]
    g_root = get_groot()

    # taking command line argument
    if cmdlargs:
        parser = argparse.ArgumentParser()
        parser.add_argument("-i", "--imp", type=int, default=imp,
                help="impurity index zero-based (int)")
        parser.add_argument("-d", "--iembeddiag",
                type=int, default=iembeddiag,
                help="choice embedding h solver (int)")

        parser.add_argument("-n", "--nval", type=int, default=nval,
                help="valence (int)")
        parser.add_argument("-m", "--mpi", type=str, default=mpi,
                help="mpirun prefix (str)")

        parser.add_argument("-p", "--path", type=str, default=g_root,
                help="executible path (str)")

        parser.add_argument("-vb", "--valbot", type=int, default=nval_bot,
                help="valence bottom (int)")
        parser.add_argument("-vt", "--valtop", type=int, default=nval_top,
                help="valence top (int)")
        parser.add_argument("-e", "--edinit", type=int, default=0,
                help="start with previous solution in ED (int):"+ \
                " 0-> yes, 1->no.")
        args = parser.parse_args()
        args.mpi = args.mpi.split()
        iembeddiag = args.iembeddiag
        g_root = args.path
        mpi = args.mpi
        nval = args.nval
        nval_bot = args.valbot
        nval_top = args.valtop
        edinit = args.edinit

    gs_driver(imp=imp,
            iembeddiag=iembeddiag,
            path=g_root,
            mpiexec=mpi,
            nval=nval,
            nval_bot=nval_bot,
            nval_top=nval_top,
            edinit=edinit,
            analysis=analysis,
            cygtail=get_cygtail(),
            )


if __name__ == "__main__":
    solve_hembed(cmdlargs=True)
