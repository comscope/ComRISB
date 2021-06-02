import h5py, json, sys, os, numpy, argparse
from pygrisb.gutz.init import screen_init


def modify_gparam():
    '''convenient scipt to modify the settings in GPARAM.h5 file.
    '''
    parser = argparse.ArgumentParser(
            description="script to modify cygutz settings.")
    parser.add_argument("--iembeddiag", type=int, default=None,
            help="hembed diag method (int)"+ \
            "    -1: Valence truncation ED (VTED); "+ \
            "    -2: VTED with Sz symmetry; "+ \
            "    -3: VTED with S=0; "+ \
            "    -4: VTED with Jz symmetry; "+ \
            "   -10: machine mearning (ml), kernel-ridge regression; "+ \
            "   -11: ml, normal mode expansion (nm) 1st order; "+ \
            "   -13: ml, nm 2nd order" + \
            "   -14: ml, nm 3rd order")
    parser.add_argument("--dc_mode", type=int, default=None,
            help="double couinting mode (int) "+ \
            "   12: FLL dc (updated at charge iteration); "+ \
            "    2: Fix dc potential; "+ \
            "    1: FLL dc self-consistent; "+ \
            "    0: No dc")
    parser.add_argument("--u_type", type=int, default=None,
            help="u-matrix type (int) "+ \
            "    1: Slater-Condo with [U,J]; "+ \
            "    3: Slater-Condo with [F0,F2,...]; "+ \
            "    2: Kanamori with [U,J]")
    parser.add_argument("--dc_nelf_list", type=str, default=None,
            help="nf list for dc of all impurities 'n1_up n1_dn n2_up...'" + \
            " (str), empty for del")
    parser.add_argument("--unique_j_ev", type=str, default=None,
            help="list of hund j for unique coor imp 'j1 j2 ...' (str)")
    parser.add_argument("--unique_u_ev", type=str, default=None,
            help="list of hubbard u for unique coor imp 'u1 u2 ...' (str)")
    parser.add_argument("--nval_bot_list", type=str, default=None,
            help="list of minimal valence 'n1 n2 ...' (str)")
    parser.add_argument("--nval_top_list", type=str, default=None,
            help="list of maximal valence 'n1 n2 ...' (str)")
    args, unknown = parser.parse_known_args()

    with h5py.File('GParam.h5', 'a') as f:
        if args.iembeddiag is not None:
            f['/'].attrs['giembeddiag'] = args.iembeddiag
        if args.dc_mode is not None:
            f['/'].attrs['dc_mode'] = args.dc_mode
        if args.dc_nelf_list is not None:
            ne_list = args.dc_nelf_list.split()
            if len(ne_list) > 0:
                ne_list = list(map(float, ne_list))
                ne_list = numpy.array(ne_list).reshape(-1, 2)
                f['/'].attrs['dc_nelf_list'] = ne_list
                # get unique_nelf_list
                unique_nelf_list = []
                imap_list = f['/'].attrs['imap_list']
                for i, imap in enumerate(imap_list):
                    if i == imap:
                        unique_nelf_list.append(ne_list[i, :].tolist())
                args.dc_nelf_list = unique_nelf_list
            else:
                del f['/'].attrs['dc_nelf_list']

        if args.nval_bot_list is not None:
            f['/'].attrs['nval_bot_list'] = [int(s) for s in
                    args.nval_bot_lis.split()]
        if args.nval_top_list is not None:
            stmp = sys.argv[sys.argv.index('-nval_top_list')+1]
            f['/'].attrs['nval_top_list'] = [int(s) for s in
                    stmp.split()]

    # change GMOTT.h5 if necessary
    if os.path.isfile('GMott.h5'):
        with h5py.File('GMott.h5', 'a') as f:
            if args.iembeddiag is not None:
                f['/'].attrs['giembeddiag'] = args.iembeddiag
    # change settings in 'ginit.h5', re-initialize if becessary
    if args.unique_j_ev is not None or args.unique_u_ev is not None:
        re_init = True
    else:
        re_init = False

    if os.path.isfile('ginit.json'):
        data = json.load(open('ginit.json', "r"))
        if args.iembeddiag is not None:
            data['gchoice']['giembeddiag'] = args.iembeddiag
        if args.dc_mode is not None:
            data['gchoice']['dc_mode'] = args.dc_mode
        if args.dc_nelf_list is not None:
            data['gchoice']['unique_nelf_list'] = args.dc_nelf_list
        if args.u_type is not None:
            data['gchoice']['u_matrix_type'] = args.u_type
        if args.unique_j_ev is not None:
            data['gchoice']['unique_j_list'] = [float(s) for s in
                    args.unique_j_ev.split()]
            if data['gchoice']['u_matrix_type'] != 2:
                data['gchoice']['u_matrix_type'] = 1
        if args.unique_u_ev is not None:
            data['gchoice']['unique_u_list'] = [float(s) for s in
                    args.unique_u_ev.split()]
            if data['gchoice']['u_matrix_type'] != 2:
                data['gchoice']['u_matrix_type'] = 1
        with open('ginit.json', "w") as f:
            json.dump(data, f, indent = 4)
    else:
        if re_init:
            raise FileNotFoundError('file ginit.json does not exist!')

    if re_init:
        if "--new-sab" not in sys.argv:
            # preserve the previous symmetry-adapted basis.
            sys.argv.append("--fixsab")
        s_init = screen_init()
        s_init.run()
