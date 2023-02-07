import h5py, numpy, sys, os
from pygrisb.estructure.dos import get_bands
import matplotlib.pyplot as plt


def get_greek_label(kname):
    '''return the possible greek label.
    '''
    if kname.upper() in ['LAMBDA', 'GAMMA', 'DELTA', 'SIGMA', 'THETA', \
            'XI', 'PI', 'UPSILON', 'PHI', 'PSI', 'OMEGA']:
        return r'\{}{}'.format(kname[0].upper(), kname[1:].lower())
    else:
        return kname


def get_k_info():
    '''get the k distance list, special k-point label and position for
    the band structure plot.
    '''
    with h5py.File('GBareH.h5', 'r') as f:
        kx = f['/'].attrs['kptx']
        ky = f['/'].attrs['kpty']
        kz = f['/'].attrs['kptz']
        kn = f['/'].attrs['kptname']
        # reciprocal lattice vectors ordered as colomn vectors
        br = f['/'].attrs['RECIP_PRIM_VEC'].T

    ktick_pos = []
    ktick_label = []

    for i in range(len(kn)):
        if i == 0:
            klist = [0.]
        else:
            dk = numpy.array([kx[i]-kx[i-1], ky[i]-ky[i-1], kz[i]-kz[i-1]])
            dl = dk.dot(br).dot(br.T).dot(dk)
            klist.append(klist[i-1]+numpy.sqrt(dl))
        klabel = kn[i].decode('ascii').strip()
        if klabel != '':
            ktick_pos.append(klist[i])
            ktick_label.append(get_greek_label(klabel))

    return klist, ktick_pos, ktick_label


def get_estimated_gap(nocc=None):
    '''Given the number of occupied bands and the band energy file,
    estimate the band gap.
    '''
    # valence band maximum
    evmax = -1.e10
    # conduction band minimum
    ecmin = 1.e10

    # if fermi level is not given, check the fermi level
    if nocc is None:
        with h5py.File('GLog.h5') as f:
            e_fermi = f['/'].attrs['efermi']

    with h5py.File('GBands.h5', 'r') as f:
        elist = f["/e_list"][()]

    for es in elist:
        for esk in es:
            if nocc is None:
                noccp = numpy.argwhere(esk < e_fermi)[-1] + 1
                nocc = noccp[0]
                print(' number of occupied bands = {}'.format(nocc))
            evmax = max(esk[nocc-1], evmax)
            ecmin = min(esk[nocc], ecmin)
    return max(ecmin - evmax, 0.)


def get_estimated_gaps(nocc=None):
    '''Given the number of occupied bands and the band energy file,
    estimate the direct/indirect band gap with associated k-index.
    '''
    # valence band maximum
    evmax = -1.e10
    kvmax = -1
    # conduction band minimum
    ecmin = 1.e10
    kcmin = -1
    # direct band gap value
    dgap = 1000.
    kdgap = -1

    # if fermi level is not given, check the fermi level
    if nocc is None:
        with h5py.File('GLog.h5') as f:
            e_fermi = f['/'].attrs['efermi']

    with h5py.File('GBands.h5', 'r') as f:
        elist = f["/e_list"][()]

    for es in elist:
        for ik, esk in enumerate(es):
            if nocc is None:
                noccp = numpy.argwhere(esk < e_fermi)[-1] + 1
                nocc = noccp[0]
                print(' number of occupied bands = {}'.format(nocc))

            if esk[nocc-1] > evmax:
                evmax = esk[nocc-1]
                # 0-based k-index
                kvmax = ik - 1
            if esk[nocc] < ecmin:
                ecmin = esk[nocc]
                kcmin = ik - 1
            if esk[nocc] - esk[nocc-1] < dgap:
                dgap = esk[nocc] - esk[nocc-1]
                kdgap = ik - 1
    # indirect band gap size
    idgap = max(ecmin - evmax, 0.)
    return idgap, kvmax, kcmin, dgap, kdgap


def driver_get_estimated_gaps():
    '''script to print estimated direct/indirect band gaps.
    '''
    msg = r'''
    Script to print estimated direct/indirect band gaps.

    inline options:

        -n n: occupied number of bands.
    '''
    nocc = None
    if '-h' in sys.argv:
        print(msg)
        sys.exit()
    if '-n' in sys.argv:
        nocc = int(sys.argv[sys.argv.index('-n')+1])
    with h5py.File('GBareH.h5', 'r') as f:
        kname = f['/'].attrs['kptname']
    idgap, kvmax, kcmin, dgap, kdgap = get_estimated_gaps(nocc=nocc)
    print(' Direct gap = {} with k = {} {}'.format(dgap, kdgap, \
            kname[kdgap]))
    print(' Inirect gap = {} with kv = {} {} kc = {} {}'.format(idgap, \
            kvmax, kname[kvmax], kcmin, kname[kcmin]))


def get_band_structure_data():
    # check struct file
    import glob
    struct_files = glob.glob("*struct")
    if len(struct_files) > 0:
        # wien2k+grisb job
        klist, ktick_pos, ktick_label = get_k_info()
        e_skn, psi_sksna, _ = get_bands()
        if os.path.isfile("GBareH.h5"):
            f = h5py.File("GBareH.h5", "r")
        else:
            f = h5py.File("GBands.h5", "r")
        symie = f["/"].attrs["SYMIE"]-1
        f.close()
        psi_skna = psi_sksna[:, :, symie, :, :]
    else:
        # comrisb
        from pygrisb.estructure.gwannier import get_bandstructure_data
        klist, ktick_pos, ktick_label, e_skn, psi_skna = \
                get_bandstructure_data()

    return klist, ktick_pos, ktick_label, e_skn, psi_skna


def plot_band_sturture(
        emin=-10.,
        emax=10.,
        scale_weight=20,
        ireps_plot=False,
        # plot the band structures with ireps-resolved weight.
        ):
    '''plot band structure with overall correlated orbital character
    in the assigned energy window [emin, emax] in the file `bands.pdf`.
    '''
    klist, ktick_pos, ktick_label, e_skn, psi_skna = \
            get_band_structure_data()
    e_skn = numpy.asarray(e_skn)

    # further weights for ireps.
    with h5py.File("GLog.h5", "r") as f:
        rmat = f["/"].attrs["RMAT"].swapaxes(1, 2)
    nasotot = rmat.shape[1]
    # get total correlated orbital weight.
    psi_skn_f = numpy.einsum('...ij,...ij->...i',
            psi_skna[...,:,:nasotot],
            psi_skna.conj()[...,:,:nasotot])
    psi_skn_f = psi_skn_f.real
    # multiply a scaling factor
    psi_skn_f *= scale_weight

    # save data
    with h5py.File("bands_plot.h5", "w") as f:
        f["/e_skn"] = e_skn
        f["/k_dist"] = klist
        f["/k_tick_pos"] = ktick_pos
        f["/k_tick_label"] = [s.encode("utf-8") for s in ktick_label]

    plot_bands(e_skn,
            klist,
            psi_skn_f,
            ktick_pos,
            ktick_label,
            emin,
            emax,
            fname="bands.pdf",
            )

    if ireps_plot:
        # block diagonal indices
        from pygrisb.math.matrix_util import get_block_diag_indices
        blk_ids = get_block_diag_indices(rmat[0, :nasotot, :nasotot])

        psi_skn_ireps = []
        for ista, iend in zip(blk_ids[:-1], blk_ids[1:]):
            psi_skn_ireps.append(
                    numpy.einsum('...ij,...ij->...i',
                    psi_skna[...,:,ista:iend],
                    psi_skna.conj()[...,:,ista:iend])/scale_weight
                    )

        for i, psi_skn_irep in enumerate(psi_skn_ireps):
            plot_bands(e_skn,
                    klist,
                    psi_skn_irep,
                    ktick_pos,
                    ktick_label,
                    emin,
                    emax,
                    fname=f"bands_{i}.pdf",
                    )


def plot_bands(e_skn,
        klist,
        psi_skn_f,
        ktick_pos,
        ktick_label,
        emin,
        emax,
        fname="bands.pdf",
        ):
    # start plotting
    fig, ax = plt.subplots()
    for n in range(e_skn.shape[2]):
        ax.plot(klist, e_skn[0, :, n], 'k-')
        ax.scatter(klist, e_skn[0,:,n], s = psi_skn_f[0, :, n], \
                c = 'r', edgecolors = 'r')
    if len(e_skn) == 2:
        for n in range(e_skn.shape[2]):
            ax.plot(klist, e_skn[1, :, n], 'k--')
            ax.scatter(klist, e_skn[1,:,n], s = psi_skn_f[1, :, n], \
                    c = 'b', edgecolors = 'b')

    ax.axhline(y = 0, ls = ':', lw = 2)
    # High-symmetry lines and labels
    for x1 in ktick_pos[1:-1]:
        ax.axvline(x = x1, ls = '--')
    ax.set_xticks(ktick_pos)
    ktick_label = [r"${}$".format(s) for s in ktick_label]
    ax.set_xticklabels(ktick_label)
    ax.set_ylabel("E (eV)")
    ax.set_xlim(klist[0], klist[-1])
    ax.set_ylim(emin, emax)
    plt.show()
    fig.savefig(fname)



def driver_plot_bands():
    import argparse
    parser = argparse.ArgumentParser(
            description="GRISB band structure plot.")
    parser.add_argument("-el", "--emin", type=float, default=None,
            help="band structure lower energy bound (float)")
    parser.add_argument("-eh", "--emax", type=float, default=None,
            help="band structure upper energy bound (float)")
    args = parser.parse_args()

    plot_band_sturture(emin=args.emin, emax=args.emax)



if __name__=='__main__':
    plot_band_sturture(emin=-3.3, emax=5.8)
