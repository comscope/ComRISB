import os, sys, h5py
import numpy as np
import matplotlib.pyplot as plt
from pygrisb.model import checkboard2


def generate_data(u_list, spindeg=True, fname='result', iembeddiag=-1):
    '''run *CyGutz* calculations for a list of U.

    Parameters:

    * u_list: real array
      list of Hubbard U parameters.
    * spindeg: boolean number
      whether to keep spin degeneracy or not.
    * fname: string
      file name to store the results.
    * iembeddiag: integer
      flag for method to solve the embedding Hamiltonian.

      * -3: valence truncation ED with S=0 (spin-singlet) constraint;
      * -1: valence truncation ED;
      * 10: Hartree-Fock.

    Result:

    save the list of energies (e_list), double occupancy (d_list),
    quasi-particle weight (z_list), magnetic moment (m_list),
    to ``fname.dat`` text file as well as the ``fname.h5`` hdf5 file.
    '''
    # set Hubbard U=0
    u = 0.

    # remove pre-existing Gutzwiller setup files.
    for f in ['GParam.h5']:
        if os.path.exists(f):
            os.remove(f)

    # generate input files with updated with u=0
    checkboard2.gutz_model_setup(u=u, spindeg=spindeg, iembeddiag=iembeddiag)

    # the import below need GParam.h5.
    from pygrisb.run.cygutz import run_cygutz
    from pygrisb.gsolver.solve_hembed import solve_hembed

    # total energy list
    e_list = []

    # quasi-particle weight z list
    z_list = []

    # double occupancy d list
    d_list = []

    # magnetic moment
    m_list = []

    # loop over the list of U.
    for u in u_list:

        print(f' working on u = {u}')

        # modify the local Coulomb matrix of two sites
        with h5py.File('GParam.h5', 'a') as f:
            # note the transposition, which is transformation
            # from Fortran convention to c-convention.
            v2e = f['/impurity_0/V2E'][()].T
            # now update the Coulom matrix
            v2e[0,0,0,0] = v2e[0,0,1,1] = v2e[1,1,0,0] = v2e[1,1,1,1] = u
            v2e[2,2,2,2] = v2e[2,2,3,3] = v2e[3,3,2,2] = v2e[3,3,3,3] = u
            f['/impurity_0/V2E'][()] = v2e.T
            f['/'].attrs['dc_u_avg'] = [u]

        # perform the *CyGutz* calculation.
        run_cygutz(cmdlargs=False)

        # get total energy
        with h5py.File('GLog.h5', 'r') as f:
            e = f['./'].attrs["etot_model"]
            e_list.append(e.real)

            # get Z = R^\dagger R
            r = f['/impurity_0/R'][0, 0]
            z = r*r.conj()
            z_list.append(z.real)

        # if unrestricted HF calculations--
        if iembeddiag == 10:
            # in the single-orbital case, double occupancy is simply
            # <n_up>_0 <n_dn>_0
            with h5py.File('GLog.h5', 'r') as f:
                nup = f['/impurity_0/NKS'][0,0]
                ndn = f['/impurity_0/NKS'][1,1]
                d = nup*ndn
                d_list.append(d.real)
                m = nup - ndn
                m_list.append(m)
        else:
            # To get double occupancy (of impurity 1), <n_up n_dn>_G,
            # we run analysis code *exe_spci_analysis*
            solve_hembed(edinit=1,
                    analysis=True)

            # double occupancy is simply the local many-body density
            # matrix element in the valence=2 block.
            with h5py.File('HEmbedAnalysis_0.h5', 'r') as f:
                if '/valence_block_2/RHO' in f:
                    d = f['/valence_block_2/RHO'][0, 0]
                else:
                    d = 0.
                d_list.append(d.real)

            with h5py.File('HEmbed.h5', 'r') as f:
                # get the density matrix of the embedding Hamiltonian
                # note the first half of the orbotals correspond to
                # the physical one-body space
                nup = f['/impurity_0/ans/DM'][0, 0]
                ndn = f['/impurity_0/ans/DM'][1, 1]
                m = nup - ndn
                m_list.append(m.real)
        # remove z=1 initial guess
        if abs(u) < 1 or (abs(u) < 6 and not spindeg):
            os.remove("GLog.h5")

    # save to text file.
    with open(fname+'.dat', 'w') as f:
        for u, e, z, d, m in zip(u_list, e_list, z_list, d_list, m_list):
            f.write('{:.2f} {:.5f} {:.3f} {:.3f} {:.3f}\n'.format( \
                    u, e, z, d, m))

    # save to hdf5 file.
    with h5py.File(fname+'.h5', 'w') as f:
        f['/u_list'] = u_list
        f['/e_list'] = e_list
        f['/z_list'] = z_list
        f['/d_list'] = d_list
        f['/m_list'] = m_list


def scan_u(spindeg=True, iembeddiag=-1, fname='result'):
    '''run *CyGutz* calculations for a list of U.

    Result:

    it will generate results for a u_list of np.arange(0.0, 5.1, 0.2).
    '''
    if os.path.isfile(fname+'.h5'):
        return

    # set range of Hubbard U.
    u_list = np.arange(0.0, 20.0, 0.5)
    generate_data(u_list, spindeg=spindeg, fname=fname, iembeddiag=iembeddiag)


def get_scan_data(fname='result'):
    with h5py.File(fname+'.h5', 'r') as f:
        u_list = f['/u_list'][()]
        e_list = f['/e_list'][()]
        z_list = f['/z_list'][()]
        d_list = f['/d_list'][()]
        m_list = f['/m_list'][()]

    return u_list, e_list, z_list, d_list, m_list


def plot_scan_u(fname='result'):
    u_list, e_list, z_list, d_list, m_list = get_scan_data(fname=fname)

    f, axarr = plt.subplots(2, 2, sharex=True)
    axarr[0, 0].plot(u_list, e_list)
    axarr[0, 0].set_ylabel('energy')
    axarr[1, 0].plot(u_list, d_list)
    axarr[1, 0].set_ylabel('double occupancy')
    axarr[1, 0].set_ylim(-0.01, 0.26)
    axarr[0, 1].plot(u_list, z_list)
    axarr[0, 1].yaxis.tick_right()
    axarr[0, 1].yaxis.set_label_position("right")
    axarr[0, 1].set_ylabel('Z')
    axarr[0, 1].set_ylim(-0.05, 1.05)
    axarr[1, 1].plot(u_list, m_list)
    axarr[1, 1].set_ylabel('local $<S_{z}>$')
    axarr[1, 1].set_ylim(-1,1)
    axarr[1, 1].yaxis.tick_right()
    axarr[1, 1].yaxis.set_label_position("right")
    axarr[1, 0].set_xlabel('U')
    axarr[1, 1].set_xlabel('U')
    axarr[1, 0].set_xlim(min(u_list), max(u_list))
    plt.tight_layout()
    plt.show()
    f.savefig(fname+'.pdf')


if __name__=='__main__':
    fname=''
    if '-sp' in sys.argv:
        spindeg = False
        fname = 'result_afm'
        iembeddiag = -1
        scan_u(spindeg=False, fname='result_afm')
    elif '-uhf' in sys.argv:
        spindeg = False
        fname = 'result_afm_uhf'
        iembeddiag = 10
    elif '-rhf' in sys.argv:
        spindeg = True
        fname = 'result_pm_rhf'
        iembeddiag = 10
    else:
        spindeg = True
        fname = 'result_pm'
        iembeddiag = -3

    scan_u(spindeg=spindeg, fname=fname, iembeddiag=iembeddiag)
    plot_scan_u(fname=fname)
