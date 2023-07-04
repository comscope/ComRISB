import os, h5py, argparse
import numpy as np
from pygrisb.model import semicir as bethlatt
from builtins import zip
import matplotlib.pyplot as plt


def generate_data(u_list, mu_list, adiabatic=True, postfix="u"):
    '''run *CyGutz* calculations for a list of U and mu.
    save the results in results.

    Parameters:

    * u_list: real array
      list of Hubbard U parameters.
    * mu_list: real array
      list of chemical potential

    Result:

    save the list of energies (e_list), double occupancy (d_list),
    quasi-particle weight (z_list) and occupation (n_list)
    to ``result.dat`` text file as well as the ``result.h5`` hdf5 file.
    '''
    # set chemical potential to 0.
    mu=0

    # set Hubbard U=0
    u = 0.

    # generate input files with updated with u=0, mu=0
    bethlatt.gutz_model_setup(u=u, nmesh=5000, mu=mu, iembeddiag=-3)

    # the import below need GParam.h5.
    from pygrisb.run.cygutz import run_cygutz
    from pygrisb.gsolver.solve_hembed import solve_hembed

    # total energy list
    e_list = []

    # quasi-particle weight z list
    z_list = []

    # double occupancy d list
    d_list = []

    # occupation list
    n_list = []

    # loop over the list of U.
    for u, mu in zip(u_list, mu_list):

        # message
        print(f' working on u = {u:.3f} mu = {mu:.3f}')

        if not adiabatic and os.path.isfile("GLog.h5"):
            os.remove("GLog.h5")

        # modify the local Coulomb matrix
        with h5py.File('GParam.h5', 'a') as f:

            # note the transposition, which is transformation
            # from Fortran convention to c-convention.
            # dataset with a name in upper cases indicates fortran-order,
            # and lower cases for c-order.
            # Coulomb matrix
            v2e = f['/impurity_0/V2E'][()].T
            # additional potential for the purpose of
            # particle-hole symmetry here.
            vext = f["/vext/impurity_0/v"][()]

            # now update the Coulom matrix
            v2e[0,0,0,0] = v2e[0,0,1,1] = v2e[1,1,0,0] = v2e[1,1,1,1] = u
            f['/impurity_0/V2E'][()] = v2e.T

            # update vext, which keeps the particle-hole symmetry of the model
            # for finite U.
            vext[0, 0] = vext[1, 1] = -u/2
            f["/vext/impurity_0/v"][()] = vext

        with h5py.File("GBareH.h5", "a") as f:
            f["/"].attrs["chempot"] = mu
            f["/"].attrs["shft_init_la1"] = mu

        # perform the *CyGutz* calculation.
        run_cygutz(cmdlargs=False)

        # get total energy
        with h5py.File('GLog.h5', 'r') as f:
            # total electron filling.
            nks = f["/impurity_0/NKS"][()]
            n = nks[0, 0] + nks[1, 1]
            n_list.append(n.real)

            # total energy (subtract vext contribution for symmetry.)
            e = f['./'].attrs["etot_model"] + u/2
            e_list.append(e.real)

            # get Z = R^\dagger R
            r = f['/impurity_0/R'][0, 0]
            z = r*r.conj()
            z_list.append(z.real)

        # To get double occupancy (of impurity 1), <n_up n_dn>_G,
        # we run analysis code *exe_spci_analysis*
        solve_hembed(
                edinit=1,
                analysis=True,
                )

        # double occupancy is simply the local many-body density matrix element
        # in the valence=2 block.
        with h5py.File('HEmbedAnalysis_0.h5', 'r') as f:
            if '/valence_block_2/RHO' in f:
                d = f['/valence_block_2/RHO'][0, 0]
            else:
                d = 0.
            d_list.append(d.real)

    with open('result.dat', 'w') as f:
        for u, mu, e, z, d, n in zip(u_list, mu_list,
                e_list, z_list, d_list, n_list):
            f.write('{:.2f} {:.2f} {:.5f} {:.3f} {:.3f} {:.3f}\n'.format( \
                    u, mu, e, z, d, n))

    with h5py.File(f'result_{postfix}.h5', 'w') as f:
        f['/mu_list'] = mu_list
        f['/u_list'] = u_list
        f['/e_list'] = e_list
        f['/z_list'] = z_list
        f['/d_list'] = d_list
        f['/n_list'] = n_list


def scan_u(mu=0.0):
    '''run *CyGutz* calculations for a list of U for a given mu.

    Parameters:

    * mu: real number
      the fixed chemical potential.

    Result:

    it will generate results for a u_list of np.arange(0.0, 5.1, 0.2)
    at fixed mu.
    '''
    if os.path.isfile('result_u.h5'):
        return

    # set range of Hubbard U.
    u_list = np.arange(0.0, 5.1, 0.2)
    mu_list = [mu for u in u_list]
    generate_data(u_list, mu_list)


def scan_mu(u=5.0):
    '''run *CyGutz* calculations for a list of mu for a given u.

    Parameters:

    * u: real number
      the fixed Hubbard U

    Result:

    it will generate results for a mu_list = np.arange(0.0, 3.1, 0.1)
    at fixed u.
    '''
    if os.path.isfile('result_mu.h5'):
        return

    # set range of chemical potential mu.
    mu_list = np.arange(0.0, 3.1, 0.2)
    u_list = [u for mu in mu_list]
    generate_data(u_list, mu_list, adiabatic=False, postfix="mu")


def plot_scan_u(skipshow=False):
    with h5py.File('result_u.h5', 'r') as f:
        u_list = f['/u_list'][()]
        e_list = f['/e_list'][()]
        z_list = f['/z_list'][()]
        d_list = f['/d_list'][()]

    f, axarr = plt.subplots(3, sharex=True, figsize=(6,6))
    axarr[0].plot(u_list, e_list)
    axarr[0].set_ylabel('E')
    axarr[0].axvline(x=3.4, ls=":")
    axarr[0].text(0.85, 0.7, "(a)", transform=axarr[0].transAxes)
    axarr[1].plot(u_list, d_list)
    axarr[1].set_ylabel(r'$d$')
    axarr[1].axvline(x=3.4, ls=":")
    axarr[1].text(0.85, 0.7, "(b)", transform=axarr[1].transAxes)
    axarr[2].plot(u_list, z_list)
    axarr[2].set_ylabel('Z')
    axarr[2].set_xlabel('U')
    axarr[2].set_xlim(min(u_list), max(u_list))
    axarr[2].axvline(x=3.4, ls=":")
    axarr[2].text(0.85, 0.7, "(c)", transform=axarr[2].transAxes)
    plt.tight_layout()
    if not skipshow:
        plt.show()
    f.savefig('result_u.png')


def plot_scan_mu(skipshow=False):
    with h5py.File('result_mu.h5', 'r') as f:
        mu_list = f['/mu_list'][()]
        e_list = f['/e_list'][()]
        z_list = f['/z_list'][()]
        d_list = f['/d_list'][()]
        n_list = f['/n_list'][()]

    f, axarr = plt.subplots(4, sharex=True, figsize=(6,6))
    axarr[0].plot(mu_list, e_list)
    axarr[0].set_ylabel('E')
    axarr[0].axvline(x=1.4, ls=":")
    axarr[0].text(0.05, 0.7, "(a)", transform=axarr[0].transAxes)
    axarr[1].plot(mu_list, d_list)
    axarr[1].set_ylabel(r'$d$')
    # axarr[1].set_ylim(0, 0.25)
    axarr[1].axvline(x=1.4, ls=":")
    axarr[1].text(0.05, 0.7, "(b)", transform=axarr[1].transAxes)
    axarr[2].plot(mu_list, z_list)
    axarr[2].set_ylabel('Z')
    axarr[2].axvline(x=1.4, ls=":")
    axarr[2].set_ylim(-0.05, 1.05)
    axarr[2].text(0.05, 0.7, "(c)", transform=axarr[2].transAxes)
    axarr[3].plot(mu_list, n_list)
    axarr[3].set_ylabel('n')
    axarr[3].set_xlabel('$\mu$')
    axarr[3].axvline(x=1.4, ls=":")
    axarr[3].set_xlim(min(mu_list), max(mu_list))
    axarr[3].text(0.05, 0.7, "(d)", transform=axarr[3].transAxes)
    plt.tight_layout()
    if not skipshow:
        plt.show()
    f.savefig('result_mu.png')



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--mu", action="store_true",
            help="job switch to scan_mu.")
    parser.add_argument("--skipshow", action="store_true",
            help="do not display figure..")
    args = parser.parse_args()

    if args.mu:
        scan_mu()
        plot_scan_mu(skipshow=args.skipshow)
    else:
        scan_u()
        plot_scan_u(skipshow=args.skipshow)
