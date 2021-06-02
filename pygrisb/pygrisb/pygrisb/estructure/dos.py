import h5py, json, argparse
from math import pi, sqrt
import numpy as np
from pygrisb.basic import units
import matplotlib.pyplot as plt


class DOS:

    def __init__(self, w_k, e_skn,  width=0.1, window=None, npts=201):
        """Electronic Density Of States object.

        width: float
          Width of guassian smearing.
        window: tuple of two float
          Use ``window=(emin, emax)``.  If not specified, a window
          big enough to hold all the eigenvalues will be used.
        npts: int
          Number of points.

        """

        self.npts = npts
        self.width = width
        self.w_k = w_k
        self.e_skn = e_skn

        if window is None:
            emin = self.e_skn.min() - 5 * self.width
            emax = self.e_skn.max() + 5 * self.width
        else:
            emin, emax = window

        self.energies = np.linspace(emin, emax, npts)


    def delta(self, energy):
        """
        Return a delta-function centered at energy.
        """
        x = -((self.energies - energy) / self.width)**2
        return np.exp(x) / (sqrt(pi) * self.width)


    def get_energies(self):
        """"
        return energy mesh.
        """
        return self.energies


    def get_dos_component(self, psi_skn_w):
        """
        Get array of DOS values.
        """
        dos_list = []

        # spin
        for e_kn, psi_kn_w in zip(self.e_skn,psi_skn_w):
            dos = np.zeros(self.npts, dtype=np.complex)
            # k-point
            for w, e_n, psi_n_w in zip(self.w_k, e_kn, psi_kn_w):
                # band index
                for e, psi_w in zip(e_n, psi_n_w):
                    dos += w * self.delta(e) * psi_w
            dos_list.append(dos)
        return np.array(dos_list)


    def get_dos_t(self):
        """
        Get array of DOS values.
        """
        dos_list = []

        # spin
        for e_kn in self.e_skn:
            dos = np.zeros(self.npts)
            # k-point
            for w, e_n in zip(self.w_k, e_kn):
                # band index
                for e in e_n:
                    dos += w * self.delta(e)
            dos_list.append(dos)
        return np.array(dos_list)


def get_full_vec(v, n, val_default=33.):
    vp = np.ones(n, dtype=v.dtype)*val_default
    n_v = len(v)
    if n_v >= n:
        vp = v[:n]
    else:
        vp[:n_v] = v
    return vp


def get_bands(
        coherentonly=False,
        nval=None,
        ):
    '''get band anergies and the correlated orbital characters.
    '''
    with h5py.File("GBareH.h5", 'r') as f:
        # number of k-points
        nkp = f["/"].attrs["kptdim"]
        # maximal number of bands calculated over the k-points.
        nbmax = f["/"].attrs["nbmax"]
        # band index list specifying the range of bands used
        # for the construction of correlated orbitals.
        bnd_ne = []
        ispin = f["/"].attrs["ispin"]
        for i in range(ispin):
            if f"/ispin_{i}/ne_list" in f:
                bnd_ne.append(f[f"/ispin_{i}/ne_list"][()])
            else:
                bnd_ne.append(np.array([[nbmax, 0, nbmax]
                        for k in range(nkp)]))

    try:
        unit = json.load(open("ginit.json", "r"))["gchoice"]["unit"]
    except:
        unit = 'eV'
    if "ryd" in unit.lower():
        use_rydberg = True
    else:
        use_rydberg = False

    with h5py.File("GLog.h5", 'r') as f:
        # Gutzwiller fermi level
        e_fermi = f["/"].attrs["efermi"]

    with h5py.File("GBands.h5", 'r') as f:
        e_skn = f["/e_list"][()]
        if nval is not None:
            vtop = np.max(e_skn[:, :, nval])
            cbot = np.min(e_skn[:, :, nval+1])
            if cbot > vtop:
                e_fermi = (vtop + cbot)/2
        e_skn -=  e_fermi

        if use_rydberg:
            e_skn *= units.rydberg_to_ev
        rmat = f["/BND_R"][()].swapaxes(1, 2)
        nasotot = rmat.shape[1]
        # expansion coefficients of the correlated orbitals in terms of the
        # band wavefunctions, i.e., <\psi_{sks, n}|\phi_{sks, \alpha}>
        # with sks := ispin, ikpt, isym.
        psi_sksna = f["/V_LIST"][()]

    # choose bands within enrgy window only.
    nbmaxin = psi_sksna.shape[3]

    # we only plot bands in the energy window of the projector
    e_skn_in = e_skn[:, :, :nbmaxin]

    # kick the missing points away.
    for isp, e1 in enumerate(e_skn):
        ispp = min(isp, ispin-1)
        for ikp, e12 in enumerate(e1):
            n1, n2 = bnd_ne[ispp][ikp, 1], bnd_ne[ispp][ikp, 0]
            e_skn_in[isp, ikp, :] = get_full_vec(e12[n1:n2], nbmaxin)

    if coherentonly:
        for isp, a in enumerate(psi_sksna):
            psi_sksna[isp,:,:,:,:nasotot] = np.tensordot(a[:,:,:,:nasotot],
                    rmat[isp,:,:].conj(), ([3],[0]))
    return e_skn_in, psi_sksna, nasotot


def h5get_dos(
        ewin=(-3., 5.),
        delta=0.05,
        npts=1001,
        coherentonly=False,
        nval=None,
        ):
    '''
    get total dos and the total correlated orbital component.
    '''
    with h5py.File("GBareH.h5", 'r') as f:
        # list of k-point weight.
        w_k = f["/"].attrs["kptwt"]

    # get band energies and the orbital characters.
    e_skn, psi_sksna, nasotot = get_bands(
            coherentonly=coherentonly,
            nval=nval,
            )

    # get total dos
    dos = DOS(w_k, e_skn,  width=delta, window=ewin, npts=npts)
    energies = dos.get_energies()
    if not coherentonly:
        dos_t = dos.get_dos_t()
    else:
        # get coherent weight
        psi_skn_coh = np.einsum('...ijk,...ijk->...j', psi_sksna[...,:,:,:], \
                psi_sksna.conj()[...,:,:,:])/psi_sksna.shape[2]
        dos_t = dos.get_dos_component(psi_skn_coh)

    # get total correlated orbital component.
    psi_skn_f = np.einsum('...ijk,...ijk->...j', psi_sksna[...,:,:,:nasotot], \
            psi_sksna.conj()[...,:,:,:nasotot])/psi_sksna.shape[2]
    dos_f = dos.get_dos_component(psi_skn_f)

    return energies, dos_t, dos_f


def plot_dos_tf(energies, dos_t, dos_f):
    '''plot total dos and total correlated component
    and save in file `dos.pdf`.
    '''
    fig, ax = plt.subplots()
    if len(dos_t) == 1:
        ax.fill_between(
            energies, 0, dos_t[0], facecolor='grey', alpha=0.5)
        ax.plot(energies, dos_t[0], color='grey', label='tot')
        ax.plot(energies, dos_f[0], color='red', label='$corr. orb.$')
        ax.set_ylim(0.)
    else:
        ax.fill_between(
            energies, 0, dos_t[0], facecolor='grey', alpha=0.5)
        ax.fill_between(
            energies, 0, -dos_t[1], facecolor='grey', alpha=0.5)
        ax.plot(energies, dos_t[0], color='grey', label='tot-up')
        ax.plot(energies, -dos_t[1], color='grey', label='tot-dn')
        ax.plot(energies, dos_f[0], color='red', label='$corr-orb-up$')
        ax.plot(energies, -dos_f[1], color='blue', label='$corr-orb-dn$')
        ax.axhline(y=0, ls='-', color='black', linewidth=0.5)
    ax.set_xlim(np.min(energies), np.max(energies))
    ax.axvline(x=0, ls='--')
    ax.set_ylabel("DOS (states/f.u.)")
    ax.set_xlabel("E (eV)")
    plt.title("DOS with correlated orbital-component")
    plt.legend()
    fig.tight_layout()
    plt.show()
    fig.savefig('dos.pdf')


def driver_plot_dos():
    parser = argparse.ArgumentParser()
    parser.add_argument("-el", "--emin", type=float, default=-5.0,
            help="energy lower bound (float)")
    parser.add_argument("-eu", "--emax", type=float, default=5.0,
            help="energy upper bound (float)")
    parser.add_argument("-d", "--delta", type=float, default=0.05,
            help="gaussian smearing factor (float)")
    parser.add_argument("-n", "--npts", type=int, default=1001,
            help="energy mesh points (int)")
    parser.add_argument("-c", "--coherentonly", action="store_true",
            help="coherent component only (bool)")
    parser.add_argument("--nval", type=int, default=None,
            help="valence top index (int) for gap system")
    args = parser.parse_args()
    # get dos
    energies, dos_t, dos_f = h5get_dos(
            ewin=(args.emin, args.emax),
            delta=args.delta,
            npts=args.npts,
            coherentonly=args.coherentonly,
            nval=args.nval,
            )

    # pick up-component to plot
    plot_dos_tf(energies, dos_t, dos_f)
    data = np.concatenate(([energies], dos_t, dos_f)).T.real
    np.savetxt("dos.dat", data, fmt="%.6f")


if __name__ == "__main__":
    '''
    Test run.
    '''
    energies, dos_t, dos_f = h5get_dos(ewin=(-0.1, 0.1), delta=0.001)
    # pick up-component.
    plot_dos_tf(energies, dos_t, dos_f)
