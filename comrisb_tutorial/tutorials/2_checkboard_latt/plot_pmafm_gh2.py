from scan_checkboard import get_scan_data
import matplotlib.pyplot as plt

def plot_scan_u_pmafm():
    '''plot a figure, which compares the PM and AFM results from Gutzwiller
    calculations.
    '''
    # Get data from precalcuated PM results.
    u_list_pm, e_list_pm, z_list_pm, d_list_pm, m_list_pm = \
            get_scan_data(fname='result_pm')

    # Get data from precalcuated Hartree-Fock PM results.
    u_list_rhf, e_list_rhf, z_list_rhf, d_list_rhf, m_list_rhf = \
            get_scan_data(fname='result_pm_rhf')

    # Get data from precalcuated AFM results.
    u_list_afm, e_list_afm, z_list_afm, d_list_afm, m_list_afm = \
            get_scan_data(fname='result_afm')

    # Get data from precalcuated Hartree-Fock AFM results.
    u_list_uhf, e_list_uhf, z_list_uhf, d_list_uhf, m_list_uhf = \
            get_scan_data(fname='result_afm_uhf')

    f, axarr = plt.subplots(2, 2, sharex=True,
            gridspec_kw={'wspace':0.05, 'hspace':0.05})
    axarr[0, 0].plot(u_list_pm[::2], e_list_pm[::2], 'o', label='PM-G')
    axarr[0, 0].plot(u_list_rhf, e_list_rhf, '-', label='PM-HF')
    axarr[0, 0].plot(u_list_afm[::2], e_list_afm[::2], 'o', label='AFM-G')
    axarr[0, 0].plot(u_list_uhf, e_list_uhf, '-', label='AFM-HF')
    axarr[0, 0].legend()
    axarr[0, 0].set_ylabel('E')
    axarr[0, 0].text(0.85, 0.7, "(a)", transform=axarr[0, 0].transAxes)
    axarr[1, 0].plot(u_list_pm[::2], d_list_pm[::2], 'o')
    axarr[1, 0].plot(u_list_rhf, d_list_rhf, '-')
    axarr[1, 0].plot(u_list_afm[::2], d_list_afm[::2], 'o')
    axarr[1, 0].plot(u_list_uhf, d_list_uhf, '-')
    axarr[1, 0].set_ylabel('$d$')
    axarr[1, 0].text(0.85, 0.7, "(c)", transform=axarr[1, 0].transAxes)
    axarr[0, 1].plot(u_list_pm[::2], z_list_pm[::2], 'o')
    axarr[0, 1].plot(u_list_rhf, z_list_rhf, '-')
    axarr[0, 1].plot(u_list_afm[::2], z_list_afm[::2], 'o')
    axarr[0, 1].plot(u_list_uhf, z_list_uhf, '-')
    axarr[0, 1].yaxis.tick_right()
    axarr[0, 1].yaxis.set_label_position("right")
    axarr[0, 1].set_ylabel('Z')
    axarr[0, 1].text(0.85, 0.7, "(b)", transform=axarr[0, 1].transAxes)
    # lattice
    axins = axarr[0, 1].inset_axes([0.0, 0.05, 0.5, 0.5])
    axins.set_aspect(1)
    axins.set_xlim(0, 1)
    axins.set_ylim(0, 1)
    from matplotlib.markers import MarkerStyle as ms
    axins.scatter(0.15, 0.15, s=150, marker=ms('o', 'none'), c="r")
    axins.scatter(0.15, 0.85, s=150, marker=ms('o', 'none'), c="r")
    axins.scatter(0.85, 0.15, s=150, marker=ms('o', 'none'), c="r")
    axins.scatter(0.85, 0.85, s=150, marker=ms('o', 'none'), c="r")
    axins.scatter(0.5, 0.5, s=150, marker=ms('o', 'none'), c="r")
    axins.plot([0.15, 0.15], [0.15, 0.85], color='k', linestyle='--')
    axins.plot([0.85, 0.85], [0.15, 0.85], color='k', linestyle='--')
    axins.plot([0.15, 0.85], [0.15, 0.15], color='k', linestyle='--')
    axins.plot([0.15, 0.85], [0.85, 0.85], color='k', linestyle='--')
    axins.text(0.65, 0.3, "$t$", transform=axins.transAxes)
    axins.text(0.5, 0.6, "$U$", transform=axins.transAxes)
    axins.set_axis_off()
    import matplotlib.patches as patches
    style = "Simple, tail_width=0.5, head_width=4, head_length=8"
    kw = dict(arrowstyle=style, color="k")
    a3 = patches.FancyArrowPatch((0.85, 0.15), (0.5, 0.5),
            connectionstyle="arc3,rad=.5", **kw)
    axins.add_patch(a3)


    axarr[1, 1].plot(u_list_pm[::2], m_list_pm[::2], 'o')
    axarr[1, 1].plot(u_list_rhf, m_list_rhf, '-')
    axarr[1, 1].plot(u_list_afm[::2], m_list_afm[::2], 'o')
    axarr[1, 1].plot(u_list_uhf, m_list_uhf, '-')
    axarr[1, 1].set_ylabel(r'$2 \times <S_z>$')
    axarr[1, 1].set_ylim(-1,1)
    axarr[1, 1].yaxis.tick_right()
    axarr[1, 1].yaxis.set_label_position("right")
    axarr[1, 1].text(0.85, 0.7, "(d)", transform=axarr[1, 1].transAxes)
    axarr[1, 0].set_xlabel('U')
    axarr[1, 1].set_xlabel('U')
    axarr[1, 0].set_xlim(min(u_list_pm), max(u_list_pm))
    plt.tight_layout()
    plt.show()
    f.savefig('result_pmafm_gh2.pdf')



if __name__=="__main__":
    plot_scan_u_pmafm()
