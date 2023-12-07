import numpy
from scipy.special import erfc
from scipy.optimize import bisect

# get_fermi_level computes and returns the Fermi level, which is defined
#   as the energy which solves num_e = \sum get_fermi_weight.  The fermi weight
#   may be governed by a Fermi-Dirac distribution if ismear=-1, or a 
#   Gaussian distribution if ismear=0.
# bnd_es = list of band energies
# wklist = list of k-point weights
# num_e = the number of electrons in the wannier space, minus a mott correction
# delta = the inverse of the temperature
# 'ismear' = ismear # control['ismear'] This specifies the 
#       integration method: fermi or gaussian smearing.
#       gaussian smearing (ismear=0) or Fermi smearing (-1) or 
#           tetra-hedron method (-5) will be used for the brillouin zone integration.
#   brillouin zone integration method: fermi or gaussian smearing
# 'iso' = iso # if control['spin_orbit'] then iso=2, otherwise iso=1
def get_fermi_level(bnd_es, wklist, num_e, delta=0.0258, \
        ismear=0, iso=1):
    # delta = 0.0258519909eV = 300K
    emin = numpy.min(bnd_es)
    emax = numpy.max(bnd_es)
    efermi = bisect(err_fun, emin, emax, \
            args=(bnd_es, wklist, num_e, delta, ismear, iso))
    return efermi


def err_fun(mu, bnd_es, wklist, num_e, delta, ismear, iso):
    ferwes = get_fermi_weight(mu, bnd_es, wklist, delta, ismear=ismear,
            iso=iso)
    _nume = numpy.sum(ferwes)
    return _nume - num_e

# get_fermi_weight computes and returns the Fermi weight.
#   The Fermi weight is defined individually as a function of orbital (including
#       Mott orbitals and also bands), spin, and k-point.
#  The fermi weight
#   may be governed by a Fermi-Dirac distribution if ismear=-1, or a 
#   Gaussian distribution if ismear=0. 
#   For Mott orbitals, the Fermi weight is defined as ne_mott/no_mott.
# mu = the Fermi energy
# bnd_es = list of band energies
# wklist = list of k-point weights
# delta = the inverse of the temperature
# 'ismear' = ismear # control['ismear'] This specifies the 
#       integration method: fermi or gaussian smearing.
#       gaussian smearing (ismear=0) or Fermi smearing (-1) or 
#           tetra-hedron method (-5) will be used for the brillouin zone integration.
#   brillouin zone integration method: fermi or gaussian smearing
# 'iso' = iso # if control['spin_orbit'] then iso=2, otherwise iso=1   
# no_mott = sum_impurities GParam.h5/mott/impurity_{i}/num_mott_orbitals
# ne_mott = sum_impurities GParam.h5/mott/impurity_{i}/num_mott_electrons
def get_fermi_weight(mu,
        bnd_es,
        wklist,
        delta=0.0258,  # ev unit, room temp.
        ismear=0,
        iso=1,
        no_mott=0,
        ne_mott=0.,
        ):
    
    ferwes = []
    
    focc = 0 if no_mott == 0 else ne_mott/no_mott
    
    # no_mott is the number of mott orbitals, and it includes spin and orbitals
    no_mott1 = no_mott if iso == 2 else no_mott//2
    
    for bnd_e in bnd_es:
        ferwes.append([])
        for bnd_ek, wk in zip(bnd_e, wklist):
            ferwes[-1].append([])
            for e in bnd_ek:
                x = (e-mu)/delta
                if ismear == -1:
                    fw = fermi_dirac(x)
                elif ismear == 0:
                    fw = gaussian(x)
                else:
                    raise ValueError("Not defined ismear = {}!".format(ismear))
                fw *= wk
                ferwes[-1][-1].append(fw)
                
            # modify due to mott
            for i in range(no_mott1):
                ferwes[-1][-1][-1-i] = focc*wk

    ferwes = numpy.asarray(ferwes)
    if len(ferwes) == 1 and iso == 1:
        ferwes *= 2
    return ferwes


def fermi_dirac(x):
    """
    Return fermi-dirac distribution weight.
    """
    if x < -200:
        f = 1.
    elif x > 200:
        f = 0.
    else:
        f = 1./(numpy.exp(x) + 1)
    return f


def gaussian(x):
    """
    Return gaussian distribution weight.
    """
    if x < -7:
        f = 2.
    elif x > 7:
        f = 0.
    else:
        f = erfc(x)
    return f/2.0
