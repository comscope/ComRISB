Checkboard lattice
------------------

In this example, we study 2-dimensional checkerboard lattice
with nearest neighbor hopping.
You will learn:

* How to set up a finite dimension (here 2D) lattice with multiple sites 
  per unit cell.
* How to setup simple antiferromagnetic (AFM) calculations.
* How the Gutzwiller method improves the descriptions
  at Hartree-Fock mean field level.

Time for this part: ~50 minutes.

GRISB calculation in PM state
=============================

In the model, we use nearest neighbor hopping :math:`t` as the energy unit. 
A function to setup the model for *CyGutz* calculation has been defined,

.. _gutz_model_setup:

.. autofunction:: pygrisb.model.checkboard.gutz_model_setup
    
For convenience, we have also predefined a function 
to run *CyGutz* for a list of Hubbard :math:`U` 
in paramagnetic (PM) or AFM phase.

.. autofunction:: pygrisb.model.test.scan_checkboard.generate_data

Let us perform a series of *CyGutz* calculations with increasing :math:`U`
in the PM phase, 
and check the behavior of the total energy, double occupancy, 
and quasi-particle weight (Z),
in a similar way as in the example of Bethe lattice.
Recall that in Gutzwiller-RISB theory,
Z=0 implies the system is in the Mott insulating phase,
where all the spectral weight becomes non-coherent. 

A script for a job of scanning :math:`U` in PM phase is defined as

.. autofunction:: pygrisb.model.test.scan_checkboard.scan_u

For a hands-on practice, from the directory 
*comrisb_tutorial/tutorials/checkboard_latt* 
type the following command::

    $ mkdir -p work && cd work
    $ python ../scan_checkboard.py

It will automatically generate the following results:

.. image:: _images/checkboard_pm.png
    :alt: e/d/z/m vs u
    :scale: 100 %
    :align: center

One can see similar generic feature of the local physical quantities
as in the previous Bethe lattice calculations, 
albeit a larger :math:`U_{c}` ~ 13 for the metal-insulator transition.

GRISB calculation in AFM state
==============================

Now let us perform a similar scan of :math:`U` for the checkerboard AFM phase,
i.e., there are two sublattices with up and down local spin moment.

Type the same command with an inline argument of '-sp' to trigger 
the setup of AFM phase::

    $ python ../scan_checkboard.py -sp

Here the user is advised to review the previous function gutz_model_setup_.

where spindeg changes to ``False`` by the inline argument '-sp',
and 'updn_full_list=[1,-1]' will be used in the ``batch_initialize``.
The following figure will be automatically generated:

.. image:: _images/checkboard_afm.png
    :alt: e/d/z/m vs u -- AFM
    :scale: 100 %
    :align: center

One can see that the local spin moment :math:`<S_{z}>` rises from zero 
with increasing :math:`U`, and saturates at 1, 
as a half-filling single-orbital system.
The kinetic energy renormalization factor :math:`Z` 
shows a non-monotonic behavior -- result of a competition 
between the effect of increasing onsite interaction, 
which will reduce :math:`Z`, and that of approaching empty/full filling 
of the single spin-orbital channels,
evidenced by the rising local spin moment, 
as it is defined as :math:`<S_{z}> = <n_{up}> - <n_{dn}>`.

To compare the results of PM and AFM phase, type::

    $ python ../plot_pmafm.py

One would get the following figure:

.. image:: _images/checkboard_pmafm.png
    :alt: e/d/z/m vs u -- PM & AFM
    :scale: 100 %
    :align: center

It is clear that the magnetically ordered AFM phase has lower energy 
than the PM phase. As compared with PM solutions, the kinetic energy
renormalization factor :math:`Z` only changes marginally in the AFM phase.
The local double occupancy decreases smoothly in a nonlinear way 
in AFM phase, in contrast to the linear decreasing in PM phase.

Hartree-Fock calculation in AFM state
=====================================

One would be interesting to see how the AFM solution 
obtained in the Gutzwiller method 
compare with the results from simple Hartree-Fock (HF) mean-field approach.
The HF-type calculations can be conveniently be performed by typing::

    $ python ../scan_checkboard.py -uhf

where the inline argument '-uhf' changes the flag of `iembeddiag` to 10,
corresponding to solve the embedding Hamiltonian of the Gutzwiller method 
at HF level. 
Upon execution of the above command, 
the following figure will be generated:

.. image:: _images/checkboard_afm_uhf.png
    :alt: e/d/z/m vs u -- AFM-HF
    :scale: 100 %
    :align: center

One can see similar features as the results based on Gutzwiller calculations,
but here the kinetic energy renormalization factor :math:`Z` is a constant 1.

To compare the three calculations results together, type::

    $ python ../plot_pmafm_gh.py

The following figure will be produced:

.. image:: _images/checkboard_pmafm_gh.png
    :alt: e/d/z/m vs u -- PM, AFM (G, HF)
    :scale: 100 %
    :align: center

As one can see, the Gutzwiller AFM solutions show much smoother variation 
of the local double occupancy, reduced local spin moment and smaller kinetic
energy renormalization factor :math:`Z`, and lower in total energy,
in comparison with HF results.
It is consistent with the respect of variational principle, 
because the Gutzwiller approach provides more variational degrees of freedom
by local Gutzwiller correlator. 

Hartree-Fock calculation in PM state
====================================

One might be actually surprised for seeing that the spin-unrestricted HF (UHF)
results are not bad, although HF method is known to not able to 
describe strong electron correlation effects. 
The reason is that here we are looking at the symmetry breaking solutions,
where the ``static correlation`` due to the presence of (almost) 
degenerate ground states has been removed. 
In order yo visualize the ``static correlcation`` error in HF method,
one can check the HF PM solutions. Typing the following commands::

    $ python ../scan_checkboard.py -rhf

The following results will be shown:

.. image:: _images/checkboard_pm_rhf.png
    :alt: e/d/z/m vs u -- PM, AFM (G, HF)
    :scale: 100 %
    :align: center

One can clearly see the constant double occupancy 
and the consequent linear increasing of the total energy 
in the spin-restricted HF (RHF) calculations.

Final Comparison
================

Finally, for comparisons with all calculations, type::

    $ python ../plot_pmafm_gh2.py

The following produced figure clearly 
shows the static correlation error of RHF, 
as the RHF total energy deviates more and more with increasing :math:`U`.

.. image:: _images/checkboard_pmafm_gh2.png
    :alt: e/d/z/m vs u -- PM, AFM (G, HF)
    :scale: 100 %
    :align: center

This concludes the tutorals for model calculations.
