#!/usr/bin/env python

import h5py, inquirer, os
import numpy as np
import pygrisb.math.matrix_basis as mb



def get_shrink_sigma(sigma):
    '''
    Shrink the structure of sigma if numbers skipped.
    '''
    elem_max = np.max(sigma)
    for elem in range(elem_max):
        if elem == np.max(sigma):
            break
        while elem + 1 not in sigma:
            spots = np.where(sigma > elem)
            sigma[spots] -= 1
    return sigma


def print_sigma(sigma, title):
    '''
    Print the structure of sigma with basis indices.
    '''
    print(title)
    print("index  " + ''.join("%4d " % (i) for i in range(len(sigma))) + '\n')
    for i, row in enumerate(sigma):
        print("%4d   " % (i) + ''.join("%4d " % (j) for j in row))


def init_mott(fname='GParam.h5'):
    '''
    Generate the input file for Mott phase calculation.
    '''
    with h5py.File("GParam.h5", "r") as f:
        val = f["/impurity_0/V2E"][0,0,0,0]
        dtype= val.dtype

    dim_hs_r_list = []
    dim_hs_l_list = []
    f = h5py.File(fname, 'a')
    num_imp = f["/"].attrs["num_imp"]
    imap_list = f["/"].attrs["imap_list"]
    if "/mott" in f:
        del f["/mott"]
    for imp in range(num_imp):
        prepath = "/impurity_" + str(imp)
        if imap_list[imp] == imp:
            sigma = f[prepath + "/symbol_matrix"][()]
            print("**********  impurity " + str(imp) + "  **********")
            print_sigma(sigma, ' symbolic matrix:')
            while True:
                orb_mott = input(
                        ' Please provide the indices of orbitals to'
                        +' be Mott localized \n (e.g., 0 2 ): ')
                yn = input(
                    ' You selected [' +
                    orb_mott +
                    '] \n to be Mott localized, right? (y/n):')
                if 'y' in yn or 'Y' in yn:
                    orb_mott = orb_mott.split()
                    ind_orb_mott = [int(s) for s in orb_mott]
                    while True:
                        ne_mott = input(
                                ' Please provide the total number of'
                                +' Mott localized electrons (per unit cell): ')
                        yn = input(
                            ' Total ' +
                            ne_mott +
                            ' electrons will be Mott localized, right? (y/n):')
                        if 'y' in yn or 'Y' in yn:
                            break
                    break

            for ind in ind_orb_mott:
                sigma[ind, :] = 0
            sigma_r = get_shrink_sigma(sigma)
            print_sigma(sigma_r, ' R structure:')
            smb = mb.dense_matrix_basis(symbol_matrix=sigma_r,
                    dtype=dtype, btype="general")
            hs_r = smb.basis
            sigma_l = np.copy(sigma_r)
            for ind in ind_orb_mott:
                sigma_l[:, ind] = 0
            sigma_l = get_shrink_sigma(sigma_l)
            print_sigma(sigma_l, ' Lambda structure:')
            smb = mb.dense_matrix_basis(symbol_matrix=sigma_l,
                    dtype=dtype, btype="herm")
            hs_l = smb.basis

        dim_hs_r_list.append(len(hs_r))
        dim_hs_l_list.append(len(hs_l))

        f[f'mott/{prepath}/symbol_matrix_r'] = sigma_r
        if len(hs_r) > 0:
            f[f'/mott/{prepath}/matrix_basis_r'] = hs_r
        f[f'/mott/{prepath}/symbol_matrix_l'] = sigma_l
        if len(hs_l) > 0:
            f[f'/mott/{prepath}/matrix_basis_l'] = hs_l
        f[f'/mott/{prepath}'].attrs['num_mott_orbitals'] = len(ind_orb_mott)
        f[f'/mott/{prepath}'].attrs['num_mott_electrons'] = int(ne_mott)
        f[f'/mott/{prepath}'].attrs['mott_orbital_indices'] = \
                np.asarray(ind_orb_mott)

    questions = [
            inquirer.List('iembeddiag',
                    message="Method to solve embedding Hamiltonian.",
                    choices=[
                            ("Valence truncation ED (VED)", -1),
                            ("VED with Sz symmetry", -2),
                            ("VED with S=0 constraint", -3),
                            ("VED with Jz symmetry", -4)],)]
    ans = inquirer.prompt(questions)
    f['/'].attrs['giembeddiag'] = int(ans["iembeddiag"])
    f[f'/mott'].attrs['matbs_r_dim_list'] = dim_hs_r_list
    f[f'/mott'].attrs['matbs_l_dim_list'] = dim_hs_l_list
    f.close()


if __name__ == "__main__":
    init_mott()
    # will have to start from scratch.
    if os.path.exist("GLog.h5"):
        os.remove("GLog.h5")
