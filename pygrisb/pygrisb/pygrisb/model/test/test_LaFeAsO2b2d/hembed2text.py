import h5py, numpy
from scipy.linalg import block_diag


tol = 1.e-8

# read data in spin-faster basis
with h5py.File(f'HEmbed.h5', 'r') as f:
    d = f['impurity_0/D'][::2, ::2].T
    h1e = f['impurity_0/H1E'][::2, ::2].T
    lam = f['impurity_0/LAMBDA'][::2, ::2].T
    v2e = f['impurity_0/V2E'][::2, ::2, ::2, ::2].T

# one spin-block is needed.
norb = h1e.shape[0]
v1e = block_diag(h1e, -lam)
v1e[:norb, norb:] = d.T
v1e[norb:, :norb] = d.conj()

with open("v1e.dat", "w") as f:
    f.write(f"# norb: {norb} spin index NOT included \n")
    for i, vi in enumerate(v1e):
        for j, vij in enumerate(vi):
            if abs(vij) > tol:
                f.write(f"{i:2d} {j:2d} {vij:15.9f}\n")

with open("v2e.dat", "w") as f:
    f.write(f"# norb: {norb} spin index NOT included \n")
    f.write(f"# chemsit's convention: v[i,j,k,l] c_i^\dag c_k^\dag c_l c_j\n")
    for i, vi in enumerate(v2e):
        for j, vij in enumerate(vi):
            for k, vijk in enumerate(vij):
                for l, vijkl in enumerate(vijk):
                    if abs(vijkl) > tol:
                        f.write(f"{i:2d} {j:2d} {k:2d} {l:2d} {vijkl:15.9f}\n")


with open("hembed.inp", "w") as f:
    # NF
    f.write(f"{norb*4}\n")
    # Ne
    f.write(f"{norb*2}\n")
    # one-body
    for i, vi in enumerate(v1e):
        for j, vij in enumerate(vi):
            f.write(f"{vij:15.9f}\n")
            if i == j:
                break
    # two-body
    v2e2 = numpy.zeros([norb*2]*4)
    v2e2[:norb, :norb, :norb, :norb] = v2e
    for i, vi in enumerate(v2e2):
        for j, vij in enumerate(vi):
            for k, vijk in enumerate(vij):
                for l, vijkl in enumerate(vijk):
                    f.write(f"{vijkl:15.9f}\n")
            if i == j:
                break
