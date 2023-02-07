import h5py, numpy


class integer_basis:
    '''basis set identified by integers.
    '''
    def __init__(self):
        raise NotImplementedError("integer_basis to be inherited!")

    def h5save(self, path="/", f=None, fname="basis.h5"):
        if f is None:
            f = h5py.File(fname, "w")
        f[f"{path}/dim"] = [len(self.v)]
        f[f"{path}/vals"] = self.v


class fock_basis(integer_basis):
    '''basis set composed of all the fock states represented by integers,
    within the valence blocks.
    '''
    def __init__(self, norb, nval_bot=0, nval_top=None):
        self._nval_bot = nval_bot
        self._nval_top = nval_top if nval_top is not None else norb
        self.norb = norb

    def calc(self):
        self.set_fock_list(self._nval_bot, self._nval_top)
        self.v = [i for group in self.fock_list for i in group]

    def set_fock_list(self, nval_bot, nval_top):
        fock_list = [[] for i in range(nval_bot, nval_top+1)]
        for i in range(2**self.norb):
            n1 = bin(i).count('1')-nval_bot
            if 0 <= n1 <= nval_top-nval_bot:
                fock_list[n1].append(i)
        self.fock_list = fock_list


class double_fock_basis(fock_basis):
    '''basis set composed of the impurity-bath fock states at half-filling.
    '''
    def calc(self):
        # set full fock basis
        self.set_fock_list(0, self.norb)
        base = 2**self.norb
        nstates = sum([len(fs)**2 for fs in self.fock_list[ \
                self._nval_bot:self._nval_top+1]])
        fock2_array = numpy.empty(nstates, dtype=numpy.int32)
        iadd = 0
        for i, fs in enumerate(self.fock_list[ \
                self._nval_bot:self._nval_top+1]):
            for fs1 in fs:
                for fs2 in self.fock_list[self.norb-i-self._nval_bot]:
                    fock2_array[iadd] = fs1 + fs2*base
                    iadd += 1
        self.v = fock2_array



if __name__ == "__main__":
    with h5py.File("fock.h5", "w") as f:
        fb = double_fock_basis(14)
        fb.h5save(f=f)
