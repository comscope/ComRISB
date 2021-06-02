import numpy, os, h5py, shutil
from pygrisb.model.almcubic_c1 import gutz_model_setup
from pygrisb.run.cygutz import driver_cygutz, get_cygtail


z_list = []
eps_c_list = numpy.arange(6.0, 8.0, 0.2)

with open("z_epsc_1.txt", "w") as frec:
    for eps_c in eps_c_list:
        gutz_model_setup(eps_c=eps_c, dtype=numpy.float, u=0.2)
        driver_cygutz(path=os.environ['WIEN_GUTZ_ROOT'],
                cygtail=get_cygtail())
        with h5py.File("GLog.h5", "r") as f:
            r = f["/impurity_0/R"][::2,::2].T
            nelect = f["/"].attrs["nelectrons"]
        with h5py.File("GIter.h5", "r") as f:
            err = f["/v_err"][()]
        maxerr = numpy.max(numpy.abs(err))
        shutil.copyfile("GLog.h5", f"GLog_{eps_c:.2f}.h5")
        z = r.T.conj().dot(r)
        w, v = numpy.linalg.eigh(z)
        frec.write(f"{eps_c:.4f} {nelect:.4f} {w[0]:.4f}" +\
                f" {v[:,0]} {maxerr:.2e}\n")
