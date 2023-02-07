import numpy, os, h5py
from pygrisb.model.almcubic import gutz_model_setup
from pygrisb.run.cygutz import driver_cygutz, get_cygtail


z_list = []
t_list = numpy.arange(-0.1, -0.00, 1.01)
u = 0.00

with open("uzv_1.txt", "w") as frec:
    for t in t_list:
        gutz_model_setup(u=u, v=.1, t=t, dtype=numpy.float)
        driver_cygutz(path=os.environ['WIEN_GUTZ_ROOT'],
                cygtail=get_cygtail())
        with h5py.File("GLog.h5", "r") as f:
            r = f["/impurity_0/R"][::2,::2].T
        z = r.T.conj().dot(r)
        w, v = numpy.linalg.eigh(z)
        frec.write(f"{t}  {w[0]} {v[:,0]}\n")
