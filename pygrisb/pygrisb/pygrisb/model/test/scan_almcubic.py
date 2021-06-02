import numpy, os, h5py
from pygrisb.model.almcubic import gutz_model_setup
from pygrisb.run.cygutz import driver_cygutz


z_list = []
u_list = numpy.arange(0, 0.08, 0.01)

for u in u_list:
    gutz_model_setup(u=u, v=.1)
    driver_cygutz(path=os.environ['WIEN_GUTZ_ROOT'])
    with h5py.File("GLog.h5", "r") as f:
        r = f["/impurity_0/R"][0,0]
    z = abs(r)**2
    z_list.append(z)

data = numpy.array([u_list, z_list]).T
numpy.savetxt(f"uzv_1.txt", data)
