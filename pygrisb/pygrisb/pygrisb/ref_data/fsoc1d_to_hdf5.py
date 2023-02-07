import pickle, h5py, numpy


alphas_DM = numpy.array([1e-12,1e-12,1e-12,1e-12,1e-12])
gammas_DM = [330.57851239669424,113.17338162064283,23.3378,400,173.1301939058]
alphas_ENE = numpy.array([1e-12,1e-10,1e-10,1e-12,1e-12])
gammas_ENE = [119.16598430985587,16.2338243161228,14.632908647345994,
        62.602717588243806,54.13360302965381]

with h5py.File("fsoc_training_data_1d.h5", "w") as f:
    F0=numpy.zeros((1,6))
    file1 = './F0/data_-0.06_-0.04_-0.71_-0.051_-0.051.h5'
    with h5py.File(file1,'r+') as h:
        cc1 = numpy.real(h['/DM'][0,0])
        cc2 = numpy.real(h['/DM'][6,6])
        ff1 = numpy.real(h['/DM'][14,14])
        ff2 = numpy.real(h['/DM'][20,20])
        cf1 = numpy.real(h['/DM'][0,14])
        cf2 = numpy.real(h['/DM'][6,20])
        ene0 = numpy.real(h['/ENE'][0])
    F0[0,0]=cc1
    F0[0,1]=cc2
    F0[0,2]=ff1
    F0[0,3]=ff2
    F0[0,4]=cf1
    F0[0,5]=cf2
    f["/f0"] = F0
    f["/ene0"] = ene0
    f["/delta0"] = [-0.06, -0.04, -0.71, -0.051, -0.051]

    print([-0.06, -0.04, -0.71, -0.051, -0.051])

    for i in range(5):
        deltas = pickle.load(open(f"./fullSetDeltas_{i+1}.pic", "rb"))
        print(deltas[0])
        print(deltas[1])
        elemavgs = pickle.load(open(f"./fullSetElemAvgs_{i+1}.pic", "rb"))
        f[f"/axis_{i}/deltas"] = deltas[:,i].reshape(-1,1)
        f[f"/axis_{i}/elemavgs"] = elemavgs
        enes = pickle.load(open(f"./fullSetENE_{i+1}.pic", "rb"))
        f[f"/axis_{i}/enes"] = enes
        f[f"/axis_{i}/alpha_gamma_dm"] = [alphas_DM[i], gammas_DM[i]]
        f[f"/axis_{i}/alpha_gamma_ene"] = [alphas_ENE[i], gammas_ENE[i]]
