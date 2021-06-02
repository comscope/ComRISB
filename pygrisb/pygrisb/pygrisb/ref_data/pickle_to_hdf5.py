import pickle, h5py

with h5py.File("fsoc_training_data.h5", "w") as f:
    gamma = pickle.load(open("gammaVal.pic", "rb"))
    f["/gamma"] = gamma
    print(gamma)
    for i in range(1,7):
        vtraining = pickle.load(open(f"Vtrain_n{i}.pic", "rb"))
        f[f"/val_{i}/v_training"] = vtraining
        weights = pickle.load(open(f"weights_laplace_n{i}.pic", "rb"))
        f[f"/val_{i}/weights_laplace"] = weights
        weights = pickle.load(open(f"weights_rbf_n{i}.pic", "rb"))
        f[f"/val_{i}/weights_rbf"] = weights
