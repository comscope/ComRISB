#!/usr/bin/env python
from pygrisb.iface.ifwannier import if_gwannier
import pickle


# get if_gwannier parameters.
with open("gwannier_params.pkl", "rb") as f:
    params = pickle.load(f)


assert "corbs_list" in params, "missing corbs_list in gwannier_params.pkl!"
#     params["corbs_list"] = control['impurity_wan']
params["delta_charge"] = params.get("delta_charge", 0.)
params["wpath"] = params.get("wpath", "../wannier")
# params["wpath"] = control['wannier_directory']
params["lpath"] = params.get("lpath", "../lattice")
params["wprefix"] = params.get("wprefix", "wannier")
params["lprefix"] = params.get("lprefix", "mdl")
#     params["lprefix"] = control['allfile']
params["lrot_list"] = params.get("lrot_list", None)
# params["lrot_list"] = lrot_list
params["iso"] = params.get("iso", 1)
#    if control['spin_orbit']:
#        params["iso"] = 2
params["ispin"] = params.get("ispin", 1)
params["ismear"] = params.get("ismear", -1)
# params["ismear"] = control['ismear']
params["delta"] = params.get("delta", 0.002)
# params["delta"] = control['delta']
params["icycle"] = params.get("icycle", 0)
# params["icycle"] = icycle
params["method"] = params.get("method", 'lda+risb')
#     params["method"] = control['method']

if_gwannier(corbs_list=params["corbs_list"],
        delta_charge=params["delta_charge"],
        wpath=params["wpath"],
        lpath=params["lpath"],
        wprefix=params["wprefix"],
        lprefix=params["lprefix"],
        # lrot_list=params["lrot_list"],
        iso=params["iso"],
        ispin=params["ispin"],
        ismear=params["ismear"],
        delta=params["delta"],
        icycle=params["icycle"],
        method=params["method"],
        )
