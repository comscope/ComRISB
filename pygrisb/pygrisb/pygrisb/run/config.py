import os, ast


def get_options():
    if os.path.isfile("cygutz.conf"):
        conf = ast.literal_eval(open("cygutz.conf", "r").read())
        return conf
    else:
        return {}
