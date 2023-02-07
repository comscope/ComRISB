# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

# here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
# with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
#    long_description = f.read()

setup(
    name='pygrisb',
    version='2019.3.18.1',
    author="Yongxin Yao, Nicola Lanata, Xiaoyu Deng",
    author_email="cygutz@gmail.com",
    description='Python libraries for CyGutz package',
    url="https://gitlab.com/yxphysics/pygrisb",
    license='GPL',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',

        # Specify the Python versions you support here.
        'Programming Language :: Python :: 3.7',
    ],
    packages=find_packages(),
    install_requires=['numpy', 'scipy', "pymatgen", "mpi4py", \
            "matplotlib", 'h5py', 'sympy', 'inquirer', 'ase', \
            "phonopy"],
)
