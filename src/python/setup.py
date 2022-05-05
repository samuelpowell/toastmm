import os
import sys
import numpy as np
from setuptools import setup, Extension
from sysconfig import get_paths
import platform

pyinc = get_paths()['include'] 
npinc = np.get_include()

# Library dependency location
distroot = '../'

# We pull the libraries from the release distribution install output. On Windows
# it would be fine to get these from the build tree, but on other OS we want the
# versions with the install RPATH in order that the libraries can be found. The
# rpath of the toast extension itself has to be set accordingly.
if sys.platform == 'win32':
    lib_dir = distroot + 'bin'
    lib_ext = '.dll'
    lib_names = ['libmath','libfe','libstoast']
    lib_dest = '././toast'
    link_dir = [distroot + 'lib']
    link_args = ''
    pkg_files = {'toast': ['*.dll']}
elif sys.platform == 'linux':
    lib_dir = distroot + 'lib'
    lib_ext = '.so'
    lib_names = ['math','fe','stoast']
    lib_dest = '/toast/lib'
    link_dir = [distroot + 'lib']
    link_args = ['-Wl,-rpath=$ORIGIN/lib']
    pkg_files = {'toast': ['lib/*.so']}
elif sys.platform == 'darwin':
    lib_dir = distroot + 'lib'
    lib_ext = '.dylib'
    lib_names = ['math','fe','stoast']
    lib_dest = '/toast/lib'
    link_dir = [distroot + 'lib']
    link_args = ['-Wl,-rpath=@loader_path/lib']
    pkg_files = {'toast': ['lib/*.dylib']}
else:
    raise Exception('Unknown platform')

# Define the toast C extension
#
# TODO: To enable an isolated installation, we need to make the paths to the source 
# headers absolute. Again, a good way to do this might be to use CMake as above, and
# to configure this file.
module1 = Extension('toast.toastmod',
                    include_dirs = [pyinc,
                                    npinc,
                                    '../..',
                                    '../../include',
                                    '../../src/libmath',
                                    '../../src/libfe',
                                    '../../src/libstoast',
                                    '../../extern/eigen-3.4.0'],
                    library_dirs = link_dir,
                    libraries = lib_names,
                    extra_link_args = link_args,
                    sources = ['toastmodule.cc'])

                         
# Define installation
setup(
    name = 'PyToast',
    version = '120529',
    description = 'Python TOAST extension',
    author = 'Martin Schweiger',
    url = 'http://www.toastplusplus.org',
    setup_requires=['wheel'],
    install_requires=["numpy", "scipy"],
    python_requires='>=3',
    extras_require={  
        "demo": ["matplotlib"]
    },
    ext_modules = [module1],
    packages=['toast'],
    package_data = pkg_files,
    include_package_data=True
    #data_files = [lib_files]
)
