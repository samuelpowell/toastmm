import os
import sys
import numpy as np
from setuptools import setup, Extension
from sysconfig import get_paths

major = "%d" % sys.version_info[0]
minor = "%d" % sys.version_info[1]

pyinc = get_paths()['include'] 
npinc = np.get_include()

module1 = Extension('toast.toastmod',
                    include_dirs = [pyinc,
                                    npinc,
                                    '../..',
                                    '../include',
                                    '../src/libmath',
                                    '../src/libfe',
                                    '../src/libstoast',
                                    '../extern/eigen-3.4.0'],
                    libraries = ['libmath','libfe','libstoast'] if "nt" in os.name else ['math','fe','stoast'],
                    library_dirs = ['../build/src/libfe/Release',
                                    '../build/src/libmath/Release',
                                    '../build/src/libstoast/Release'],
                    runtime_library_dirs = None if "nt" in os.name else ['../lib'],
                    sources = ['toastmodule.cc'])

# Install library files, on windows these must go alongside the library, on Linux and
# MacOS the rpath will be set such that they are sought in the lib subdirectory
if "nt" in os.name:
    lib_files = ('lib/site-packages/toast', ['../../build/src/libfe/Release/libfe.dll',
                      '../../build/src/libmath/Release/libmath.dll',
                      '../../build/src/libstoast/Release/libstoast.dll'])
else:
    lib_files = ('lib/site-packages/toast', ['../../build/src/libfe/Release/libfe.dll',
                         '../../build/src/libmath/Release/libmath.dll',
                         '../../build/src/libstoast/Release/libstoast.dll'])
                         

setup(
    name = 'PyToast',
    version = '120529',
    description = 'Python TOAST extension',
    author = 'Martin Schweiger',
    url = 'http://www.toastplusplus.org',
    setup_requires=['wheel'],
    install_requires=["numpy", "scipy", "glumpy"],
    python_requires='>=3',
    extras_require={  
        "demo": ["matplotlib"]
    },
    ext_modules = [module1],
    packages=['toast'],
    data_files = [lib_files]
)
