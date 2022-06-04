import os
import sys
import numpy as np
from setuptools import setup, Extension
from sysconfig import get_paths
from wheel.bdist_wheel import bdist_wheel
import platform

pyinc = get_paths()['include'] 
npinc = np.get_include()

# Debug flag
debug_flag = False

# Threading flag
thread_flag = True

# We pull the libraries from the release distribution install output. On Windows
# it would be fine to get these from the build tree, but on other OS we want the
# versions with the install RPATH in order that the libraries can be found. The
# rpath of the toast extension itself has to be set accordingly.
if sys.platform == 'win32':
    lib_names = ['libstoast', 'libraster', 'libfe', 'libmath']
    compile_args = ['-DTOAST_STATIC']
    link_args = ''
    pkg_files = ''
    lib_dirs = ['build/src/libmath/Release','build/src/libfe/Release','build/src/libraster/Release','build/src/libstoast/Release']
elif sys.platform == 'linux':
    lib_names = ['pthread', 'stoast', 'raster','fe', 'math']
    compile_args = ['-std=c++11','-Wno-comment', '-DTOAST_STATIC', '-fPIC']
    link_args = ''
    pkg_files = ''
    lib_dirs = ['build/src/libmath','build/src/libfe','build/src/libraster','build/src/libstoast']
elif sys.platform == 'darwin':
    lib_names = ['stoast', 'raster', 'fe', 'math']
    compile_args = ['-DTOAST_STATIC']
    link_args = ''
    pkg_files = ''
    lib_dirs = ['build/src/libmath','build/src/libfe','build/src/libraster','build/src/libstoast']
else:
    raise Exception('Unknown platform')

inc_dirs = [pyinc, npinc, 'include','src/common','src/libmath','src/libfdot','src/libfe','src/libraster','extern/eigen-3.4.0','src/libstoast']

if debug_flag:
    compile_args.extend(['-O0','-g'])

if thread_flag:
    compile_args.extend(['-DTOAST_THREAD'])

# Try to get ABI3 compat
# Following https://github.com/joerick/python-abi3-package-sample/blob/main/setup.py
#
class bdist_wheel_abi3(bdist_wheel):
    def get_tag(self):
        python, abi, plat = super().get_tag()

        if python.startswith("cp"):
            # on CPython, our wheels are abi3 and compatible back to 3.6
            return "cp36", "abi3", plat

        return python, abi, plat



# Define the toast C extension
#
module1 = Extension('toast.toastmod',
                    include_dirs = inc_dirs,
                    library_dirs = lib_dirs,
                    libraries = lib_names,
                    extra_link_args = link_args,
                    extra_compile_args = compile_args,
                    define_macros = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
                                     ("Py_LIMITED_API", "0x03060000")],
                    py_limited_api=True,
                    sources = ['src/python/toastmodule.cc',
                               'src/common/calc_jacobian.cc',
                               'src/common/calc_gradient.cc',
                               'src/common/calc_mesh.cc'])

                         
# Define installation
setup(
    name = 'PyToast',
    version = '0.1.0',
    description = 'Python TOAST extension',
    author = 'The Toast Authors',
    url = 'http://www.toastplusplus.org',
    setup_requires=['wheel'],
    install_requires=["numpy", "scipy"],
    python_requires='>=3',
    extras_require={  
        "demo": ["matplotlib"]
    },
    ext_modules = [module1],
    packages=['toast'],
    package_dir = {'toast': 'src/python/toast'},
    cmdclass={"bdist_wheel": bdist_wheel_abi3}
)
