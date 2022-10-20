# Toast--

A toolbox for image reconstruction in diffuse optical tomography.

## Purpose

Toast-- is a fork of the [Toast++](https://github.com/toastpp/toastpp) project which aims to: 

 1. Provide a drop-in replacement of the Toast++ MATLAB interface, and an updated Python interface.
 2. Enhance maintainability by simplifying the codebase, and removing features not exposed by the MATLAB/Python interfaces.
 3. Reduce external dependencies and update the build system to support compilation on modern compilers.
 4. Prioritise performance on larger problems on modern multi-core architectures.

To achieve these goals, the original codebase has been extensively modified. Details can be found in the project [history](https://github.com/samuelpowell/toastmm/blob/master/HISTORY.md).

## Compatibility

The MATLAB interface is fully compatible with the Toast++ object-oriented interface, with the
exception that OpenGL / Mesas based projection of images to a mesh is not supported. In
certain circumstances Toast-- may employ different solvers that Toast++, and hence results
may not be exactly equiavelent. Any significant differences, or breakage of existing code, will
be treated as a bug.

The Python interface provides the same functionality as the PyToast interface, though the API
differs in places. This interface is under development and should not be considered stable.

## Installation

Pre-built binaries are provided for a number of platforms.

### MATLAB

The Matlab interface is built using MATLAB R2020a on Windows x64, Linux x64, and MacOS.
It is expected that the interface will work on newer versions (please report if this is not
the case), but compilation from source will be required for older installations. Versions
prior to R2018a are not supported (owing to changes in the complex numeric API).

To install, download a release for the appropriate platform and unzip the folder in a
convenient location. Start MATLAB, change to the correct directory, then run:

```
>> toast_setup
```

This script will configure the appropriate include paths. Save these paths to make the
installation permanent. Check your installation by running one of the included demos:

```
>> demo toolbox toast
```

### Python

The Python interface is built for Python >= v3.6.0 on Windows x64, Linux x64, and MacOS.

To install, download a release for the appropriate platform an unzip the foler in a 
convenient location. Each release includes both a Python wheel and some example scripts. It is recommended that you install the python module in a suitable
virtual environment:

```
> python -m pip install toastmm-<ver>-cp36-abi3-<platform>.whl
```

The Python interface depends upon numpy and scipy. In order to execute the inlcuded examples, you will also need to install matplotlib in your environment. Check your installation by running one of the included examples:

```
> cd examples
> python recon1.py
```

## Building from source

Toast-- uses the CMake build system, and requires only a modern C++ compiler.

## Library

To build the distribution, unzip the source distribution to a convenient location,  change to this directory, and then:

```
mkdir build
cd build
cmake ../
cmake --build . --target install --config Release
```

The result of the buld process will be a set of utility executables and static library artifacts located under the `dist` subdirectory of the source tree.

### MATLAB

If a suitable version of MATLAB is located on the sytem, the MEX interface will be built automatically after compilation fo the library. The output will be located in `dist/matlab` and can be used in-place if desired.

### Python

Building the install target of the library prepares the Python build, but the actual compialtion is driven by Python. After executing the library build, change to the root directory of the distribution and build the wheel:

```
python -m pip wheel .
```

The result will be a wheel file called `toastmm-<ver>-cp36-abi3-<platform>.whl` in the root directory. The additional files included with the CI build distribution are located in `script/python`.






