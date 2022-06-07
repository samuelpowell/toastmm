# Toast--

Less toast, fewer problems

## Installation

Pre-built binaries are provided for a number of platforms.

### MATLAB

The Matlab interface is built using MATLAB R2020a on Windows x64, Linux x64, and MacOS (>=10.9)
It is expected that the interface will work on newer versions (please report if this is not
the case), but compilation from source will be required for older installations. Versions
prior to R2018a are not supported (owing to changes in the complex numeric API).

To install, download a release for the appropriate platform and unzip the folder in a
convenient location. Start matlab and change to the folder, then run:

```>> toast_setup```

This script will configure the appropriate include paths. Save these paths to make the
installation permanent. Check your installation by running one of the included demos:

```>> demo toolbox toast```

### Python

The Python interface follows the Python stable API with a minimum version of Python v3.6.0, 
with builds available for MacOS (>= 10.9), Windows x64, and Linux x64 (glibc, musl).

To install, download a release for the appropriate platform an unzip the foler in a 
convenient location. It is recommended that you install the python module in a suitable
virtual environment:

```> python -m pip install PyToast-<ver>-cp36-abi3-<platform>.whl```

The Python interface depends upon numpy and scipy, and matplotlib will need to be installed
to execute some examples. Check your installation by running one of the included examples:

```> cd examples
   > python recon1.py
```

## Building from source

### Requirements

### MATLAB

### Python