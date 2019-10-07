# Installation of FairShip on MacOs systems

## Premilinary requirements
FairShip inherits the main framework and many packages for the alibuild and alidist framework.
Is therefore useful to follow the first steps explained here, to ensure you have the needed prerequisites:
https://alice-doc.github.io/alice-analysis-tutorial/building/prereq-macos.html. Basically, you will need:
* Homebrew for package installation
* Similarly, pip for python package installation
* gfortran
System Integrity Protection is recommended to be disabled to avoid complains from the MacOs during the software installation

### Other packages which are required but are not listed in the alibuild guide:
* `pip install scipy sklearn`
* `brew install openssl davix`
* `brew install autconf automake texinfo gettext libtool`

## Installation with  FairShip/alibuild.sh. 
It is recommended to check that Python, Python-modules and GCC-Toolchain are found from the system, because their build usually leads to errors.
If necessary, the required version of GCC-Toolchain can be installed with brew install llvm@8 and brew install gcc@8.
Note: MacOs uses clang as internal cc and c++ compiler. It is needed to set the links "cc -> clang" and "c++ -> clang++".

## Encountered issues:
1. pythia6.sh: a line requires to copy libpythia6.so, which does not exist in Mac, while the equivalent libpythia6.dylib is already copied by default.
The line must be commented;
2. The system cannot find libtoolize. This is due to the fact that it is called glibtoolize. A symbolic link solves this issue;
3. gsl/gsl_version.h is not found, even after installing it with `brew install gsl`. The installation path needs to be included to the CPATH env variable
