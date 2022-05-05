[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2651100.svg)](https://doi.org/10.5281/zenodo.2651100)[![Twitter Follow](https://img.shields.io/twitter/follow/spinw4.svg?style=social&label=Follow)](https://twitter.com/intent/user?screen_name=spinw4) [![Old Github All Releases](https://img.shields.io/github/downloads/tsdev/spinw/total.svg)](https://github.com/tsdev/spinw/releases)[![Github All Releases](https://img.shields.io/github/downloads/spinw/spinw/total.svg)](https://github.com/spinw/spinw/releases)


<img src="spinw3_logo.png" width="150">

**SpinW** (*spin-double-u*) is a Matlab library that can optimize magnetic structures using mean field theory and calculate spin wave dispersion and spin-spin correlation function for complex crystal and magnetic structures. For details check http://www.spinw.org

# Current Status
We are currently in a period of change. **SpinW will be moving to python/C++ (with a Matlab interface)**. I'm sure you can appreciate that this will be a lot of work as all of the code will be completely re-written and updated. In this period the Matlab version will be stabilized at v3.1.1 with bug fixes and reviewed external pull requests. More details of the new version will follow. [**For Q&A we are testing GitHub Discussions.**](https://github.com/SpinW/spinw/discussions)

# Documentation
* experimental and under construction, the address can change in the future
* documentation of the master branch
* use `swdoc`/`swhelp` instead of the Matlab built-in `doc`/`help` functions to get help on SpinW
* can be also accessed from the browser: https://spinw.github.io/spinwdoc/

# Build Status
We are currently testing on Centos 7.3, macOS and Windows using MATLAB R2018b. It should be noted that MATLAB symbolic calculation changed post R2018a and as such symbolic results may be differ with a relative tolerance of < 0.03 %.

Testing can be pulled from the [testing](https://www.github.com/spinw/testing) repository and run with the `runspinwFunctionalityTests` command from the `Testing` directory.
<!---
### MacOS
## SpinW
### Linux - CentOS 7.3
[![Build Status](https://jenkins.esss.dk/spinw/job/SpinW-CentOS/badge/icon)](https://jenkins.esss.dk/spinw/job/SpinW-CentOS)
### MacOS
[![Build Status](https://jenkins.esss.dk/spinw/job/SpinW-OSX/badge/icon)](https://jenkins.esss.dk/spinw/job/SpinW-OSX)
### Windows 10
[![Build Status](https://jenkins.esss.dk/spinw/job/SpinW-Windows/badge/icon)](https://jenkins.esss.dk/spinw/job/SpinW-Windows/)

## pySpinW
### Linux - CentOS 7.3
[![Build Status](https://jenkins.esss.dk/spinw/job/pySpinW-CentOS-Compile/badge/icon)](https://jenkins.esss.dk/spinw/job/pySpinW-CentOS-Compile)
### MacOS
[![Build Status](https://jenkins.esss.dk/spinw/job/pySpinW-OSX-Compile/badge/icon)](https://jenkins.esss.dk/spinw/job/pySpinW-OSX-Compile)
### Windows 10
[![Build Status](https://jenkins.esss.dk/spinw/job/pySpinW-Windows-Compile/badge/icon)](https://jenkins.esss.dk/spinw/job/pySpinW-Windows-Compile/)
-->
