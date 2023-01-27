#!/bin/bash
cd ext/cantera/
scons build prefix=../ext_install/cantera python_package=none
