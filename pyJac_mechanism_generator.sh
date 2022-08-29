#!/bin/bash

pyjac --lang c --input ext_install/cantera/share/cantera/data/gri30.cti -b src/Mechanism/GPU
pyjac --lang cuda --input ext_install/cantera/share/cantera/data/gri30.cti -b src/Mechanism/GPU