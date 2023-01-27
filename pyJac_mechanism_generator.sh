#!/bin/bash

pyjac --lang c --input resources/mechanisms/gri30.cti -b src/Mechanism/CPU
pyjac --lang cuda --input resources/mechanisms/gri30.cti -b src/Mechanism/GPU