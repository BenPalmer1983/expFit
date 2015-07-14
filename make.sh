#!/bin/bash
cd ${0%/*}
mpif90 -g -Wall -Wno-unused-function -O3 \
-fbounds-check -fcheck=all -mtune=native \
expFit.f90 -o expFit.x