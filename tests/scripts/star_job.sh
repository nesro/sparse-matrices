#!/bin/bash
#
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
# https://github.com/nesro/sparse-matrices
#
# This is a helper file for running jobs on our school server STAR.
# Please do NOT to try hack this in any way.

make KAT_N=4 tests
time ./tests/bin/test_mat_vec