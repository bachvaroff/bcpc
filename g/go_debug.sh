#!/bin/sh

rm -f problem
cc -D_DEBUG_ -D_DEBUG_CHECK_ problem.c -o problem -lm

rm -f output_*
./problem < input_easy > output_easy
./problem < input_real > output_real
