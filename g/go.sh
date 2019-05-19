#!/bin/sh

rm -f problem
cc problem.c -o problem -lm

rm -f output_*
./problem < input_easy > output_easy
./problem < input_real > output_real
