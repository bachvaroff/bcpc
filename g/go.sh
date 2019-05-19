#!/bin/sh

cc problem.c -o problem -lm
./problem < input_easy > output_easy
./problem < input_real > output_real
