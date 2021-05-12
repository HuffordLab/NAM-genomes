#!/bin/bash
input=input.vcf

snphylo.sh -v $input -p 70 -l 0.4 -M 0.7 -A -b &> default-p70-missing-0.7-ld-0.4
