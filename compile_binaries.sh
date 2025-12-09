#!/bin/bash
set -e  # stop if any command fails

# Build variants
variants=(
  ""                 # default (no macro → ORd-static)
  "-DCIPAORDV1_0"    # CiPAORdv1.0
  "-DTOR_ORD"        # ToR-ORd
  "-DTOR_ORD_DYNCL"  # ToR-ORd-dynCl
  "-DGRANDI"  # Grandi
)

# for the full binary
for flag in "${variants[@]}"; do
  make CXXFLAGS="$flag" clean
  make all CXXFLAGS="-std=c++11 $flag"
done

# for the postprocessing
for flag in "${variants[@]}"; do
  make CXXFLAGS="-DPOSTPROCESSING $flag" clean
  make all CXXFLAGS="-std=c++11 -DPOSTPROCESSING $flag"
done
