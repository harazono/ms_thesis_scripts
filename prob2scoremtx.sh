#! /bin/bash
for i in 1 2 3 4 5 6 ; do
  ./scorerize.rb ../Ecoli/2452_${i}mer_mod.csv > ./kenpin_${i}.h
done

