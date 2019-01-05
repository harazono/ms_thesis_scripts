#! /bin/bash
for i in 1 2 3 4 5 ; do
  ./scorerize.rb ../Ecoli/2452_${i}mer_mod.csv > ../minimi/src/score_${i}mer_mod_3digit.h
done

