#!/bin/bash

# Mafft each .fasta 
for i in 1_raw/*.fasta; do mafft $i > ${i%}.aligned.fas; mv 1_raw/*.fas 2_aligned/; done
done
