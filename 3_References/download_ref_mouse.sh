#!/bin/bash
 
VERSION=113 # curl -s 'http://rest.ensembl.org/info/software?content-type=application/json' | grep -o '"release":[0-9]*' | cut -d: -f2
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/gtf/mus_musculus/Mus_musculus.GRCm39.$VERSION.gtf.gz

