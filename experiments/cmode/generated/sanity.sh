#!/bin/bash
python3 ../generate_beds.py remap2022_crm_macs2_hg38_v1_0.bed

realpath remap2022_crm_macs2_hg38_v1_0_*.bed > remap_beds.txt 
realpath remap2022_crm_macs2_hg38_v1_0.bed > remap_prime.txt

python3 ../../lib/bed_sourmash_1m.py remap_beds.txt remap_prime.txt -p 50 -o ./remap2022_generated --mode B

realpath $(ls remap*mode?_[2,3,4].bed) > remap_beds_short.txt 

python3 ../../lib/bed_sourmash_1m.py remap_beds_short.txt remap_prime.txt -p 30 -o ./remap2022_generated_short --mode C 
