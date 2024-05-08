#!/bin/bash - 
#===============================================================================
#
#          FILE: generate.sh
# 
#         USAGE: ./generate.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (), 
#  ORGANIZATION: 
#       CREATED: 05/06/2024 18:32
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

python3 ../../lib/generate_detbed.py interval_values.csv

ls -d -1 $PWD/*.bed > files.txt
