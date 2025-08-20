#!/bin/bash
./assignBarcodes -v  -w  /mnt/pikachu/processing/genome/10x_version3_whitelist.txt -o 20 -m 2 -t 4 -s 1 -i 0 -u 12 --limit_search 2 -c 30   --max_barcode_mismatches 1 -f /mnt/pikachu/ref_feature_geneBC.csv --translate_NXT -d test_output/30_KO_ES /storage/gene/30_KO_ES
