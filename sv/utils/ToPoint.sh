#!/usr/bin/env bash

bioawk -c hdr 'BEGIN{OFS="\t";} { if ($svType == "insertion" && NR > 1) { $tEnd=$tStart+1;} print;}'
