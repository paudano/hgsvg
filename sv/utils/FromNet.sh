#!/usr/bin/env bash

bioawk -c hdr 'BEGIN{OFS="\t";} { if ($svType == "deletion" && NR > 1) { $svLen=-1*$svLen;} print;}'
