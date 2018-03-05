#!/usr/bin/env bash

bioawk -c hdr '{ if ($svType == "insertion" && NR > 1) { $tEnd=$tStart+$svLen;} print;}'
