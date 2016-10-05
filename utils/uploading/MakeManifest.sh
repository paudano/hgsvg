#!/usr/bin/env bash
if [ "$#" -ne 2 ]; then
   echo "Usage: MakeManifest.sh file_or_directory manifest.tsv"
   exit 1
fi
~mchaisso/software/gca-tools/submissions/file-manifest $1 > $2
