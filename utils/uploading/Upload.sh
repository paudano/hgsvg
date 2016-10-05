#!/usr/bin/env bash
module load aspera/latest
if [ "$#" -lt 1 ]; then
  echo "Usage: Upload.sh file1 [file2 file3 ...]"
  exit 1
fi

export ASPERA_SCP_PASS=3fnMsgrU
while test $# -gt 0; 
 do 
		 ascp -q -k 1 -QT -l 300M $1 drop-sv@fasp.sra.ebi.ac.uk:./
		 ts=`date`
		 echo $ts	$PWD $1 >> /net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/Uploads/upload_log.txt
shift
done

