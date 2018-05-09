#!/usr/bin/env bash

BASE="$( cd "$(dirname "$0")" ; pwd -P )"
usage()
{
cat << EOF
		RunTrioTiledAssemblyOnRegions.sh regions paramfile 
		-j job  Will submit jobs under this name to sge.
		-a asm  The full path to the assembler makefile to use.
    -d dir  Run in this directory
EOF
exit 1

}


while true; do 
     case "$1" in
         -h | --help)
             usage
             exit 1
             ;;
         --vcf )
             VCF=$2
						 shift
						 shift
             ;;
				 --fa )
						 FA=$2
						 shift
						 shift
						 ;;
				 --mo )
						 MO=$2
						 shift
						 shift
						 ;;
				 --child )
						 CH=$2
						 shift
						 shift
						 ;;
				 --faBed )
						 faBed=$2
						 shift
						 shift
						 ;;
				 --moBed )
						 moBed=$2
						 shift
						 shift
						 ;;				 
         -? )
						 usage
						 shift
             exit
             ;;
				 * )
						 break
     esac
done

echo $faBed
$BASE/DetermineInheritance.py --vcf $VCF --fa $FA --mo $MO --child $CH --faBed $faBed --moBed $moBed

bgzip $faBed
tabix $faBed.gz
bgzip $moBed
tabix $moBed.gz




