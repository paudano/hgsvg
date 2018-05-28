#!/usr/bin/env bash

BASE="$( cd "$(dirname "$0")" ; pwd -P )"
usage()
{
cat << EOF
		DetermineInheritance.sh [all options required]
		--vcf Trio phased vcf.
		--fa  ID of father.
    --mo  ID of mother.
    --child ID of child.
    --faBed output bed file for father.
    --moBed outpub bed file for mother.    
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




