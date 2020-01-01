#!/bin/bash
nam=$1
cd /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}
icd /iplant/home/shared/NAM/Transcript_Assemblies
RESULT=$?
if [ $RESULT -eq 0 ]; then
  iput -a -K -P -r -T ${nam}/
else
  echo "upload failed: cyverse directory change unsuccessful"
fi
