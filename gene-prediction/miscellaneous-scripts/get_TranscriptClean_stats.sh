# Script that uses the output files from TranscriptClean to obtain
# summary stats about the run

inputSam=$1
prefix=$2

outputSam=${prefix}_clean.sam
logTE=${prefix}_clean.TE.log

# Count the total number of transcripts in the input SAM
inputTranscripts=$(grep -v "^@" $inputSam | wc -l)
echo "Number of transcripts in input: ${inputTranscripts}"

# Count the total number of reads processed by
# TranscriptClean (ie not multimapped or unmapped
processedTranscripts=$(grep -v "^@" $outputSam | wc -l)
echo "Number of transcripts processed: ${processedTranscripts}"

percentReduction=$(echo "($processedTranscripts - $inputTranscripts)*100.0 / $inputTranscripts" | bc)
echo "Percent change in numbr of transcripts: ${percentReduction}%"

# Count the number of deletions in the data before and after
deletionsBefore=$(grep "Deletion" ${logTE} | wc -l)
deletionsAfter=$(grep "Deletion" ${logTE} | grep "Uncorrected" | wc -l)
echo "Number of deletions in the data before and after TranscriptClean: ${deletionsBefore}, ${deletionsAfter}"

percentReduction=$(echo "($deletionsAfter - $deletionsBefore)*100.0 / $deletionsBefore" | bc)
echo "Percent reduction in deletions: ${percentReduction}%"

# Count the number of insertion in the data before and after
insertionsBefore=$(grep "Insertion" ${logTE} | wc -l)
insertionsAfter=$(grep "Insertion" ${logTE} | grep "Uncorrected" | wc -l)
echo "Number of insertions in the data before and after TranscriptClean: ${insertionsBefore}, ${insertionsAfter}"

percentReduction=$(echo "($insertionsAfter - $insertionsBefore)*100.0 / $insertionsBefore" | bc)
echo "Percent reduction in insertions: ${percentReduction}%"

# Count the number of mismatches in the data before and after
mismatchesBefore=$(grep "Mismatch" ${logTE} | wc -l)
mismatchesAfter=$(grep "Mismatch" ${logTE} | grep "Uncorrected" | wc -l)
echo "Number of mismatches in the data before and after TranscriptClean: ${mismatchesBefore}, ${mismatchesAfter}"

percentReduction=$(echo "($mismatchesAfter - $mismatchesBefore)*100.0 / $mismatchesBefore" | bc)
echo "Percent reduction in msimatches: ${percentReduction}%"

# Count the number of noncanonical splice junctions in the data before and after
ncBefore=$(grep "NC_SJ" ${logTE} | wc -l)
ncAfter=$(grep "NC_SJ" ${logTE} | grep "Uncorrected" | wc -l)
echo "Number of noncanonical splice junctions in the data before and after TranscriptClean: ${ncBefore}, ${ncAfter}"

percentReduction=$(echo "($ncAfter - $ncBefore)*100.0 / $ncBefore" | bc)
echo "Percent reduction in noncanonical splice junctions: ${percentReduction}%"
