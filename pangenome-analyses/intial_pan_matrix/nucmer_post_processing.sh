# get every file listed in the DIR folling prefix
for FILE in $(ls *_c1000.delta); do
    NAME=${FILE%%.*}
    DELTANAME="$NAME.fil.delta"
    TAB=$'\t'

    # Delta filter for low matching syntenic
    delta-filter \
        -g \
        -u 75 \
        $FILE \
        > $DELTANAME

    # convert to a table set of coordinates
    show-coords \
        -c -l \
        -r -T \
        $DELTANAME \
         > $NAME.fil.coords

    # remove first three lines with sed
    sed -i '1,3d' $NAME.fil.coords
    # change the first line so it can be read by pandas
    sed -i "1 s/^.*$/S1${TAB}E1${TAB}S2${TAB}E2${TAB}LEN1${TAB}LEN2${TAB}%IDY${TAB}LENR${TAB}LENQ${TAB}COVR${TAB}COVQ${TAB}TAGS${TAB}TAGS2/" $NAME.fil.coords
done
