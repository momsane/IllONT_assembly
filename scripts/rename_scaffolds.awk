#!/bin/awk -f

BEGIN {
    RS = ">"
    FS = "\n"
}

NR == 1 {
    # skip empty record created by initial ">"
    next
}

NF > 1 {
    #split header 
    n=split($1, header, "_")
    # rename contig contigid_length_xxx_circular or contigid_length_xxx_linear
    if (header[n] == "true") {
        print ">"header[1]"_len_"header[3]"_circular"
    } else {
        print ">"header[1]"_len_"header[3]"_linear"
    }
    for (i=2; i<=NF; i++) {
      print $i
    }
}

