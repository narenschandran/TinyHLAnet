# Inputs:
# [1] - Mandatory input: should either be a file with MixMHCpred results or
#       should be the string "header" to get only the output header
# [2] - Mandatory input: Allele name
# [3] - Optional input: if input [1] is a file, any input [3] will trigger
#       the parsed result to also have a header.

f="$1"
allele="$2"

if [ "$f" = "header" ]; then
    # We can request the header row separately if necessary
    echo -e "allele\tpeptide\tbindscore"
else
    # The MixMHCpred results start in tabular format after a header.
    # We track when the header start (with the column "Peptide")
    # and parse the results until the end of the file.
    if [ -z "$3" ]; then
        hd=0
    else
        hd=1
    fi
    awk -vallele="$allele" -vhd="$hd" '
    BEGIN {
        # We print the header if the option is toggled.
        if (hd == 1) {
            print "allele\tpeptide\tbindscore"
        }


        start_next = 0
        parse = 0
    }

    (start_next == 1) {
        start_next = 0
        parse = 1
    }

    $1 ~ /^Peptide/ {
        start_next = 1
    }

    (parse == 1) {
        print allele"\t"$1"\t"$5
    }

    ' $f
fi
