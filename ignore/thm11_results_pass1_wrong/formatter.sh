awk '{for (i=2; i<=NF && i<=18; i+=2) printf "%s ", $i; print ""}' $1 
