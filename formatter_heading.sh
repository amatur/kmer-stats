awk '{for (i=1; i<=NF && i<=23; i+=2) printf "%s ", $i; print ""}' $1 
