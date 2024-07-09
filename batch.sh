for i in `seq 1 100`; do
  n=$(echo $i | awk '{printf "%.2f", $1 / 2}')
  ./a.out -a $n -r 0.5
done > r0.5.txt

