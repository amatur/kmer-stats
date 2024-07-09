
r=0.01
for i in `seq 1 100`; do
  n=$(echo $i | awk '{printf "%.2f", $1 / 2}')
  ./a.out -a $n -r $r
done > rt_$r.txt

