r=$1
l=100
reps=100
l2=$(echo $l | awk '{printf "%.2f", $1 / 2}')
t=1
for i in `seq 1 $l2`; do
  k=$(echo $i | awk '{printf "%.2f", $1 * 2}')
  ./a.out -l $l -k $k -r $r -t $t -c $reps
done > rt_$r.txt
