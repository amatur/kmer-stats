r=$1
l=100
reps=100
lx=100
#l2=$(echo $lx | awk '{printf "%.2f", $1 / 5}')
t=0
#for i in 500 600 700; do
for i in `seq 1 50`; do
  k=$(echo $i | awk '{printf "%.2f", $1 * 2}')
  ./a.out -l $l -k $k -r $r -t $t -c $reps
done > rt_$r.txt
