
### Theorem 12 Verification

To compile Hi Precision / arbitrary precision version   
```
cd thm12
make
```

(The one without any high-precision library is in folder `thm12_low_prec`)

For 100 k-mers, kmer size 40, mutation rate 0.01 and 100 replicates
```
./thm12 -l 100 -k 40 -r 0.01 -c 100
```


Link to results sheet
```
https://docs.google.com/spreadsheets/d/1PQViX9GX8P0pSakuZIdcRMEufI_URUm9JNwdqk6-lso/edit?usp=sharing
```


How to get results for k = 2, 4, 6, ..., 100, type = 0 (random)

Contents of `thm12.sh`
```sh
r=$1
l=100
reps=100
t=0
for i in `seq 1 50`; do
  k=$(echo $i | awk '{printf "%.2f", $1 * 2}')
  ./thm12 -l $l -k $k -r $r -t $t -c $reps
done > rt_$r.txt
```

```sh
sh thm12.sh 0.01
```




### Theorem 3 Verification

You need to compile and run
`thm3/thm3.cpp`

Link to results 
```
https://docs.google.com/spreadsheets/d/1N6YMjxQiQEvTJs_nvwuPI4gxseqN1T2BwlPfQHgqdV4/edit?usp=sharing
```




### Fast Hamming Distance, Sketching Code

You need to compile and run `hd_histogram_sketch/hamming.cpp`

Results are in 
```
https://docs.google.com/spreadsheets/d/1CFozuonz5w_aTjx0C7sIWCAan_u5IHa8vjecYXfg-0E/edit?usp=sharing
```

