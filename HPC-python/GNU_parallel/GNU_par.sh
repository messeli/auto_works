# start=`date +%s.%N`
# python3 write-Omeg_range.py "5.01,7.02,0.1"
# end=`date +%s.%N`  #$(date +%s.%N)
# runtime=$($end-$start)
# echo Time taken `expr $end-$start` #$runtime

date +%s.%N
python3 write-Omeg_range.py "5.01,7.02,0.1"
date +%s.%N
parallel python3 script_GNU_parallel.py :::: Omeg_range.txt
date +%s.%N
python3 plot.py
date +%s.%N

