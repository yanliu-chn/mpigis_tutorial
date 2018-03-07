#!/bin/bash
t1=`date +%s`
bindir=`dirname $0`
for job in 1 2 3 4 5 6 7 8
do
	$bindir/dosth.sh $job 3
done
t2=`date +%s`
let "t=$t2-$t1"
echo "Total time: $t seconds."
