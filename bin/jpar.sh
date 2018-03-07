#!/bin/bash
t1=`date +%s`
bindir=`dirname $0`
plist=""
for job in 1 2 3 4 5 6 7 8
do
	$bindir/dosth.sh $job 3 &
	plist="$plist $!"
done
t2=`date +%s`
let "t=$t2-$t1"
echo "Jobs submitted in $t seconds."
for p in $plist
do
	wait $p
done
t3=`date +%s`
let "t=$t3-$t2"
let "tt=$t3-$t1"
echo "Jobs completed in $t seconds."
echo "Total time: $tt seconds."
