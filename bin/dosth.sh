#!/bin/bash
myid=$1
if [ -z $myid ]; then
	myid=0
fi
duration=$2
if [ -z $duration ]; then
	duration=10
fi
t1=`date +%s`
echo "$myid: [`date`] start"
echo "$myid: doing something on `hostname -f` for `whoami` for $duration seconds..."
sleep $duration
t2=`date +%s`
let "t=$t2-$t1"
echo "$myid: [`date`] done in $t seconds."
