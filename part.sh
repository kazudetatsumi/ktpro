#!/bin/bash
j=0
string=""
for i in `grep grid_point ir_grid_points.yaml | grep -|grep -o '[0-9]*'` ; do
    if [ $((j%10)) = 0 -a $i != 0 ];then
       echo $string
       mkdir ./gpjob_$j
	   cd ./gpjob_$j
       cp  ../../for_gpjob/* .
       sed s/GP/"$string"/ ../job-ph3-part.sh > job.sh    
       qsub ./job.sh
       cd ../ 
       string=""
    #echo $j
    fi
    string=$string" "$i
    j=$((j+1))
done 

if [ $((j%10)) != 1 ];then
       echo $string
       mkdir ./gpjob_$j
	   cd ./gpjob_$j
       cp  ../../for_gpjob/* .
       sed s/GP/"$string"/ ../job-ph3-part.sh > job.sh    
       qsub ./job.sh
       cd ../ 
fi
