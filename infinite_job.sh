#!/bin/bash

while true ; do
   if squeue -u ssilvest | grep 'satori' &> logfile; then
    	    echo "job running"
   else
    	    echo "job stopped, restarting"
    	    cp  run_prototype.sh run_file.sh 
	    new_restart=$(ls -Art | grep RealisticOcean_checkpoint | tail -n 1 | tr -dc '0-9')
	    new_restart=${new_restart:1}
	    new_restart=${new_restart:0:-1}
	    sed -i -e "s/RESTARTHERE/$new_restart/g" run_file.sh 
    	    echo "restaring simulation"
	    ./run_file.sh
    fi	
    sleep 100
done
