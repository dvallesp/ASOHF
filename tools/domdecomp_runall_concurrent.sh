#!/bin/bash

echo "Folders to visit:"
for DOM_FOLDER in domain_*;
do echo "$DOM_FOLDER"
done

for DOM_FOLDER in domain_*;
do
	echo ""
	echo "************************"
	echo "$DOM_FOLDER"
	echo "************************"
	echo ""
       	cd $DOM_FOLDER 
	./run.sh
	cd ..	
done
