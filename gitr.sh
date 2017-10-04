#! /bin/bash

cd $1
#echo $1
#echo $2
#echo $3
#echo "'$3'"
git add $2
#git commit -m "'$3'"
git commit -m $3
git push -u origin master

