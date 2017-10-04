#! /bin/bash

cd $1
#echo $1
#echo $2
#echo $3
git add $2
#msg="'$3'"
#echo $msg
#echo $3
#git commit -m "'$3'"
git commit -m $3
git push -u origin master

