#! /bin/bash

cd $1
#echo "'$3'"
git add $2
#git commit -m "'$3'"
git commit -m $3
git push -u origin master

