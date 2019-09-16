#!/bin/bash
TIMESTAMP=`date +"%m-%d-%Y"`
if [ -e $HOME/$TIMESTAMP ]
 then
  echo "directory $HOME/$TIMESTAMP exists"
 else
  mkdir $HOME/$TIMESTAMP
fi
