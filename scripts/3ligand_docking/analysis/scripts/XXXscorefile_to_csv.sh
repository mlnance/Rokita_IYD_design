#!/bin/bash

#Example
# ./scorefile_to_csv.sh scorefile.sc > scorefile.csv

#Given an input scorefile.sc file $1,
#replace all spaces with commas
#condense multiple commas to single comma
#May also need after run is an emacs command
#to remove commas that trail newline characters
#Esc-x replace-string , Ctl-q Ctl-j, Ctl-q Ctl-j

#then use Pandas to import into a DataFrame
# df = pandas.DataFrame().from_csv("scorefile.csv", index_col=None, header=1)


out="$(echo $1 | cut -d'.' -f1)"
cat $1 | tr ' ', ',' | tr -s ','

#link to website for commands
#https://stackoverflow.com/questions/25908258/substitute-space-with-comma-in-shell-scripting
