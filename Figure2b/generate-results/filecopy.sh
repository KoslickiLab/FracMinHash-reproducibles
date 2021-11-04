F=/data/dmk333/repos/Metalign/data/organism_files

for line in $(cat file_list)
do
   cp $F/$line .
done
