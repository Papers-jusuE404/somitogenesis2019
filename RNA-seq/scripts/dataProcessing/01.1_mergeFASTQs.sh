### merge FASTq files for each samples
## there are five files for each sample coming from the 5 lanes used for sequencing.
## data is paired-end

# /ark03/repository/do26259-do26350

r1=p1.fq.gz
r2=p2.fq.gz

for i in $(seq 26259 26350)
do
  dir=/ark03/repository/do$i
  
  files=$(ls -d -1 $dir/* | grep $r1)
  cat $files > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2019/RNA-seq/data/FASTq/do$i"_1.fq.gz"
  
  files=$(ls -d -1 $dir/* | grep $r2)
  cat $files > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2019/RNA-seq/data/FASTq/do$i"_2.fq.gz"
done