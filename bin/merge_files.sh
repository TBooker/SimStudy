
for i in $1/*
   do
   zcat $i >> $2
   echo None >>$2
   echo CUT_HERE >>$2
   done
gzip $2
