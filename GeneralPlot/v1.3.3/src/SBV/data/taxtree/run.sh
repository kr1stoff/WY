# prepare the input file 
perl /home/aipeng/work/develepment/SBV/src/tax2nhx.pl all.tsv

# change the config file to set the 'file' and 'percent' path which created before

# draw figure
perl /home/aipeng/work/develepment/SBV/bin/sbv.pl taxtree -conf ../../etc/taxtree.conf -out taxtree --strict
