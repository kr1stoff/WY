#!/bin/sh 

head=`head -1 FBS-vs-RZSXF.P.xls | awk '{print "Class\t"$0}'`
sed 1d FBS-vs-RZSXF.P.xls | awk '{print "BP\t"$0}' >  FBS-vs-RZSXF.GO.xls
sed 1d FBS-vs-RZSXF.F.xls | awk '{print "MF\t"$0}' >> FBS-vs-RZSXF.GO.xls
sed 1d FBS-vs-RZSXF.C.xls | awk '{print "CC\t"$0}' >> FBS-vs-RZSXF.GO.xls

sed -i 1i"$head" FBS-vs-RZSXF.GO.xls

perl ../../enrichPlot.pl -o FBS-vs-RZSXF.GO.enrich FBS-vs-RZSXF.GO.xls
perl ../../enrichPlot.pl -d ../FBS-vs-RZSXF.glist -o FBS-vs-RZSXF.GO.enrich.updown FBS-vs-RZSXF.GO.xls
