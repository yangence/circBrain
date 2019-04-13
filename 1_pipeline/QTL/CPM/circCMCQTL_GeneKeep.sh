#!/usr/bin/sh
for i in {1..22}
do
{
/opt/microsoft/ropen/3.4.2/lib64/R/bin/Rscript /media/data3/circCMC/script/circCMCQTL_GeneKeep.R  $i >/media/data3/circCMC/data/qtl_C/GeneKeep/log/$i.log 2>&1
}&
done
