#! /bin/sh

filelist="pdblist.inp"

rm $filelist

for file in ../*.cif

do
name=`basename ${file} .cif`
fileout="${name}.dat"
#echo "i $file ${file%.cif}.jpg $fileout"


# do not want TER...
#grep -E "ATOM|HETATM|TER" $file    | grep -vE "REMARK|MASTER|TITL|KEY|JRNL"> $fileout
grep -E "ATOM|HETATM" $file    | grep -vE "REMARK|MASTER|TITL|KEY|JRNL|REVDAT|CAVEAT|HETNAM|MDLTYP|COMPND|MODRES|AUTHOR|SOURCE"> $fileout

echo "done for $file in $fileout"

fileprintlist="${file##*/}"

echo "${fileprintlist%.cif}" >>$filelist


done



echo "DONE, find list in $filelist"  
