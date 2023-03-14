#snv
mkdir -p snv
for i in `ls FBC*.vcf`; do grep "#" $i > header; id=`basename $i .vcf`;cat header <( grep -v "#" $i |awk '{if(length($4)==1 && length($5)==1)print}') > snv/${id}.snv.vcf;done

#indel
mkdir -p indel
for i in `ls FBC*.vcf`; do grep "#" $i > header; id=`basename $i .vcf`; cat header <(grep -v "#" $i |grep "INDEL") > indel/${id}.indel.vcf;done

rm -f header
