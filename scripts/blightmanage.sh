kff_one_sec=$(realpath $1)
DIR=$(dirname $kff_one_sec)
K=$2
M=$3

cd $DIR
bz=$(basename $kff_one_sec .kff)
echo $bz
kff-tools outstr -i $kff_one_sec >  $bz.outstr
cat $bz.outstr | cut -f1 -d" " | awk '{print ">contig\n" $0}' > $bz.fa
essCompress -i $bz.fa -k $K -u 
./snippet $bz.outstr $bz.fa.essd $K $bz.instr

