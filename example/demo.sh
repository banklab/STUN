# author: Li Juan

popSize=5000
add=0.01
numLoci=15
driftThreshold=3

for epi in 0 0.01 0.1
do
	../target/release/stun -N $popSize --neutralsfs $driftThreshold 1 --replicates 1 --landscapes 10 RMF -L $numLoci -m $add -s 0 -S $epi
done

Rscript ./plot.R ./data/

exit 0
