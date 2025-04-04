for file in $(ls *.bam)
do
    base=$(basename $file ".bam")
    /scratch/monoallelic/clt54/Nanopore/dist_modkit_v0.4.1_cec0a0b/modkit pileup \
    --ref /data/monoallelic/christian/Genomes/Nvit_psr_1.1/Nvit_psr_1.1/GCF_009193385.2_Nvit_psr_1.1_genomic.fa --cpg \
    --threads 16 \
    --prefix "$base" \
    # --combine-strands \ #remove this hash before this to make data unstranded
    "$file" "$base"_meth.bed
done
