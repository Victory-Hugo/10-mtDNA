#########################Examples on single individual ###################################################################

python hGrpr2.py -v data/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf -t data/mt_phyloTree_b17_Tree2.txt -l data/mt_phyloTree_b17_Mutation.txt -f data/rCRS.fasta -i docs/HG00096.txt -o docs/HG00096_mt_hg.txt -x docs/HG00096_mt_allScores.txt


python hGrpr2.py -v data/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf -t data/mt_phyloTree_b17_Tree2.txt -l data/mt_phyloTree_b17_Mutation.txt -f data/rCRS.fasta -w data/mt_phyloTree_b17_Mutation_Wt.txt -i docs/HG00096.txt -o docs/HG00096_mt_hg_wt.txt -x docs/HG00096_mt_allScores_wt.txt


python hGrpr2.py -v data/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf -i docs/HG00096.txt -t data/chrY_hGrpTree_isogg2016.txt -l data/chrY_locusFile_b37_isogg2016.txt -o docs/HG00096_y_hg.txt -x docs/HG00096_y_hg_allScores.txt

############################Phase 3 mitochondrial complete sequence ################################################################

python hGrpr2.py -v data/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf -t data/mt_phyloTree_b17_Tree2.txt -l data/mt_phyloTree_b17_Mutation.txt -f data/rCRS.fasta -o docs/ALL.chrMT.phase3_callmom-v0_4_hg.txt -x docs/ALL.chrMT.phase3_callmom-v0_4_allScores.txt


###########################Phase 3 mitochondrial control region ##########################

python hGrpr2.py -v data/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf -t data/mt_phyloTree_b17_Tree2.txt -l data/mt_phyloTree_b17_Mutation.txt -f data/rCRS.fasta -r 0-600,16024-16365 -o docs/ALL.chrMT.phase3_CR_callmom-v0_4_hg.txt -x docs/ALL.chrMT.phase3_CR_callmom-v0_4_allScores.txt


###########################Phase 3 mitochondrial control region with weights ##########################

python hGrpr2.py -v data/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf -t data/mt_phyloTree_b17_Tree2.txt -l data/mt_phyloTree_b17_Mutation.txt -f data/rCRS.fasta -w data/mt_phyloTree_b17_Mutation_Wt.txt -r 0-600,16024-16365 -o docs/ALL.chrMT.phase3_CR_Wt_callmom-v0_4_hg.txt -x docs/ALL.chrMT.phase3_CR_Wt_callmom-v0_4_allScores.txt


###############################Phase 1 mitochondrial complete sequence#######################


python hGrpr2.py -v data/ALL.chrMT.phase1_samtools_si.20101123.snps.low_coverage.genotypes.vcf -t data/mt_phyloTree_b17_Tree2.txt -l data/mt_phyloTree_b17_Mutation.txt -f data/rCRS.fasta -o docs/ALL.chrMT.phase1_hg.txt -x docs/ALL.chrMT.phase1_hg_allScores.txt


###############################Phase 3 y chromosome #######################
python hGrpr2.py -v data/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf -t data/chrY_hGrpTree_isogg2016.txt -l data/chrY_locusFile_b37_isogg2016.txt -o docs/ALL.chrY.phase3_y_hg.txt -x docs/ALL.chrY.phase3_y_hg_allScores.txt

###############################Phase 3 y chromosome - GSA SNPs #######################

python hGrpr2.py -v data/ALL.chrY.phase3_integrated_v2a_GSA_Y.recode.vcf -t data/chrY_hGrpTree_isogg2016.txt -l data/chrY_locusFile_b37_isogg2016.txt -o docs/ALL.chrY.phase3_integrated_v2a_GSA_y_hg.txt -x docs/ALL.chrY.phase3_integrated_v2a_GSA_y_hg_allScores.txt

###############################Phase 3 y chromosome - using ISOGG 2019  #######################
python hGrpr2.py -v data/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf -t data/treeFileNEW_isogg2019.txt -l data/snpFile_b37_isogg2019.txt -o docs/ALL.chrY.phase3_y_hg_ISOGG2019.txt -x docs/ALL.chrY.phase3_y_hg_ISOGG2019_allScores.txt

