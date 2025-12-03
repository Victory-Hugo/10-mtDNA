HaploGrouper
========

A software to classify haplotypes into haplogroups on the basis of a known phylogenetic tree.

## Requirements

HaploGrouper has been run successfully using both Python 2.7 and 3.6. It will likely work using similar versions.

It requires the `numpy` library.


## Running HaploGrouper 

```sh
# clone gitlab repository
git clone https://gitlab.com/bio_anth_decode/haploGrouper.git

# run HaploGrouper
cd haploGrouper/

python hGrpr2.py -h
usage: hGrpr2.py [-h] -v VCFFILE -o OUTFILE -l HGRPLOCUSFILE -t HGRPTREEFILE
                 [-i IDLISTFILE] [-r REGIONS] [-c CHROM] [-f REFERENCEFASTA]
                 [-m] [-w WEIGHTFILE] [-x VERBOSEFILE]

Determine haplogroup for list of individuals based on VCF file and haplogroup
tree. ----

optional arguments:
  -h, --help            show this help message and exit
  -v VCFFILE, --vcfFile VCFFILE
                        path of vcfFile to be converted (can be gzipped vcf
                        file)
  -o OUTFILE, --outFile OUTFILE
                        path of output file
  -l HGRPLOCUSFILE, --hGrpLocusFile HGRPLOCUSFILE
                        path of file information about loci used to assign
                        individuals to haplogroups and the branches in the
                        haplogroup tree
  -t HGRPTREEFILE, --hGrpTreeFile HGRPTREEFILE
                        path of file information all branches in the
                        haplogroup tree
  -i IDLISTFILE, --IDListFile IDLISTFILE
                        path of file with subset of IDs from VCF file that are
                        to be used for haplogroup assignment
  -r REGIONS, --regions REGIONS
                        Only base haplogroup assignment on positions from the
                        specified regions (format: startPos-stopPos,startPos-
                        stopPos) Multiple regions can be specified
  -c CHROM, --chrom CHROM
                        Only process loci from vcf file with this chromosome
                        name. This must be used when the vcf file contains
                        loci from multiple chromosomes.
  -f REFERENCEFASTA, --referenceFasta REFERENCEFASTA
                        Read reference sequence from fasta file. This is
                        useful when the vcfFile is based on full sequence
                        data, but only reports polymorphic positions or
                        differences from the reference sequence. In such
                        cases, many phylogenetically informative positions
                        would be ignored. When the reference sequence is
                        provided, positions not reported in the vcfFile are
                        assumed to have the reference state for all
                        individuals in the file.
  -m, --mismatchLoc     Report genotypes of mismatching loci in outFile
  -w WEIGHTFILE, --weightFile WEIGHTFILE
                        Give mutations differing weights read from a separate
                        file. This is only useful for loci with a high rate of
                        recurrent mutations - like the mtDNA control region.
                        The file should be tab-delimited with four columns:
                        pos ancAl derAl weight. Mutations not included in the
                        file default to a weight of 1
  -x VERBOSEFILE, --verboseFile VERBOSEFILE
                        path of file for full matrix of scores for each node
                        in the tree (lines) for each individual (columns)
```


## Example 1: running on a single individual to assign mitochondrial haplogroup

```sh
echo "HG00096" > docs/HG00096.txt 

python hGrpr2.py -v data/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf \
                 -t data/mt_phyloTree_b17_Tree2.txt \
                 -l data/mt_phyloTree_b17_Mutation.txt \
                 -f data/rCRS.fasta \
                 -i docs/HG00096.txt \
                 -o docs/HG00096_mt_hg.txt \
                 -x docs/HG00096_mt_allScores.txt  
``` 
 
The -o parameter is mandatory, and specifies where the haplogroup label file is written to (see `docs/HG00096_mt_hg.txt`). 
 
 
### Example output 1: mitochondrial haplogroup assigned for HG00096

<table>
<thead>
<tr>
<th>ID</th>
<th>HaploGroup</th>
<th>netScore</th>
<th>matchScore</th>
<th>mismatchScore</th>
<th>mismatchLoci</th>
<th>backMutCnt</th>
<th>pruning</th>
<th>allMaxNetScore</th>
</tr>
</thead>
<tbody>
<tr>
<td>HG00096</td>
<td>H16a1</td>
<td>45</td>
<td>48</td>
<td>3</td>
<td> </td>
<td>3</td>
<td> </td>
<td>H16a1[48-3]</td>
</tr>
</tbody>
</table>


## Example 2: running on a single individual to assign mitochondrial haplogroup using the weight option

```sh
python hGrpr2.py -v data/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf \
                 -t data/mt_phyloTree_b17_Tree2.txt \
                 -l data/mt_phyloTree_b17_Mutation.txt \
                 -f data/rCRS.fasta \
                 -w data/mt_phyloTree_b17_Mutation_Wt.txt \
                 -i docs/HG00096.txt \
                 -o docs/HG00096_mt_hg_wt.txt \
                 -x docs/HG00096_mt_allScores_wt.txt
```

### Example output 2: mitochondrial haplogroup assigned for HG00096 using the weight option

<table>
<thead>
<tr>
<th>ID</th>
<th>HaploGroup</th>
<th>netScore</th>
<th>matchScore</th>
<th>mismatchScore</th>
<th>mismatchLoci</th>
<th>backMutCnt</th>
<th>pruning</th>
<th>allMaxNetScore</th>
</tr>
</thead>
<tbody>
<tr>
<td>HG00096</td>
<td>H16a1</td>
<td>460.1</td>
<td>23.1</td>
<td> </td>
<td>3</td>
<td> </td>
<td>H16a1[460.1-23.1]</td>
</tr>
</tbody>
</table>

## Example 3: running on a single individual to determine Y-chromosome haplogroup

```sh
python hGrpr2.py -v data/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf \
                 -i docs/HG00096.txt \
                 -t data/chrY_hGrpTree_isogg2016.txt \
                 -l data/chrY_locusFile_b37_isogg2016.txt \
                 -o docs/HG00096_y_hg.txt \
                 -x docs/HG00096_y_hg_allScores.txt
 ```

### Example output 3: Y-chromosome haplogroup assigned for HG00096 

<table>
<thead>
<tr>
<th>ID</th>
<th>HaploGroup</th>
<th>netScore</th>
<th>matchScore</th>
<th>mismatchScore</th>
<th>mismatchLoci</th>
<th>backMutCnt</th>
<th>pruning</th>
<th>allMaxNetScore</th>
</tr>
</thead>
<tbody>
<tr>
<td>HG00096</td>
<td>R1b1a2a1a2c1k1</td>
<td>698</td>
<td>698</td>
<td>0</td>
<td> </td>
<td>0</td>
<td> </td>
<td>R1b1a2a1a2c1k1[698-0]</td>
</tr>
</tbody>
</table>

## Example 4: running on a single individual to assign Y-chromosome haplogroup using ISOGG 2019 files

```sh
python hGrpr2.py -v data/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf \
                 -t data/treeFileNEW_isogg2019.txt \
                 -l data/snpFile_b37_isogg2019.txt \
                 -i docs/HG00096.txt \
                 -o docs/HG00096_y_hg_ISOGG2019.txt \
                 -x docs/HG00096_y_hg_allScores_ISOGG2019.txt
 ```

### Example output 4: Y-chromosome haplogroup assigned for HG00096 using ISOGG 2019

<table>
<thead>
<tr>
<th>ID</th>
<th>HaploGroup</th>
<th>netScore</th>
<th>matchScore</th>
<th>mismatchScore</th>
<th>mismatchLoci</th>
<th>backMutCnt</th>
<th>pruning</th>
<th>allMaxNetScore</th>
</tr>
</thead>
<tbody>
<tr>
<td>HG00096</td>
<td>R1b1a1b1a1a2c1a1f1</td>
<td>952</td>
<td>952</td>
<td>0</td>
<td> </td>
<td>0</td>
<td> </td>
<td>R1b1a1b1a1a2c1a1f1[952-0]</td>
</tr>
</tbody>
</table>

## Support

If you experience problems using HaploGrouper, please [report an issue on the GitLab repository here](https://gitlab.com/bio_anth_decode/haploGrouper/-/issues).

If you cannot use the GitLab issues tracker for whatever reason, you can also email `agnar [at-symbol] decode [dot] is` and `kristjanm [at-symbol] decode [dot] is`.