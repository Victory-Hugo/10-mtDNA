#!/usr/bin/env python

import sys
import argparse
from time import clock, time
import os
import numpy as np
import hGrpUtil


'''
treeParFile=${workDir}/mt_phyloTree_b17_Tree2.txt

inFile1=/mnt/c/ritvinn/Greinar/haplogrouper/code/haplogrouper/mt_1000genomes_truthSet.txt

inFile2=/mnt/c/ritvinn/Greinar/haplogrouper/code/nov2019/data/ALL.chrMT.phase3_callmom-v0_4.20130502_parent_wRef_mtHGrps.txt
outFile=/mnt/c/ritvinn/Greinar/haplogrouper/code/nov2019/data/mtDNA_dist_truth_v_wRef_full.txt

inFile2=/mnt/c/ritvinn/Greinar/haplogrouper/code/nov2019/data/ALL.chrMT.phase3_callmom-v0_4.20130502_parent_wRef_HVR12_mtHGrps.txt
outFile=/mnt/c/ritvinn/Greinar/haplogrouper/code/nov2019/data/mtDNA_dist_truth_v_wRef_HVR12.txt


python /mnt/c/ritvinn/Greinar/haplogrouper/code/haplogrouper/hGrpDist.py -i1 ${inFile1} -i2 ${inFile2} -o ${outFile} -t ${treeParFile}


'''
    

def calcDist(hGrp1, hGrp2, hGrpTree):
  if hGrp1 in hGrpTree.nodeDict and hGrp2 in hGrpTree.nodeDict:
   hGrp1Idx = hGrpTree.nodeDict[hGrp1]
   hGrp2Idx = hGrpTree.nodeDict[hGrp2]
  
   hGrp1PathLen = hGrpTree.nodeList[hGrp1Idx].pathArr.size
   hGrp2PathLen = hGrpTree.nodeList[hGrp2Idx].pathArr.size
   pathLen = min(hGrp1PathLen, hGrp2PathLen)
  
 
   for pathIdx in range(pathLen):
     if hGrpTree.nodeList[hGrp1Idx].pathArr[pathIdx] != hGrpTree.nodeList[hGrp2Idx].pathArr[pathIdx]:
       break
   mrcaIdx = hGrpTree.nodeList[hGrp1Idx].pathArr[pathIdx]
  
   hGrp1PathLeft = hGrp1PathLen - (pathIdx + 1)
   hGrp2PathLeft = hGrp2PathLen - (pathIdx + 1)
   hGrp12Dist = hGrp1PathLeft + hGrp2PathLeft

   return hGrpTree.nodeList[mrcaIdx].name, pathIdx, hGrp1PathLeft, hGrp2PathLeft, hGrp12Dist
  else:
   return 99999,99999,99999,99999,99999

 

def main(args):  

  startTime = clock()

 
  if os.path.isfile(args.hGrpTreeFile) == False:
    sys.stdout.write("Aborting! hGrpTreeFile not found: %s\n" % (args.hGrpTreeFile))
    return 0
  
  if os.path.isfile(args.inFile1) == False:
    sys.stdout.write("Aborting! inFile1 not found: %s\n" % (args.inFile1))
    return 0

  if os.path.isfile(args.inFile2) == False:
    sys.stdout.write("Aborting! inFile2 not found: %s\n" % (args.inFile2))
    return 0
    
    
    
  sys.stdout.write("## Reading haplogroup tree and locus files\n")
  hGrpTree = hGrpUtil.readHGrpTreeFile(args.hGrpTreeFile)

  f1ColList = [int(colNr) for colNr in args.colSTR1.split(":")]
  f2ColList = [int(colNr) for colNr in args.colSTR2.split(":")]
  
  outStream = open(args.outFile, 'w')
  outStream.write("id\thgrp1\thgrp2\tmrcahgrp\tmrcaRootDist\thGrp1tailDist\thGrp2tailDist\ttailDistSum\n")

  # get IDs and hGrps from file 1
  f1Dict = {}
  inStream1 = open(args.inFile1)
  for line in inStream1:
    line = line.rstrip('\r\n')
    f = line.split('\t')
    id = f[f1ColList[0]]
    hgrp = f[f1ColList[1]]
    f1Dict[id] = hgrp
  inStream1.close()
  
  inStream2 = open(args.inFile2)
  for line in inStream2:
    line = line.rstrip('\r\n')
    #print(line)
    f = line.split('\t')
    id = f[f2ColList[0]]
    hgrp = f[f2ColList[1]]
    if id in f1Dict:
      mrcaName, mrcaRootDist, hGrp1PathLeft, hGrp2PathLeft, hGrp12Dist = calcDist(f1Dict[id], hgrp, hGrpTree)
      outStream.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n" % (id, f1Dict[id], hgrp, mrcaName, mrcaRootDist, hGrp1PathLeft, hGrp2PathLeft, hGrp12Dist))
  
  inStream2.close()
  

  # process haplogroup tree and branch mutations
  
  sys.stderr.write("Finished in %.2f seconds\n" % (clock()-startTime))
    


if __name__ == '__main__':

  parser = argparse.ArgumentParser(
    description="Calculate the branch distance between two haplogroups based on a tree file.\n"
                "----\n")

  parser.add_argument("-i1", "--inFile1",
                      help="path of file1 with haplogroup labels",
                      required=True,
                      type=str,
                      default="")
                      
  parser.add_argument("-c1", "--colSTR1",
                      help="column numbers of ID and hGrp in file 1 (default 0:1)",
                      type=str,
                      default="0:1")
                      
  parser.add_argument("-i2", "--inFile2",
                      help="path of file2 with haplogroup labels",
                      required=True,
                      type=str,
                      default="")
                      
  parser.add_argument("-c2", "--colSTR2",
                      help="column numbers of ID and hGrp in file 2 (default 0:1)",
                      type=str,
                      default="0:1")
                      
  parser.add_argument("-o", "--outFile",
                      help="path of output file",
                      required=True,
                      type=str,
                      default="")
                      
  parser.add_argument("-t", "--hGrpTreeFile",
                      help="path of file information all branches in the haplogroup tree",
                      required=True,
                      type=str,
                      default="")

                      
  args = parser.parse_args()


  main(args)

