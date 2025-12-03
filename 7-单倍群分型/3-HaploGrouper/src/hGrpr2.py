#!/usr/bin/env python

import sys
import argparse
#from time import perf_counter, time
#from time import perf_counter
#from time import perf_counteror as tick
#from time import process_time as timer
from time import perf_counter
import os
import numpy as np
import hGrpUtil


    

def main(args):  

  startTime = perf_counter()

  sys.stdout.write("\n############ Running haplogrouper ############\n")
  sys.stdout.write("vcf file: %s\n" % (args.vcfFile))
  sys.stdout.write("Haplogroup assignment will be based on\ntreeFile: %s\nlocusFile: %s\n" % (args.hGrpTreeFile, args.hGrpLocusFile))
  
  if os.path.isfile(args.vcfFile) == False:
    sys.stdout.write("Aborting! vcfFile not found: %s\n" % (args.vcfFile))
    return 0
    
  if os.path.isfile(args.hGrpTreeFile) == False:
    sys.stdout.write("Aborting! hGrpTreeFile not found: %s\n" % (args.hGrpTreeFile))
    return 0
  
  if os.path.isfile(args.hGrpLocusFile) == False:
    sys.stdout.write("Aborting! hGrpLocusFile not found: %s\n" % (args.hGrpLocusFile))
    return 0
    
  if args.referenceFasta != "":
    sys.stdout.write("GTs for variants not in vcfFile inferred as reference allele based on reference sequence: %s\n" % (args.referenceFasta))
    sys.stdout.write("If the vcf file contains genotypes for only a subset of positions, then you specify the regions in question using the -r option. Otherwise haplogroup assignments are liable to be biased towards the reference.\n")
    if os.path.isfile(args.referenceFasta) == False:
      sys.stdout.write("Aborting! file not found: %s\n" % (args.referenceFasta))
      return 0
    else:
      sys.stdout.write("Reading reference sequence: %s\n" % (args.referenceFasta))
      seqNameList, seqList, maxSeqLen = hGrpUtil.readFastaFile(args.referenceFasta)
      refSeq = seqList[0]
      refSeqName = seqNameList[0]
      del seqNameList, seqList
    
      
      
  sys.stdout.write("Results will be written to file: %s\n" % (args.outFile))
  
  # check whether limitations have been placed on the positions to be processed - these will be applied to the haplogroup locus file
  regionList = []
  if args.regions != "":
    # used to limit the processing of positions from vcfFile and refSeq
    regionList, regionCnt = hGrpUtil.parseSeqRegion(args.regions)
    regionListLen = len(regionList)
    if regionCnt > regionListLen:
      sys.stderr.write("Aborting! Not all regions sepcified with -r %s were valid\n" % (args.regions))
      return 0
      
    regPosSum = 0
    for (startPos, stopPos) in regionList:
      regPosSum += (stopPos + 1) - startPos
    sys.stdout.write("Processing will be limited to %d positions in %d genomic region(s) based on %s\n" % (regPosSum, regionListLen, args.regions))
  else:
    regionList = [(0, 1e20)]
    regionListLen = 1
    sys.stdout.write("All positions from haplogroup locus file will be used in the analysis\n")


  sys.stdout.write("\n")

  # process haplogroup tree and branch mutations
  sys.stdout.write("## Reading haplogroup tree and locus files\n")
  hGrpTree = hGrpUtil.readHGrpTreeFile(args.hGrpTreeFile)
  hGrpUtil.readLocusFile(args.hGrpLocusFile, hGrpTree, regionList)
  hGrpUtil.pruneHGrpTree(hGrpTree)
  sys.stdout.write("Read %d mutations at %d positions tagging %d haplogroup labels from tree and snp files\n" % (hGrpTree.mutCnt, hGrpTree.posCnt, hGrpTree.nodeUseCnt))

  
  ### fix the weights stuff
  if args.weightFile != "":
    sys.stdout.write("\n## Reading weights from file: %s\n" % (args.weightFile))
    if os.path.isfile(args.weightFile) == False:
      sys.stdout.write("Aborting! weightFile not found: %s\n" % (args.weightFile))
      return 0
    weightStream = open(args.weightFile)
    weightCnt = 0
    for lineNr, line in enumerate(weightStream):
      f = line.rstrip('\r\n').split('\t')
      pos = int(f[0])
      ancAl = f[1]
      derAl = f[2]
      weight = float(f[3])
      if pos in hGrpTree.mutPosDict:
        for mutIdx in hGrpTree.mutPosDict[pos]:
          if hGrpTree.mutList[mutIdx].derAl == derAl and  hGrpTree.mutList[mutIdx].ancAl == ancAl:
            hGrpTree.mutList[mutIdx].weight = weight
            weightCnt += 1
    sys.stdout.write("Processed a total of %d weights for mutations. They overlapped %d times with mutations in the haplogroup tree\n" % (lineNr, weightCnt))
    weightStream.close()
  
  

  ## determine which IDs from vcfFile will be processed
  if args.IDListFile != "":
    # read user-defined list of IDs to process
    userIDList = hGrpUtil.getListFromFile(args.IDListFile)
  else:
    userIDList = []  
  
  sys.stdout.write("\n## Scoring based on vcfFile\n")
  vcfStream, vcfHeader = hGrpUtil.openVCFStream(args.vcfFile, userIDList)
  #maxGTCnt = float(vcfHeader.IDUseIdxListLen)

  if args.IDListFile != "":
    sys.stdout.write("Found %d IDs that overlap between %d from IDListFile and %d from VCF file\n" % (vcfHeader.IDUseIdxListLen, len(userIDList), vcfHeader.IDListLen))
    if vcfHeader.IDUseIdxListLen == 0:
      sys.stdout.write("None of the IDs in IDListFile found in VCF file\nAborting!\n")
      return
  else:
    sys.stdout.write("Using all %d IDs from VCF file\n" % (vcfHeader.IDListLen))
  
  hGrpInMtx = np.zeros((hGrpTree.nodeCnt, vcfHeader.IDUseIdxListLen), dtype=float)
  hGrpOutMtx = np.zeros((hGrpTree.nodeCnt, vcfHeader.IDUseIdxListLen), dtype=float)

  # write code to check matching and mis-matching alleles for specific haplogroups
  posGTDict = {}
   
  buffPosCnt = 0
  buffPosDict = {}
  buffPosList = []
  buffGTList = []
  hGrpNrUsedList = []

  motifMtx = [[] for i in range(vcfHeader.IDUseIdxListLen)]

  sys.stderr.write("Processing GTs from vcfFile[. for every 1000 loci]")
  for lineNr, line in enumerate(vcfStream):
    line = line.rstrip('\r\n')
    vcfGT = hGrpUtil.procGTLine(line, vcfHeader)
    if (args.chrom != "" and vcfGT.chrom == args.chrom) or args.chrom == "":
      # check if position is from the correct chromosome - if user-defined chromosome criterion is defined args.chrom
      startPos = int(vcfGT.pos)
      stopPos = startPos + (len(vcfGT.refAl) - 1)
      if args.regions == "" or (any(lower <= startPos <= upper for (lower, upper) in regionList) == True and any(lower <= stopPos <= upper for (lower, upper) in regionList) == True):
        # if args.region has been specified, then check is locus falls within the specified range of positions
            
        '''
        Check for overlap of positions covered by reference allele with those in hGrpPosDict. 
        # This is complicated slightly by the fact that reference alleles can cover >1 position - so check overlap for each position covered by reference allele from single line of the vcfFile. 
        This ensures that we deal appropriately with indels that might overlap other SNPs (lines) in vcfFile
        '''
        usePosIdxList = []
        for i, refBase in enumerate(vcfGT.refAl):
          pos = startPos + i
          if pos in hGrpTree.mutPosDict:
            usePosIdxList.append([i, pos])
          
        if len(usePosIdxList) > 0:
          vcfGT = hGrpUtil.procGTSTRList(vcfGT, vcfHeader)
          gtArr = np.array([tmpList[0] for tmpList in vcfGT.gtList])
          for i, pos in usePosIdxList:
            posAlArr = []
            for a in vcfGT.alList[:-1]:
              if len(a) > i:
                posAlArr.append(a[i])
              else:
                posAlArr.append("")
            posAlArr.append("N")  # Agnar added this
            posAlArr = np.array(posAlArr)
        
            if pos in buffPosDict:
              # current position overlaps with a previous locus that had a reference allele with length >1
              buffPosIdx = buffPosDict[pos]
              gtIdxUpdateArr = np.where((gtArr>0) & (gtArr<vcfGT.missAL))[0]  # find the indexes of individuals who have non-reference alleles in gtArr
              #print "######", vcfGT.missAL, gtIdxUpdateArr
              buffGTList[buffPosIdx][gtIdxUpdateArr] = posAlArr[gtArr[gtIdxUpdateArr]]
            else:
              # current position does not overlap with a previous locus that had a reference allele with length >1
              buffPosDict[pos] = buffPosCnt
              buffGTList.append(posAlArr[gtArr])
              buffPosList.append([pos, posAlArr[0]])
              buffPosCnt += 1
        
        ## now remove positions <= startPos from buffer
        keepPosIdxBuffList = []
        remPosIdxBuffList = []
        
        for i, (pos, refAl) in enumerate(buffPosList):
          if pos <= startPos:
            # check this if statement
            remPosIdxBuffList.append(i)
            ###go through buffGTList and update 
            
            '''
            Need to take back mutations into account - particulary for mtDNA 
            for example - on branch L1'2'3'4'5'6 there are two mutations C146T and C182T that mutate back later for most haplotypes (T182C on L3'4) - but they get against scores for these branches
            The haplogroups assignments should not 
            
            could be done ahead of time - when processing the hGrpLocusFile - could assign all muts to each node - but this would entail much repitition 
            or here - by updating hGrpCntForMtx and hGrpCntAgainstMtx here instead of later - i.e. doing the path tracing on a locus by locus basis. Might slow things down. 
            '''
            for mutIdx in hGrpTree.mutPosDict[pos]:
              mut = hGrpTree.mutList[mutIdx]
              hGrpNrUsedList.append(mut.nodeNr)
              hGrpInMtx[mut.nodeNr] += (buffGTList[i] == mut.derAl) * mut.weight
              hGrpOutMtx[mut.nodeNr] += (buffGTList[i] == mut.ancAl) * mut.weight
            posGTDict[pos] = buffGTList[i]
          else:      
            keepPosIdxBuffList.append(i)
        
        if len(remPosIdxBuffList) > 0:
          # refresh info in position buffer variables
          buffGTList = [buffGTList[i] for i in keepPosIdxBuffList]
          buffPosList = [buffPosList[i] for i in keepPosIdxBuffList]
          buffPosCnt = len(buffPosList)
          buffPosDict = {}
          for i, (pos, refAl) in enumerate(buffPosList):
            buffPosDict[pos] = i  
        
    if lineNr % 1000 == 0:
      sys.stderr.write('.')      
      
  vcfStream.close()
  usedPosCnt = len(posGTDict)
  
  posGTDictLen = len(posGTDict)
  sys.stdout.write("Done in %.2f seconds.\nUsed GTs from %d of %d positions listed in hGrpLocusFile (%d in the vcfFile)\n" % (perf_counter()-startTime, usedPosCnt, hGrpTree.posCnt, lineNr))
  
  if usedPosCnt == 0:
    sys.stdout.write("No overlap of positions in the vcfFile and those in the hGrpLocusFile. Aborting!\n")
    return 0
  
  
  #print "posGTDictLen:", posGTDictLen
  #print posGTDict[73]
  

  if args.referenceFasta != "":
    '''
    Use information for positions that do not occur in the vcfFile
    - based on assumption that those positions were sequenced and all individuals carry the reference allele such positions
    '''
    #seqNameList, seqList, maxSeqLen = hGrpUtil.readFastaFile(args.referenceFasta)
    #refSeq = seqList[0]
    sys.stdout.write("\n## Scoring based on reference sequence [%s] of length %d\n" % (refSeqName, maxSeqLen))
    #del seqNameList, seqList
    for pos in hGrpTree.mutPosDict:
      # it is assumed that the reference sequence is from the correct chromosome. If not there will obviously be some problems!
      if pos not in posGTDict:
        # if position is in hGrpLocusFile and not in vcfFile - then set all GTs to reference allele
        refAl = refSeq[pos - 1]
        posGTDict[pos] = np.full(vcfHeader.IDUseIdxListLen, refAl)
        for mutIdx in hGrpTree.mutPosDict[pos]:
          mut = hGrpTree.mutList[mutIdx]
          hGrpNrUsedList.append(mut.nodeNr)
          if refAl == mut.derAl:
            hGrpInMtx[mut.nodeNr] += mut.weight
          elif refAl == mut.ancAl:
            hGrpOutMtx[mut.nodeNr] += mut.weight

    refSeqLocCnt = len(posGTDict) - usedPosCnt  
    sys.stdout.write("GTs were inferred for %d additional loci not reported in vcfFile\n" % (refSeqLocCnt))
          
  
  
  hGrpNrUsedList = list(set(hGrpNrUsedList))
  # now make sure all the ancestral nodes are included in this list
  hGrpNrUsedWAncCntArr = np.zeros(hGrpTree.nodeCnt, dtype=int)
  for idxHGrp in hGrpNrUsedList:
    hGrpNrUsedWAncCntArr[hGrpTree.nodeList[idxHGrp].pathArr] += 1
  hGrpNrUsedWAncArr = np.where(hGrpNrUsedWAncCntArr > 0)[0]
  
  sys.stdout.write("\n## Making assignments based on %d of %d haplogroup labels that were encountered for scoring\n" % (hGrpNrUsedWAncArr.size, hGrpTree.nodeUseCnt))


  # calculate for, against and net difference matrices - but only for hGrpNrUsedList (and also their ancestral nodes?)
  hGrpCntForMtx = np.zeros((hGrpTree.nodeCnt, vcfHeader.IDUseIdxListLen), dtype=float)
  hGrpCntAgainstMtx = np.zeros((hGrpTree.nodeCnt, vcfHeader.IDUseIdxListLen), dtype=float)
  #for idxHGrp in hGrpNrUsedList:
  for idxHGrp in hGrpNrUsedWAncArr:
    hGrpCntForMtx[idxHGrp] += np.sum(hGrpInMtx[hGrpTree.nodeList[idxHGrp].pathArr,:], axis=0)
    hGrpCntAgainstMtx[idxHGrp] += np.sum(hGrpOutMtx[hGrpTree.nodeList[idxHGrp].pathArr,:], axis=0)
  hGrpCntDiffMtx = hGrpCntForMtx - hGrpCntAgainstMtx


  # get topN highest net difference scores 
  topN = 10
  kPos = hGrpTree.nodeCnt - topN
  hGrpCntDiffTop10Mtx = np.partition(hGrpCntDiffMtx.T, kPos)[:,-topN:]
  maxDiffArr = np.max(hGrpCntDiffTop10Mtx, axis=1)


  ## write outFile
  outStream = open(args.outFile,'w')
  outStream.write("ID\tHaplogroup\tnetScore\tmatchScore\tmismatchScore\tmismatchLoci\tbackMutCnt\tpruning\tallMaxNetScore\n")
  for i, indIdx in enumerate(vcfHeader.IDUseIdxList):
    pn = vcfHeader.IDList[indIdx]
    idxHGrpArr = np.where((hGrpCntDiffMtx[:,i] == maxDiffArr[i]) & (hGrpInMtx[:,i] > 0))[0]
    # determine the minimum number of loci that contradict the haplogroup assignment
    maxDiffHGrpList = ["%s[%.5g-%.5g]" % (hGrpTree.nodeList[idxHGrp].name, hGrpCntForMtx[idxHGrp, i], hGrpCntAgainstMtx[idxHGrp, i]) for idxHGrp in idxHGrpArr]
    
    pruneList = ["",""]
    minAgainst = np.min(hGrpCntAgainstMtx[idxHGrpArr,i])
    idxHGrpFinalArr = idxHGrpArr[hGrpCntAgainstMtx[idxHGrpArr, i] == minAgainst]
    if idxHGrpArr.size > idxHGrpFinalArr.size:
      # at least one haplogroup was removed from the list because it had > minAgainst
      pruneList[0] = "A"
    if idxHGrpFinalArr.size > 1:
      # still >1 node has highest match score - find their mrca
      pruneList[1] = "M"
      idxHGrpFinalArr = hGrpUtil.getMRCAnode2(idxHGrpFinalArr, hGrpTree)
    idxHGrpFinal = idxHGrpFinalArr[0]
    #print "################ idxHGrpFinal", pn, idxHGrpFinal, hGrpPathList[idxHGrpFinal][0]
    ## identify the GTs that mismatch the mutations underlying the haplogroup assignment
    mutMisMatchList = []
    backMutList = []
    tmpDict = {}
    
    for nodeNr in np.flip(hGrpTree.nodeList[idxHGrpFinal].pathArr,axis=0):
      for mutIdx in hGrpTree.nodeList[nodeNr].mutIdxList:
        mut = hGrpTree.mutList[mutIdx]
        if mut.pos in posGTDict:
          gt = posGTDict[mut.pos][i]
          if gt != mut.derAl and gt == mut.ancAl:
            # GT is mismatch with mutation on branch leading to this node
            outSTR = "%s[%s]%s%d%s" % (gt, hGrpTree.nodeList[mut.nodeNr].name, mut.ancAl, mut.pos, mut.derAl)
            if mut.pos not in tmpDict:
              mutMisMatchList.append(outSTR)
            else:  
              backMutList.append(outSTR)
          if mut.pos not in tmpDict:
            tmpDict[mut.pos] = gt
      
    backMutListCnt = len(backMutList)
    outStream.write("%s\t%s\t%.5g\t%.5g\t%.5g\t%s\t%d\t%s\t%s\n" % (pn, hGrpTree.nodeList[idxHGrpFinal].name, hGrpCntDiffMtx[idxHGrpFinal, i], hGrpCntForMtx[idxHGrpFinal, i], hGrpCntAgainstMtx[idxHGrpFinal, i], ":".join(mutMisMatchList), backMutListCnt, "".join(pruneList), ":".join(maxDiffHGrpList)))
    #print pn, hGrpTree.nodeList[idxHGrpFinal].name, hGrpCntDiffMtx[idxHGrpFinal, i], hGrpCntForMtx[idxHGrpFinal, i], hGrpCntAgainstMtx[idxHGrpFinal, i], hGrpNrUsedArr[idxHGrpFinal]
  outStream.close()
    
  if args.verboseFile != "":
    verboseStream = open(args.verboseFile, 'w')  
    tmpSTR = "\t".join([vcfHeader.IDList[indIdx] for i, indIdx in enumerate(vcfHeader.IDUseIdxList)])
    verboseStream.write("hGrp\t%s\n" % (tmpSTR))
    for idxHGrp in hGrpNrUsedList:
      tmpSTR = "\t".join(["%.5g-%.5g" % (hGrpCntForMtx[idxHGrp, i], hGrpCntAgainstMtx[idxHGrp, i]) for i in range(vcfHeader.IDUseIdxListLen)])
      verboseStream.write("%s\t%d\t%s\n" % (hGrpTree.nodeList[idxHGrp].name, len(hGrpTree.nodeList[idxHGrp].mutIdxList), tmpSTR))

    verboseStream.close()
  
  sys.stderr.write("Finished in %.2f seconds\n" % (perf_counter()-startTime))
    


if __name__ == '__main__':

  parser = argparse.ArgumentParser(
    description="Determine haplogroup for list of individuals based on VCF file and haplogroup tree.\n"
                "----\n")

    
  parser.add_argument("-v", "--vcfFile",
                      help="path of vcfFile to be converted (can be gzipped vcf file)",
                      required=True,
                      type=str,
                      default="")
  
  parser.add_argument("-o", "--outFile",
                      help="path of output file",
                      required=True,
                      type=str,
                      default="")

  parser.add_argument("-l", "--hGrpLocusFile",
                      help="path of file information about loci used to assign individuals to haplogroups and the branches in the haplogroup tree",
                      required=True,
                      type=str,
                      default="")
                      
  parser.add_argument("-t", "--hGrpTreeFile",
                      help="path of file information all branches in the haplogroup tree",
                      required=True,
                      type=str,
                      default="")

  parser.add_argument("-i", "--IDListFile",
                      help="path of file with subset of IDs from VCF file that are to be used for haplogroup assignment",
                      type=str,
                      default="")
                      
  parser.add_argument("-r", "--regions",
                      help="Only base haplogroup assignment on positions from the specified regions (format: startPos-stopPos,startPos-stopPos) Multiple regions can be specified",
                      type=str,
                      default="")
                      
  parser.add_argument("-c", "--chrom",
                      help="Only process loci from vcf file with this chromosome name. This must be used when the vcf file contains loci from multiple chromosomes.",
                      type=str,
                      default="")
                      
  parser.add_argument("-f", "--referenceFasta",
                      help="Read reference sequence from fasta file. This is useful when the vcfFile is based on full sequence data, but only reports polymorphic positions or differences from the reference sequence. In such cases, many phylogenetically informative positions would be ignored. When the reference sequence is provided, positions not reported in the vcfFile are assumed to have the reference state for all individuals in the file.",
                      type=str,
                      default="")    

  parser.add_argument("-m", "--mismatchLoc",
                      help="Report genotypes of mismatching loci in outFile",
                      default=False,
                      action="store_true")

  parser.add_argument("-w", "--weightFile",
                      help="Give mutations differing weights read from a separate file. This is only useful for loci with a high rate of recurrent mutations - like the mtDNA control region. The file should be tab-delimited with four columns: pos ancAl derAl weight. Mutations not included in the file default to a weight of 1",
                      type=str,
                      default="")  

  parser.add_argument("-x", "--verboseFile",
                      help="path of file for full matrix of scores for each node in the tree (lines) for each individual (columns)",
                      type=str,
                      default="")                        
                      
  args = parser.parse_args()


  main(args)

