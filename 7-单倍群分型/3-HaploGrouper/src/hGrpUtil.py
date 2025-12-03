#!/usr/bin/env python


import gzip
import sys
import numpy as np
import os

class vcfHeaderClass:
  def __init__(self):
    self.firstIDColNr = None
    self.IDList = []  # all IDs in the vcfFile
    self.IDListLen = 0
    self.IDUseIdxList = [] # indexes in IDlist of the IDs for whom GTs will be examined
    self.extIDListIdxList = [] # indexes in user supplied IDlist - same individuals and in same order as IDUseIdxList
    self.IDUseIdxListLen = 0 
    self.headerList = []
    self.infoList = []
    self.formatList = []
    self.otherList = []
    self.valTypeDict = {"ID":0, "Number":1, "Type":2, "Description":3}
    self.headerLineCnt = 0
    # entries in infoDict and formatDict are key: ID and value: [Number, Type, Description]


class vcfGTClass:
  def __init__(self):
    self.chrom = None
    self.pos = None
    self.refAl = ""
    self.ID = ""
    self.qual = 0.0
    self.infoSTR = ""
    self.filtList = []
    self.GTColNr = None
    self.ADColNr = None
    self.GQColNr = None
    self.missAL = None
    self.fGTList = [] # FORMAT
    self.missGTCnt = 0 
    self.alList = []
    self.indel = 0
    self.alCnt = 0
    self.gtList = [] # GTs
    #self.gtcharList = [] # GTs string 
    self.gtSTRList = [] # list of unparsed strings with genotype info for individuals
    self.totRCnt = 0
    self.rDepthAl = []
    self.rDepthIndSum = []
    self.rDepthAlSum = []


class hGrpMutClass:
  def __init__(self):
    self.pos = None
    self.name = ""
    self.derAl = ""
    self.nodeNr = None
    self.ancAl = ""
    self.weight = 1.0
    
class hGrpNodeClass:
  def __init__(self):
    self.name = ""
    self.parentName = ""
    self.pathArr = None
    self.idx = None
    self.parentIdx = None
    self.depth = 0
    self.mutIdxList = []

class hGrpClass:
  def __init__(self):
    self.nodeDict = {}
    self.nodeList = []
    self.nodeUseList = []
    self.nodeDropList = []
    self.nodeCnt = None
    self.mutPosDict = {}
    self.mutList = []
    self.mutCnt = None
    self.posCnt = None


def openVCF(vcfFile, zip=0):
  '''
  open vcfFile and return vcfStream. If zip = 1, then open with gzip
  '''  
  if zip == 1:
    vcfStream = gzip.open(vcfFile,'r')
  else:
    vcfStream = open(vcfFile, 'r')
  return vcfStream



  
def openVCFStream(vcfFile, useIDList=[], verbose=True, idDelim=""):
  # check if vcf file is gzipped, open and return vcfStream and vcfHeader. Specification of useIDList is optional
  try:
    # try opening vcf files as gzipped  
    vcfStream = openVCF(vcfFile, 1)
    #print "trying as zipped"
    vcfHeader = checkVCFFile(vcfStream, useIDList, idDelim)
    if verbose == True:
      sys.stderr.write("VCF file is gzipped\n")
  except:  
    vcfStream = openVCF(vcfFile, 0)
    #print "trying as unzipped"
    vcfHeader = checkVCFFile(vcfStream, useIDList, idDelim)
    if verbose == True:
      sys.stderr.write("VCF file is not gzipped\n")

  return vcfStream, vcfHeader

  
def procGTLine(GTLine, vcfHeader):
  '''
  Initial processing of line with genotype data from vcf - returns locus specific information and list of genotype strings
  # CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT
  #chr21	5030736	.	T	C,<NON_REF>	13.95	.	AC=1,0;AN=2;BaseQRankSum=0.135;ClippingRankSum=0.807;DP=16;MLEAC=1,0;MLEAF=0.500,0.00;MQ0=0;MQ=43.42;MQRankSum=-1.614;ReadPosRankSum=-0.942;SF=6	GT:GQ:DP:PL:SB:AD	.	.	.	.	.	.	0/1:42:16:42,0,421,80,430,510:4,9,1,2:13,3,0	.	.	.	.	.
  # chr21	5030912	.	CA	C,<NON_REF>	0.01	.	AC=1,0;AN=2;BaseQRankSum=-1.613;ClippingRankSum=-0.594;DP=15;MLEAC=1,0;MLEAF=0.500,0.00;MQ0=0;MQ=35.23;MQRankSum=1.104;ReadPosRankSum=0.639;SF=8	GT:GQ:DP:PL:SB:AD	.	.	.	.	.	.	.	.	0/1:12:15:12,0,314,51,320,372:9,4,1,1:13,2,0	.	.	.

  might want to write a more versatile version that does not parse all the info from the locus - used for filtering 
  '''

  vcfGT = vcfGTClass()
  f = GTLine.split('\t')
  vcfGT.chrom = f[0]  
  vcfGT.pos = f[1]  
  vcfGT.ID = f[2]
  
  vcfGT.refAl = f[3]
  if vcfGT.refAl == ".":
    vcfGT.refAl = "N"
  vcfGT.alList = [vcfGT.refAl]
  tmpList = f[4].split(',')
  for al in tmpList:
    if al != "<NON_REF>":
      if al != ".":
        vcfGT.alList.append(al)
      else:  
        vcfGT.alList.append("N")
  # check for indels
  for al in vcfGT.alList:
    if len(al) > 1:
      vcfGT.indel = 1

  try:
    vcfGT.qual = float(f[5])  
  except:
    vcfGT.qual = None
  filtSTR = f[6]
  vcfGT.infoSTR = f[7]
  if filtSTR != ".":
    vcfGT.filtList = [filt for filt in filtSTR.split(';')]
  vcfGT.alCnt = len(vcfGT.alList)
  vcfGT.missAL = vcfGT.alCnt
  vcfGT.alList.append("N") # missing allele

  # find GT and AD in genotype string
  vcfGT.fGTList = f[vcfHeader.firstIDColNr - 1].split(':')
  try:
    vcfGT.GTColNr = vcfGT.fGTList.index("GT")
  except:
    #print "no GT present in line!"
    vcfGT.GTColNr = None

  try:
    vcfGT.ADColNr = vcfGT.fGTList.index("AD")
  except:
    #print "no AD present in line!"
    vcfGT.ADColNr = None

  try:
    vcfGT.GQColNr = vcfGT.fGTList.index("GQ")
  except:
    #print "no GQ present in line!"
    vcfGT.GQColNr = None

  vcfGT.gtSTRList = f[vcfHeader.firstIDColNr:]

  return vcfGT


def procGTSTRList(vcfGT, vcfHeader):
  '''
  Extract genotypes from unparsed list of genotype strings from vcf genotype line (vcfGT.gtSTRList)
  Example: [GT:AD:DP:GQ:PL]
  0/0:98,6:107:14:0,14,438	0/1:80,5:85:64:64,0,296	0/0:96,8:104:33:0,33,377	0/1:71,6:77:13:13,0,228	0/1:108,14:123:22:22,0,433	0/0:73,5:78:39:0,39,440	0/0:73,5:78:7:0,7,368	0/0:124,7:131:33:0,33,661	0/0:83,8:93:14:0,14,491	0/0:70,4:74:8:0,8,419	0/0:82,5:88:60:0,60,741	0/0:132,9:142:27:0,27,538
  ./.:.:.:.:.	0/0:11,1:12:6:0,6,63	0/0:2,0:2:3:0,3,22	./.:.:.:.:.	0/0:14,3:18:15:0,15,194	0/0:6,0:6:9:0,9,113	0/0:7,0:9:6:0,6,73	0/1:22,3:26:99:103,0,122	1/1:5,2:7:3:26,3,0	./.:.:.:.:.	0/0:2,0:3:3:0,3,26	0/0:13,1:15:15:0,15,183
  GT should always be present
  '''

  vcfGT.gtList = [[] for i in range(vcfHeader.IDUseIdxListLen)]
  if vcfGT.ADColNr != None:
    #vcfGT.rDepthList = [[] for i in xrange(vcfHeader.IDListLen)]
    #vcfGT.rDepthList = np.array([[0 for j in xrange(vcfHeader.IDListLen)] for i in xrange(vcfGT.alCnt)])
    vcfGT.rDepthAl = np.zeros((vcfHeader.IDUseIdxListLen, vcfGT.alCnt), dtype=int)


  for i, idx in enumerate(vcfHeader.IDUseIdxList):    
  #for idx, vcfGTSTR in enumerate(vcfGT.gtSTRList):
    vcfGTSTR = vcfGT.gtSTRList[idx]
    vcfGTList = vcfGTSTR.split(':')
    gtSTR = vcfGTList[vcfGT.GTColNr]
    vcfGT.gtList[i] = parseGT(gtSTR, vcfGT.missAL)
    if vcfGT.gtList[i] == vcfGT.missAL:
      vcfGT.missGTCnt += 1
    if vcfGT.ADColNr != None:
      # AD - read counts for alleles
      ADstr = vcfGTList[vcfGT.ADColNr]
      if ADstr != ".":
        ADList = ADstr.split(',')[:vcfGT.alCnt]
        for alNr, rCntSTR in enumerate(ADList):
          #vcfGT.rDepthAl[alNr][idx] = int(rCntSTR)
          vcfGT.rDepthAl[i][alNr] = int(rCntSTR)

          #except:
      #  print vcfGT.gtSTRList
      #  print vcfGT.pos, vcfGT.fGTList, vcfGTSTR, vcfGTList
      #  print "Error", vcfGT.GTColNr, gtSTR, vcfGT.ADColNr, ADstr

  return vcfGT
  

def calcrDepthSum(vcfGT):
  vcfGT.rDepthIndSum = np.sum(vcfGT.rDepthAl, axis=1, dtype=float)
  vcfGT.rDepthAlSum = np.sum(vcfGT.rDepthAl, axis=0, dtype=float)
  vcfGT.totRCnt = np.sum(vcfGT.rDepthAlSum)

def checkRDepthThresh(vcfGT, rDepthThresh):
  # check whether all individuals have read depth >= rDepthThresh 
  passThresh = True 
  calcrDepthSum(vcfGT)
  rdFailCnt = (vcfGT.rDepthIndSum < rDepthThresh).sum()
  if rdFailCnt > 0:
    #print "fail rd"
    passThresh = False
  return passThresh
  
  
def parseGT(gtSTR, missAL):
  '''
  gtSTR is string with format [a1][delim][a2], where delim is either [/] or [|]
  a1 and a2 can be either integer, where 0 is reference allele, 1 is first alternative allele, etc.
  or a1 and a2 can be [.] - indicating missing allele. What is the integer for missing? 0?
  '''
  ploid = 0
  phased = 0

  gtList = gtSTR.split('/')
  ploid = len(gtList)
  if ploid == 1:
    # genotype is not diploid unphased, check if diploid phased
    gtList = gtSTR.split('|')
    ploid = len(gtList)
    if ploid > 1:
      # genotype is diploid phased
      phased = 1
  #[x+1 if x >= 45 else x+5 for x in l]   
  return [int(a) if a != "." else missAL for a in gtList]



def procHeaderLine(headerLine, vcfHeader):
  ##INFO=<ID=ALL_AF,Number=1,Type=Float,Description="ALT Allele Frequency for samples from all population together">
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">
  ##fileformat=VCFv4.1
  # example of info string: VT=SNP;SNPSOURCE=LOWCOV;ALL_AF=0.0078;AFR_AF=0.0020;AMR_AF=0.0083;ASN_AF=0.0175;EUR_AF=0.0040
  if headerLine[:8] == "##FORMAT":
    #print headerLine
    tmpList = procFormatInfoLine(headerLine, vcfHeader.valTypeDict)
    vcfHeader.formatList.append(tmpList)
    #print "format:", 
  elif headerLine[:6] == "##INFO":
    #print headerLine
    tmpList = procFormatInfoLine(headerLine, vcfHeader.valTypeDict)
    vcfHeader.infoList.append(tmpList)
  else:
    vcfHeader.otherList.append(headerLine.lstrip('#'))


def procFormatInfoLine(headerLine, valTypeDict):
  '''
  Cannot currently handle: ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
  Should only have 4 
  '''
  headerLine = headerLine.lstrip('#')[:-1]
  #print "******\n", headerLine
  #print headerLine, len(headerLine.split('=<'))
  headType, propSTR = headerLine.split('=<')
  propList = propSTR.split(',')
  propListLen = len(propList)
  if propListLen > 4:
    propList[3] = ",".join(propList[3:])
    propList = propList[:4]
  outList = [None for i in range(4)]
  for prop in propList:
    #print prop
    #name, val = prop.split('=')
    propList = prop.split('=')
    name = propList[0]
    val = "".join(propList[1:])
    val = val.strip('"')
    if name in valTypeDict:
      valIdx = valTypeDict[name]
      outList[valIdx] = val
    else:
      sys.stderr.write("%s not recognized from vcf header" % (prop))
  return outList


def checkVCFFile(vcfStream, IDUseList=[], idDelim=""):
  # check for vcf format and return list of IDs from file
  vcfBool = True
  IDList = []

  #deal with file structure - header info. Should I store this in a dictionary?
  vcfHeader = vcfHeaderClass()

  varInfoBool = True
  while varInfoBool == True:
    line = vcfStream.readline()
    line = line.strip()
    vcfHeader.headerLineCnt += 1
    if line[:2] == "##":
      #print line
      procHeaderLine(line, vcfHeader)
    else:
      varInfoBool = False
    # should add code to properly parse and use this information
  
  # deal with header - CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  HG00096  HG00097  HG00099  HG00100  ... and so on
  # note - INFO column is optional
  line = line.lstrip('#')
  headerList = line.split('\t')
  vcfHeader.firstIDColNr = headerList.index("FORMAT") + 1

  if len(headerList) <= vcfHeader.firstIDColNr:
    vcfBool = False
  if headerList[:vcfHeader.firstIDColNr] != ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]:
    vcfBool = False

  if vcfBool == False:
    sys.stderr.write("%s does not appear to be a VCF file with genotypes" % (vcfFile))
  else:
    IDList = headerList[vcfHeader.firstIDColNr:]
    if idDelim != "":
      IDList = [id.split(idDelim)[0] for id in IDList]
    #print line
  vcfHeader.headerList = headerList
  vcfHeader.IDList = IDList
  vcfHeader.IDListLen = len(IDList)
  if len(IDUseList) != 0:
    # specified subset of IDs to be used
    IDUseDict = {}
    for i, ID in enumerate(IDUseList):
      IDUseDict[ID] = i
    #IDUseSet = set(IDUseList)
    for i, ID in enumerate(vcfHeader.IDList):
      if ID in IDUseDict:
        vcfHeader.IDUseIdxList.append(i)
        IDUSeNr = IDUseDict[ID]
        vcfHeader.extIDListIdxList.append(IDUSeNr)
    vcfHeader.IDUseIdxListLen = len(vcfHeader.IDUseIdxList)
  else:      
    # all IDs to be used
    vcfHeader.IDUseIdxList = [i for i in range(vcfHeader.IDListLen)]
    vcfHeader.IDUseIdxListLen = vcfHeader.IDListLen

    #print vcfHeader.IDUseIdxList
  #print "formatList", vcfHeader.formatList
  #print "infoList", vcfHeader.infoList
  #print "otherList", vcfHeader.otherList
  #return IDList, vcfHeader
  return vcfHeader


def makeInfoDict(vcfHeader):
  ##INFO=<ID=ALL_AF,Number=1,Type=Float,Description="ALT Allele Frequency for samples from all population together">
  # example of info string: VT=SNP;SNPSOURCE=LOWCOV;ALL_AF=0.0078;AFR_AF=0.0020;AMR_AF=0.0083;ASN_AF=0.0175;EUR_AF=0.0040
  infoDict = {}
  for infoItem in vcfHeader.infoList:
    infoDict[infoItem[0]] = None
  return infoDict


def parseSeqRegion(regionSTR):
  # assumes input startPos-stopPos,startPos-stopPos
  
  regionList = []
  tmpList = regionSTR.split(',')
  tmpListLen = len(tmpList)
  for tmpSTR in tmpList:
    f = tmpSTR.split("-")
    flen = len(f)
    if flen == 2:
      startPos = int(f[0])
      stopPos = int(f[1])
      if stopPos >= startPos:
        regionList.append((startPos, stopPos))
      else:
        sys.stderr.write("Position range %s does not make sense and will be ignored. startPos must be smaller than stopPos.\n" % (tmpSTR))      

  return regionList, tmpListLen

 

def getListFromFile(inFile, delim = "\t", colNr=0):
  # if multiple columns (tab-delimited) - then read from first column
  outList = []
  if os.path.isfile(inFile):
    inStream = open(inFile)
    for line in inStream:
      line = line.strip('\r\n')
      f = line.split(delim)
      outList.append(f[colNr])
    inStream.close()  
  else:
    sys.stderr.write("File: %s does not exist!\n" % (inFile))
  return outList


 
def readHGrpPathFile(hGrpPathFile):
  hGrpPathDict = {}
  hGrpPathList = []
  '''
  # Format:
  # K1e	mt-MRCA:L1'2'3'4'5'6:L2'3'4'5'6:L2'3'4'6:L3'4'6:L3'4:L3:N:R:U:U2'3'4'7'8'9:U8:U8b:K:K1:K1_sub:K1e
  # H3k1	mt-MRCA:L1'2'3'4'5'6:L2'3'4'5'6:L2'3'4'6:L3'4'6:L3'4:L3:N:R:R0:HV:H:H3:H3_sub:H3k:H3k1
  '''
  inStream = open(hGrpPathFile)
  for hGrpPathCnt, line in enumerate(inStream):
    line = line.rstrip('\r\n')
    f = line.split('\t')
    hGrpPathList.append([f[0],f[1]])
    hGrpPathDict[f[0]] = hGrpPathCnt
  inStream.close()  
  hGrpPathCnt += 1

  return hGrpPathCnt, hGrpPathDict, hGrpPathList


 
def readHGrpTreeFile(hGrpPathFile):

  '''
  Tree file has two columns, where first column contains the name of a node and the second the name of the node's parent. 
  From this it is possible to derive the path of each node to the root of the tree (the only node with no parent)
  # Format:
  YRoot		0
  A0000	YRoot	1
  A000-T	YRoot	1
  A000	A000-T	2
  '''
  hGrpTree = hGrpClass()
  inStream = open(hGrpPathFile)
  for lineNr, line in enumerate(inStream):
    line = line.rstrip('\r\n')
    f = line.split('\t')
    hGrpNode = hGrpNodeClass()
    hGrpNode.name = f[0]
    hGrpNode.idx = lineNr
    hGrpNode.parentName = f[1]
    hGrpTree.nodeList.append(hGrpNode)
    hGrpTree.nodeDict[hGrpNode.name] = lineNr
  inStream.close()  
  hGrpTree.nodeCnt = len(hGrpTree.nodeList)

  ## get idx from hGrpTree.nodeDict for each parent and path to root
  for node in hGrpTree.nodeList:
    if node.parentName in hGrpTree.nodeDict:
      node.parentIdx = hGrpTree.nodeDict[node.parentName]
    else:
      sys.stderr.write("Node %s in line %d has parent %s that is not present in the tree. This should only happen for root node\n" % (node.name, node.idx, node.parentName))            
      node.parentIdx = -999
      
  ## get path to root for each node 
  for node in hGrpTree.nodeList:
    ancNode = hGrpTree.nodeList[node.idx]
    pathList = [node.idx]
    while ancNode.parentIdx >= 0:
      pathList.append(ancNode.parentIdx)
      ancNode = hGrpTree.nodeList[ancNode.parentIdx]
    node.depth = len(pathList) - 1  
    pathList.reverse()
    node.pathArr = np.array(pathList)
  return hGrpTree


 
def readLocusFile(hGrpLocusFile, hGrpTree, regionList): 
  '''
  ## read the locus file for haplogroup path file 
  Format:
  SNP	pos	forAl	hapGrp	againstAl
  G263A	263	A	L0	G
  
  Also, if the --region option was used, then prune all positions that do not intersect with the specified regions 

  
  hGrpNodeDict is a dictionary with hGrp as key and hGrpIdx as value derived from the function readHGrpPathFile
  '''
  #hGrpMutList = []
  #hGrpPosDict = {}
  #mutByHGrpList = [[] for hGrp in hGrpNodeDict]

  hGrpLocusStream = open(hGrpLocusFile, 'r')
  hdr = hGrpLocusStream.readline() # read and ignore header line
  hGrpMutCnt = 0
  for line in hGrpLocusStream:
    f = line.strip('\r\n').split('\t')
    if f[1].isdigit() == True:
      pos = int(f[1])
      if any(lower <= pos <= upper for (lower, upper) in regionList):
        nodeName = f[3]
        if nodeName in hGrpTree.nodeDict:
          nodeNr = hGrpTree.nodeDict[nodeName]
          # populate hGrpLocClass object
          tmpLoc = hGrpMutClass()
          tmpLoc.pos = pos
          tmpLoc.name = f[0]
          tmpLoc.derAl = f[2].upper()
          tmpLoc.nodeNr = nodeNr
          tmpLoc.ancAl = f[4].upper()
          hGrpTree.mutList.append(tmpLoc)
          if pos not in hGrpTree.mutPosDict:
            hGrpTree.mutPosDict[pos] = []
          hGrpTree.mutPosDict[pos].append(hGrpMutCnt)
          hGrpTree.nodeList[nodeNr].mutIdxList.append(hGrpMutCnt)
          hGrpMutCnt += 1
          
        else:
          sys.stderr.write("nodeName %s not found in path file. Should not happen!\n" % (nodeName))      
    else:
        sys.stderr.write("position %s is not numeric and will be ignored\n" % (f[1]))      
  hGrpLocusStream.close()   
  hGrpTree.mutCnt = len(hGrpTree.mutList)
  hGrpTree.posCnt = len(hGrpTree.mutPosDict)



def pruneHGrpTree(hGrpTree):
  '''
  make list of node indexes that have at least one mutation - only these nodes can be assigned to an individuals
  Nodes with no mutations will be removed from hGrpTree
  
  '''
  for nodeIdx, node in enumerate(hGrpTree.nodeList):
    if len(node.mutIdxList) > 0:
      hGrpTree.nodeUseList.append(nodeIdx)
    else:
      hGrpTree.nodeDropList.append(nodeIdx)
  hGrpTree.nodeUseCnt = len(hGrpTree.nodeUseList)
  dropListLen = len(hGrpTree.nodeDropList)
  if dropListLen == 0:
    sys.stderr.write("All %d nodes have at least one mutation\n" % (hGrpTree.nodeCnt))
  else:
    sys.stderr.write("%d of %d nodes have at least one mutation and will be used (%d have none and will be ignored)\n" % (hGrpTree.nodeUseCnt, hGrpTree.nodeCnt, dropListLen))
      
  

 

def readSNPFile(hGrpLocusFile, hGrpNodeDict): 
  '''
  ## read the locus file for haplogroup path file 
  Format:
  SNP	pos	forAl	hapGrp	againstAl
  G263A	263	A	L0	G
  
  hGrpNodeDict is a dictionary with hGrp as key and hGrpIdx as value derived from the function readHGrpPathFile
  '''
  hGrpMutList = []
  hGrpPosDict = {}
  mutByHGrpList = [[] for hGrp in hGrpNodeDict]

  hGrpLocusStream = open(hGrpLocusFile, 'r')
  hdr = hGrpLocusStream.readline() # read and ignore header line
  hGrpMutCnt = 0
  for line in hGrpLocusStream:
    f = line.strip('\r\n').split('\t')
    if f[1].isdigit() == True:
      pos = int(f[1])
      nodeName = f[3]
      if nodeName in hGrpNodeDict:
        nodeNr = hGrpNodeDict[nodeName]
        # populate hGrpLocClass object
        tmpLoc = hGrpMutClass()
        tmpLoc.pos = pos
        tmpLoc.name = f[0]
        tmpLoc.derAl = f[2].upper()
        tmpLoc.nodeNr = nodeNr
        tmpLoc.ancAl = f[4].upper()
        hGrpMutList.append(tmpLoc)
        if pos not in hGrpPosDict:
          hGrpPosDict[pos] = []
        hGrpPosDict[pos].append(hGrpMutCnt)
        mutByHGrpList[nodeNr].append(hGrpMutCnt)
        hGrpMutCnt += 1
        
      else:
        sys.stderr.write("nodeName %s not found in path file. Should not happen!\n" % (nodeName))      
    else:
        sys.stderr.write("position %s is not numeric and will be ignored\n" % (f[1]))      
  hGrpLocusStream.close()   
  
  return hGrpMutList, hGrpPosDict, mutByHGrpList

 


def getMRCAnode(idxHGrpFinalArr, hGrpPathAncList):
  ancDict = {}
  for idxHGrp in idxHGrpFinalArr:
    for idxAnc in hGrpPathAncList[idxHGrp]:
      if idxAnc not in ancDict:
        ancDict[idxAnc] = 0
      ancDict[idxAnc] += 1
  for idxAnc in np.flip(hGrpPathAncList[idxHGrpFinalArr[0]],axis=0):
    if ancDict[idxAnc] == idxHGrpFinalArr.size:
      idxHGrpFinalArr = np.array([idxAnc])
      break
  return idxHGrpFinalArr 



def getMRCAnode2(idxHGrpFinalArr, hGrpTree):
  ancDict = {}
  for idxHGrp in idxHGrpFinalArr:
    for idxAnc in hGrpTree.nodeList[idxHGrp].pathArr:
      if idxAnc not in ancDict:
        ancDict[idxAnc] = 0
      ancDict[idxAnc] += 1
  for idxAnc in np.flip(hGrpTree.nodeList[idxHGrpFinalArr[0]].pathArr,axis=0):
    if ancDict[idxAnc] == idxHGrpFinalArr.size:
      idxHGrpFinalArr = np.array([idxAnc])
      break
  return idxHGrpFinalArr 


  
def readFastaFile(fasFile):
  # assumes that muliple sequences are aligned
  fasStream = open(fasFile)
  seq = ""
  seqNameList = []
  seqList = []
  maxSeqLen = 0
  for line in fasStream:
    line = line.rstrip('\r\n')
    if line != "":
      if line[0] == ">":
        #seqname
        seqNameList.append(line[1:])
        if seq != "":
          seqList.append(seq)
          seqLen = len(seq)
          #print seqLen
          if seqLen > maxSeqLen:
            maxSeqLen = seqLen
          seq = ""
        else:
          if len(seqList) > 0:
            sys.stderr.write("%s has empty sequence!" % (seqNameList[-1]))
      else:
        #Process seq  
        seq = seq + line
  seqList.append(seq)
  seqLen = len(seq)

  if seqLen > maxSeqLen:
    maxSeqLen = seqLen
  fasStream.close()    
  return (seqNameList, seqList, maxSeqLen)
