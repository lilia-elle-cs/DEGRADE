import csv
from striprtf.striprtf import rtf_to_text
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np 
import pandas as pd
import scipy
from collections import OrderedDict
import math
import gseapy as gp


kegg = {} 
ChEA1 = {} #prelim dict used to consolidate line name
ChEA = {} #tf: [targets]
wikipathway = {} #pathway: [genes]
mycpaths = [] #list of pathways containing myc
usf1paths = [] #list of pathways containing usf1
myctargets = [] #list of targets of myc
usf1targets = [] #list of targets of usf1
M0MYC = {} #gene: [logFC, p]
M0MYCup = [] #list of upregulated genes in M0MYC
M0MYCdown = [] #list of downregulated genes in M0MYC
M0USF = {} #gene: [logFC, p]
M0USFup = [] #upregulated
M0USFdown = [] #downregulated
M1USF = {} #gene: [logFC, p]
M1USFup = [] #upregulated
M1USFdown = [] #downregulated
downaged = []
downathero = []
downruptured = []
downstable = []
upaged = []
upathero =[]
upruptured = []
upstable = []
goidict = {}
ageddict = {}
testingup = []
testingdown = []
atherodict = {}
testatheroup = []
testatherodown =[]
deg_list = []
allgenes = {}


goi = ['MMP8', 'FAXC', 'ADH1B', 'HSBP7', 'OGN', 'GAD1', 'SIRT4', 'RGS5', 'FNDC1', 'SDK1', 'ACTA2-AS1', 'GRIA1', 'SCN3A', 'PNMA6A', 'SPEG', 'ITGB1BP2', 'CLEC18C', 'SPTB', 'MYO5B', 'HMNC1', 'TSPAN11', 'SLC28A3', 'AMOTL1', 'LIMCH1', 'SFRP4', 'KIAA1217', 'GYXLT2', 'SGIP1', 'HMCN1', 'ANK3', 'JAM3', 'LAMA2', 'LRP6', 'WDFY3-AS2', 'H2BC15', 'ADGRG1']
#genes of interest (top 30%)
goistats = {} #gene:[logFC, padj]

def openfile(filename, dictname):
  with open(filename) as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
      pathway = line[0]
      genes = []
      for x in line[1:]:
        if x != '':
          genes.append(x)
      dictname[pathway] = genes
#opens the file and seperates out based on tabs
#creates a dictionary of genes in each pathway in the file
def opencsv(filename, dictname):
  with open(filename) as f:
    reader=csv.reader(f,delimiter=',')
    for line in reader:
      if line[0] != '':
        if line[6] != 'NA':
          if float(line[6]) < 0.05:
            dictname[line[9]] = float(line[3]), float(line[6])

def opentsv(filename, dictname):
  with open(filename) as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
      if line[0] != '' and line[0] != 'GeneID':
        if line[1] != 'NA':
          if float(line[1]) < 0.05:
            dictname[line[7]] = float(line[5]), float(line[1])
        #gene = [logFC, p]
        #uses padj to determine if gene is significant
#opens the file and seperates out based on commas, making a dictionary of the data in each line

def opentsv2(filename, dictname):
  with open(filename) as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
      if line[0] != '' and line[0] != 'ID':
        if line[1] != 'NA':
          if float(line[1]) < 0.05:
            dictname[line[6]] = float(line[5]), float(line[1])

def opentsv3(filename, dictname):
  with open(filename) as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
      dictname[line[0]] = float(line[1]), float(line[2])

def openrtf(filename, listname):
  with open(filename) as f:
    for line in f:
      stripped = rtf_to_text(line)
      if stripped.startswith('{') == False and stripped != '' and stripped != '\n':
        listname.append(stripped)
opentsv3('ALLM0genes.csv', allgenes)    
openrtf('DOWN aged.rtf', downaged)
openrtf('UP aged.rtf', upaged)
openrtf('DOWN athero.rtf', downathero)
openrtf('UP athero.rtf', upathero)
openrtf('DOWN ruptured.rtf', downruptured)
openrtf('UP ruptured.rtf', upruptured)
openrtf('DOWN stable.rtf', downstable)
openrtf( 'UP stable.rtf', upstable)
openfile('KEGG_2021_Human.txt', kegg)
openfile('ChEA_2022.txt', ChEA1)
openfile('WikiPathway_2023_Human.txt', wikipathway)
opencsv('M0MYC_allgenes.csv', M0MYC)
opencsv('M0USF1_allgenes.csv', M0USF)
opencsv('M1USF1_allgenes.csv', M1USF)
opentsv('aged.top.table.tsv', ageddict)
opentsv2('athero.top.table.tsv', atherodict)


sorted_items = sorted(ageddict.items(), key=lambda x: abs(x[1][0]), reverse=True)
sortedageddict = OrderedDict(sorted_items)
top250aged = OrderedDict(list(sortedageddict.items())[:250])

def findtfs(dict):
  for x in dict:
    for y in dict[x]:
      if y == 'MYC':
        mycpaths.append(x)
      if y == 'USF1':
        usf1paths.append(x)
      #may be duplicate pathways as adding multiple files to one list
      #could add another parameter to make different lists for different files if i want, but this bit of code isnt really necessary or referenced back to

for x in ChEA1:
  linename = x.split(' ')
  tf = linename[0]
  ChEA[tf] = ChEA1[x]
#sets the dictionary up to be tf : [targets]

for x in ChEA:
  if x == 'MYC':
    for y in ChEA[x]:
      myctargets.append(y)
  elif x == 'USF1' or x == 'FCHL' or x == 'UEF'  or x == 'MLTF'  or x == 'MLTFI' or x == 'HYPLIP1' or x == 'BHLHB11'  or x == 'FCHL1' or x == 'USF':
    for y in ChEA[x]:
      usf1targets.append(y)

findtfs(kegg)
findtfs(wikipathway)

def indict(dict, glist):
  for x in glist:
    for y in dict:
      if x in dict[y]:
        print(x, y)

'''for x in goi:
  if x in myctargets:
    print(x)'''

aaa = [M0MYC, M0USF, M1USF] #to iterate through all datasets

'''for x in goi:
  for y in aaa:
    if x in y:
      goistats[x] = y[x][0], y[x][1]
      #gene: logFC, padj'''

def DE(dataset, DElistup, DElistdown):
  for x in dataset:
    if dataset[x][0] == 'NA':
      continue
    elif float(dataset[x][0]) > 0:
      DElistup.append(x)
    elif float(dataset[x][0]) < 0:
      DElistdown.append(x)
    else:
      continue
      #finds differentially expressed genes
#take out filter to find more still sig hits 

DE(M0MYC, M0MYCup, M0MYCdown)
DE(M0USF, M0USFup, M0USFdown)
DE(M1USF, M1USFup, M1USFdown)
DE(ageddict, testingup, testingdown)
DE(atherodict, testatheroup, testatherodown)

def intersection2(dataset1, dataset2):
  a = []
  for x in dataset1:
    if x in dataset2:
      a.append(x)
  print(len(a))
  print(a)
  #performs an intersection of two lists

def intersection3(d1, d2, d3):
  a = []
  for x in d1:
    if x in d2 and x in d3:
      a.append(x)
  print(len(a))
  print(a)
  #performs an intersection of three lists

#intersection2(downstable, M0USFup)

# creating a DataFrame from the dictionaries
'''
df = pd.DataFrame({
  'M0MYC': [M0MYC.get(gene, [None, None])[0] for gene in goi],
  'M0USF': [M0USF.get(gene, [None, None])[0] for gene in goi],
  'M1USF': [M1USF.get(gene, [None, None])[0] for gene in goi],
}, index=goi)
df.fillna(df.mean(), inplace=True)
df.replace([np.inf, -np.inf], np.nan, inplace=True)
df.fillna(1e6, inplace=True)
# Create a clustered heatmap
clustered_heatmap = sns.clustermap(df, cmap="coolwarm", method="average", col_cluster=True)

# Add a title to the plot
plt.title('Clustered Heatmap of Gene Expression (logFC)')

# Show the plot
plt.show()

# Plot the heatmap
def makeheat(df):
  sns.heatmap(df, annot=True, cmap="coolwarm", cbar_kws={'label': 'logFC'})
  plt.title('Heatmap of Gene Expression (logFC)')
  plt.show()
  '''

'''
for x in usf1paths:
  if x == 'FAM26F' or x == 'CAHLM6':
    print(x)
'''

'''
print("down in knockouts")
intersection2(ChEA['MYC'], M0MYCdown)
print(' ')
print("down in aged")
intersection2(ChEA['MYC'], testingdown)
print(' ')
print("down in aged and knockout")
intersection3(ChEA['MYC'], testingdown, M0MYCdown)
print(' ')
print("up in knockout") 
intersection2(ChEA['MYC'], M0MYCup)
print(' ')
print("up in aged")
intersection2(ChEA['MYC'], testingup)
print(' ')
print("up in knockout and aged")
intersection3(ChEA['MYC'], testingup, M0MYCup)

print('-----------------')
print("This is me checking i havent mixed up and down")
print("up in aged down in knockouts")
intersection3(ChEA['MYC'], testingup, M0MYCdown)
print("down in aged up in knockouts")
intersection3(ChEA['MYC'], testingdown, M0MYCup)

'''
'''
print('-----------------')
print("athero genes")
print(' ')
print("down in athero")
intersection2(ChEA['MYC'], downruptured)
print(' ')
print("up in athero")
intersection2(ChEA['MYC'], upruptured)
print(' ')
print("down in athero and down in knockouts")
intersection3(ChEA['MYC'], downruptured, M0MYCdown)
print(' ')
print("up in athero and down in knockouts")
intersection3(ChEA['MYC'], upruptured, M0MYCdown)
print(' ')
print("down in athero and up in knockouts")
intersection3(ChEA['MYC'], downruptured, M0MYCup)
print(' ')
print("up in athero and up in knockouts")
intersection3(ChEA['MYC'], upruptured, M0MYCup)

print('-----------------')
print("rupture genes")
print(' ')
print("down in rupture")
intersection2(ChEA['MYC'], downruptured)
print(' ')
print("up in rupture")
intersection2(ChEA['MYC'], upruptured)
print(' ')
print("down in rupture and down in knockouts")
intersection3(ChEA['MYC'], downruptured, M0MYCdown)
print(' ')
print("up in rupture and down in knockouts")
intersection3(ChEA['MYC'], upruptured, M0MYCdown)
print(' ')
print("down in rupture and up in knockouts")
intersection3(ChEA['MYC'], downruptured, M0MYCup)
print(' ')
print("up in rupture and up in knockouts")
intersection3(ChEA['MYC'], upruptured, M0MYCup)

print('-----------------')
print("stable genes")
print(' ')
print("down in stable")
intersection2(ChEA['MYC'], downstable)
print(' ')
print("up in stable")
intersection2(ChEA['MYC'], upstable)
print(' ')
print("down in stable and down in knockouts")
intersection3(ChEA['MYC'], downstable, M0MYCdown)
print(' ')
print("up in stable and down in knockouts")
intersection3(ChEA['MYC'], upstable, M0MYCdown)
print(' ')
print("down in stable and up in knockouts")
intersection3(ChEA['MYC'], downstable, M0MYCup)
print(' ')
print("up in stable and up in knockouts")
intersection3(ChEA['MYC'], upstable, M0MYCup) 

'''

'''
print('LogFC and p of genes up/downreg')
for x in ChEA['MYC']:
  gene = x
  alreadydone = []
  if x in testingdown and x in M0MYCdown:
    print("------------------------")
    print("*** " + x + " ***")
    alreadydone.append(x)
    print("Down in aged and knockout")
    print("Knockout stats: ")
    print( M0MYC[x])
    print("pre mono stats: ")
    print(ageddict[x])
    print("LogFC Difference:")
    diff = M0MYC[x][0] - ageddict[x][0]
    if diff > 0:
      print("     Knockout logFC bigger by " + str(diff))
    else:
      print("     pre-mono logFC bigger by " + str(diff*-1))
    print("Present in wikipathways: ")
    for key, valuelist in wikipathway.items(): 
      for item in valuelist:
        if item == gene:
          print("     " + key)
    print("Present in kegg pathways: ")
    for key, valuelist in kegg.items(): 
      for item in valuelist:
        if item == gene:
          print("     " + key)

  if x in testingup and x in M0MYCdown:
    print("------------------------")
    print("*** " + x + " ***")
    alreadydone.append(x)
    print("Up in aged and down in knockout")
    print("Knockout stats: ")
    print( M0MYC[x])
    print("pre mono stats: ")
    print(ageddict[x])
    print("LogFC Difference:")
    diff = M0MYC[x][0] - ageddict[x][0]
    if diff > 0:
      print("     Knockout logFC bigger by " + str(diff))
    else:
      print("     pre-mono logFC bigger by " + str(diff*-1))
    print("Present in wikipathways: ")
    for key, valuelist in wikipathway.items(): 
      for item in valuelist:
        if item == gene:
          print("     " + key)
    print("Present in kegg pathways: ")
    for key, valuelist in kegg.items(): 
      for item in valuelist:
        if item == gene:
          print("     " + key)


  if x in testingdown and x in M0MYCup:
    print("------------------------")
    print("*** " + x + " ***")
    alreadydone.append(x)
    print("Down in aged and up in knockout")
    print("Knockout stats: ")
    print( M0MYC[x])
    print("pre mono stats: ")
    print(ageddict[x])
    print("LogFC Difference:")
    diff = M0MYC[x][0] - ageddict[x][0]
    if diff > 0:
      print("     Knockout logFC bigger by " + str(diff))
    else:
      print("     pre-mono logFC bigger by " + str(diff*-1))
    print("Present in wikipathways: ")
    for key, valuelist in wikipathway.items(): 
      for item in valuelist:
        if item == gene:
          print("     " + key)
    print("Present in kegg pathways: ")
    for key, valuelist in kegg.items(): 
      for item in valuelist:
        if item == gene:
          print("     " + key)

      
  if x in testingup and x in M0MYCup:
    print("------------------------")
    print("*** " + x + " ***")
    alreadydone.append(x)
    print("Up in aged and in knockout")
    print("knockout stats: ")
    print(M0MYC[x])
    print("pre mono stats: ")
    print(ageddict[x])
    print("LogFC Difference:")
    diff = M0MYC[x][0] - ageddict[x][0]
    if diff > 0:
      print("     Knockout logFC bigger by " + str(diff))
    else:
      print("     pre-mono logFC bigger by " + str(diff*-1))
      print("Present in wikipathways: ")
      for key, valuelist in wikipathway.items(): 
        for item in valuelist:
          if item == gene:
            print("     " + key)
      print("Present in kegg pathways: ")
    for key, valuelist in kegg.items(): 
      for item in valuelist:
        if item == gene:
          print("     " + key)
'''
if x in testatherodown and x in M0MYCdown:
    print("------------------------")
    if x not in alreadydone:
      print("*** " + x + " ***")
      alreadydone.append(x)
    print("Down in athero and knockout")
    print("Knockout stats: ")
    print(M0MYC[x])
    print("athero stats: ")
    print(atherodict[x])
    print("LogFC Difference:")
    diff = M0MYC[x][0] - atherodict[x][0]
    if diff > 0:
      print("     Knockout logFC bigger by " + str(diff))
    else:
      print("     athero logFC bigger by " + str(diff*-1))
    print("Present in wikipathways: ")
    for key, valuelist in wikipathway.items(): 
      for item in valuelist:
        if item == gene:
          print("     " + key)
'''


  if x in testatheroup and x in M0MYCup:
    print("------------------------")
    if x not in alreadydone:
      print("*** " + x + " ***")
      alreadydone.append(x)
    print("Up in athero and in knockout")
    print("Knockout stats: ")
    print(M0MYC[x])
    print("athero stats: ")
    print(atherodict[x])
    print("LogFC Difference:")
    diff = M0MYC[x][0] - atherodict[x][0]
    if diff > 0:
      print("     Knockout logFC bigger by " + str(diff))
    else:
      print("     athero logFC bigger by " + str(diff*-1))
    print("Present in wikipathways: ")
    for key, valuelist in wikipathway.items(): 
      for item in valuelist:
        if item == gene:
          print("     " + key)



print("-----------------------------------------")
print(" ")
for x in finalcandiates:
  print(" ")
  print(x) 
  if x in M0MYC:
    print("knockout stats:")
    print(M0MYC[x])
  else:
    print("not in knockout")
  if x in atherodict:
    print("athero stats: ")
    print(atherodict[x])
  else:
    print("not in athero")'''

    
datasetdescriptions = {}
datasetdescriptions['aged'] = "mRNA profiles of human macrophages (CD45+, CD163+, CD11b+, Side-Scatter-Low) from central marrow (CM) and spicule associated cells (SAC) isolated from the marrow aspirate of healthy human donors by Illumina HiSeq 2500 sequencing."
datasetdescriptions['athero'] = "Plaques were obtained from 4 patients agand dissected in stable and unstable regions. RNAs were extracted from the two sections and sent for RNAseq (Genewiz)."
datasetdescriptions['knockout'] = "Charlottes M0 MYC knockdown ageing model"

finalcandidates = ['FCHO2', 'H1-2', 'EMILIN1', 'EPB41L2', 'ATP6V1H', 'RNF44', 'ZNF589', 'YWHAG', 'UTP3', 'SLCO5A1', 'TMED9', 'DAPK1', 'RRAS']
genepvalues = []

for x in finalcandidates:
  genepvalues.append(M0MYC[x][1])

def calculate_fdr(p_values):

#Calculate False Discovery Rate (FDR) using Benjamini-Hochberg procedure.
#Args:
#- p_values: List or array of p-values
#Returns:
#- fdr: List of False Discovery Rates corresponding to input p-values
  
  m = len(p_values)  # Total number of tests
  sorted_indices = sorted(range(m), key=lambda i: p_values[i])  # Indices sorted by p-value
  sorted_p_values = [p_values[i] for i in sorted_indices]  # Sorted p-values

  fdr = [0] * m
  for i in range(m):
    fdr[i] = sorted_p_values[i] * m / (i + 1)

  # Adjust for monotonicity
  for i in range(m - 2, -1, -1):
    fdr[i] = min(fdr[i], fdr[i + 1])

# Back to original order
  fdr_by_index = [0] * m
  for i, index in enumerate(sorted_indices):
    fdr_by_index[index] = fdr[i]

  print(fdr_by_index)




def opencsvNOP(filename, dictname):
  with open(filename) as f:
    reader=csv.reader(f,delimiter=',')
    for line in reader:
      if line[0] != '':
        if line[6] != 'NA':
          dictname[line[9]] = float(line[3]), float(line[6])


'''    
M0MYCNOP = {}
opencsvNOP('M0MYC_allgenes.csv', M0MYCNOP)
    
MYCtargetsanalysis = (ChEA['MYC'])
deletelist = []
for x in M0MYCNOP.keys():
  if x not in MYCtargetsanalysis:
    deletelist.append(x)

for y in deletelist:
  del(M0MYCNOP[y])


neglog10dict = {}
for key, value in M0MYCNOP.items():
  neglog10dict[key] = [value[0], -math.log10(value[1])]
  
for x in neglog10dict:
  #print(x)
  print(str(neglog10dict[x][1]))'''
    
  #str(neglog10dict[x][0]))
        
   #     + "\t" + str(neglog10dict[x][1]))


'''
print("Pathways enriched in list of upregulated genes in aged")
enrichment_results = gp.enrichr(gene_list=upaged, gene_sets='KEGG_2019_Human')

# Extract and print the names of the top 5 enriched pathways
top_pathways = enrichment_results.res2d.head(10)  # Select the top 5 enriched pathways
for index, row in top_pathways.iterrows():
    print(f"{row['Term']}: p-value = {row['P-value']}")

print("Pathways enriched in list of downregulated genes in aged")
enrichment_results = gp.enrichr(gene_list=downaged, gene_sets='KEGG_2019_Human')

# Extract and print the names of the top 5 enriched pathways
top_pathways = enrichment_results.res2d.head(10)  # Select the top 5 enriched pathways
for index, row in top_pathways.iterrows():
    print(f"{row['Term']}: p-value = {row['P-value']}")

print("Pathways enriched in MYC targets ")
enrichment_results = gp.enrichr(gene_list=ChEA['MYC'], gene_sets='KEGG_2019_Human')

# Extract and print the names of the top 5 enriched pathways
top_pathways = enrichment_results.res2d.head(10)  # Select the top 5 enriched pathways
for index, row in top_pathways.iterrows():
    print(f"{row['Term']}: p-value = {row['P-value']}")
'''

#for x in M0MYC:
  #if x in ChEA['MYC']:
    #print(x, M0MYC[x][0], M0MYC[x][1])

targetslist = []

for x in allgenes:
  if x in ChEA['MYC']:
    targetslist.append(x)
    #print(x, allgenes[x])

#print(len(allgenes))
#print(len(targetslist))
#print(len(myctargets))

'''counter = 0
for x in ChEA:
  counter += 1
  if x == 'MYC':
    print(counter)
    break'''

for x in ChEA1:
  linename = x.split(' ')
  tf = linename[0]
  if tf == 'MYC':
    print(linename)


