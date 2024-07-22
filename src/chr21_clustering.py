## This is a program that explores the functional categorization of chr21 genes. 

## Libraries
import gseapy
import re
import matplotlib.pyplot as plt


## read in the chromosome annotation file form migsigdb
chr_gmt = gseapy.read_gmt('/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/data/migsigdb/c1.all.v2023.2.Hs.symbols.gmt')

## find chromosome 21 genes
## use regex for 'chr21' key 
pattern = re.compile('chr21')

# Search for keys that match the pattern
chr21_dict = {key:value for key, value in chr_gmt.items() if pattern.match(key)}
chr21_gene_list = list(chr21_dict.values())
chr21_gene_list = [item for sublist in chr21_gene_list for item in sublist]

## read in the hallmark and curated pathway gene sets from migsigdb
hallmark_gmt = gseapy.read_gmt('/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/data/migsigdb/h.all.v2023.2.Hs.symbols.gmt')

curated_gmt = gseapy.read_gmt('/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/data/migsigdb/c2.all.v2023.2.Hs.symbols.gmt')

## find the overlaps between the list of chr21 genes and each pathway

def find_overlaps(target_list, source_dict):
    # Initialize a dictionary to store the overlaps
    overlaps = {}
    
    # Convert the target list to a set for efficient comparison
    target_set = set(target_list)
    
    # Iterate over the dictionary items
    for key, values in source_dict.items():
        # Find the intersection between the list and the current dictionary values
        overlap = target_set.intersection(values)
        # Store the overlap in the dictionary
        overlaps[key] = list(overlap)
    
    return overlaps

# Find the overlaps

## Hallmark
hallmark_overlaps = find_overlaps(chr21_gene_list, hallmark_gmt)

    # count the overlaps
hallmark_overlaps_count = {}

# Calculate the length of each list value in the dictionary
for key, value in hallmark_overlaps.items():
    hallmark_overlaps_count[key] = len(value)

hallmark_overlaps_count = {k: v for k, v in sorted(hallmark_overlaps_count.items(), key=lambda item: item[1], reverse=True)}

## reoder hallmark overlaps by key count keys
hallmark_overlaps = {key: hallmark_overlaps[key] for key in hallmark_overlaps_count}



# Filter the dictionary to include only pairs with values greater than the threshold
hallmark_overlaps_top = dict(list(hallmark_overlaps_count.items())[:10])

# Create a bar plot
plt.clf()
plt.barh(range(len(hallmark_overlaps_top)), list(hallmark_overlaps_top.values()), align='center')
plt.yticks(range(len(hallmark_overlaps_top)), list(hallmark_overlaps_top.keys()), fontsize = 8)
plt.xlabel('Overlap')
plt.title('Chr21 Overlap of Hallmark Gene Sets')
plt.show()

for key, overlaps in hallmark_overlaps.items():
    print(f"{key}: {overlaps}")



## Curated
curated_overlaps = find_overlaps(chr21_gene_list, curated_gmt)

    # count the overlaps
curated_overlaps_count = {}

# Calculate the length of each list value in the dictionary
for key, value in curated_overlaps.items():
    curated_overlaps_count[key] = len(value)

curated_overlaps_count = {k: v for k, v in sorted(curated_overlaps_count.items(), key=lambda item: item[1], reverse=True)}

## reoder curated overlaps by key count keys
curated_overlaps = {key: curated_overlaps[key] for key in curated_overlaps_count}



# Filter the dictionary to include only pairs with values greater than the threshold
curated_overlaps_top = dict(list(curated_overlaps_count.items())[:10])

# Create a bar plot
plt.clf()
plt.barh(range(len(curated_overlaps_top)), list(curated_overlaps_top.values()), align='center')
plt.yticks(range(len(curated_overlaps_top)), list(curated_overlaps_top.keys()), fontsize = 8)
plt.xlabel('Overlap')
plt.title('Chr21 Overlap of Curated Gene Sets')
plt.show()

for key, overlaps in curated_overlaps.items():
    print(f"{key}: {overlaps}")
