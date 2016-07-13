# Definitions of the relative file names

import config

# CSV file on Taxonomy classification
def TAXONOMY_FILE():
    return config.WORK_FILES_DIR() + "taxonomy.csv"

# List of all directories with Prokaryot genomes
def PROKARYOT_DIRS_FILE():
    return config.PROKARYOTS_DIR() + "genome_dirs.txt"

# Dump of the dictionary directory -> ProkGenome
def PROK_GENOME_DICT():
    return config.WORK_FILES_DIR() + "genome_dict.json"

# Dump ProkDna key -> ProkDna dictionary
def PROK_DNA_DICT():
    return config.WORK_FILES_DIR() + "dna_dict.json"

# Dump of the dictionary directory -> ProkDnaSet
def PROK_CLEAN_GENOME_DICT():
    return config.WORK_FILES_DIR() + "chromosome_dict.json"

# Dump of CogInst list
def COG_INST_LIST():
    return config.WORK_FILES_DIR() + "cog_inst_list.json"

# Genome dir -> set of COG names
def COG_DICT():
    return config.WORK_FILES_DIR() + "cog_dict.json"

# Dump of sample CogInst list
def SAMPLE_COG_INST_LIST():
    return config.WORK_FILES_DIR() + "sample_cog_inst_list.json"

# dir1, dir2 -> COG distance
def COG_DIST_DICT():
    return config.WORK_FILES_DIR() + "cog_dist_dict.json"

# Dump of Cog list
def COG_LIST():
    return config.WORK_FILES_DIR() + "cog_list.json"

# Dump of list of (genome, COG count) tuples
def GENOME_COG_CNT_LIST():
    return config.WORK_FILES_DIR() + "genome_cog_cnt_list.json"

# Dump of ProkDnaSet.key -> Taxa
def PROK_TAXA_DICT():
    return config.WORK_FILES_DIR() + "taxa_dict.json"

# Set of organism names from the Taxonomy file that could not be matched
# against the names in teh Prokaryotic database
def UNMATCHED_TAXA_SET():
    return config.WORK_FILES_DIR() + "unmatched_taxa_set.json"

# Set of cleaned ProkDna directories that have not matched anything in
# taxonomy
def UNMATCHED_PROC_DNA_SET():
    return config.WORK_FILES_DIR() + "unmatched_proc_dna_set.json"

# Dict of genome dirs -> list of counts for their taxonomy distances
def GENOME_TAX_DIST_CNT_DICT():
    return config.WORK_FILES_DIR() + "genome_tax_dist_cnt_dict.json"

# Dict of TaxType key -> list of counts for their taxonomy distances
def TAXTYPE_TAX_DIST_CNT_DICT():
    return config.WORK_FILES_DIR() + "taxtype_tax_dist_cnt_dict.json"

def COG_WEIGHTS_DICT_LIST():
    return config.WORK_FILES_DIR() + "cog_weights_dict_list.json"

# List of UtilObject's for reclassified genomes with old neighbors
# in the taxonomy classification
def RECLASSIFIED_DIR_LIST():
    return config.WORK_FILES_DIR() + "reclassified_dir_list.json"
