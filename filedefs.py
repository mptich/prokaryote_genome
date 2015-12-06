# Definitions of the relative file names

import config

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

# Dump of CogInst set
def COG_INST_SET():
    return config.WORK_FILES_DIR() + "cog_inst_set.json"

# Dump of Cog list
def COG_LIST():
    return config.WORK_FILES_DIR() + "cog_list.json"

# Dump of list of (genome, COG count) tuples
def GENOME_COG_CNT_LIST():
    return config.WORK_FILES_DIR() + "genome_cog_cnt_list.json"

# Dump of ProkDnaSet.key() -> Taxa
def PROK_TAXA_DICT():
    return config.WORK_FILES_DIR() + "taxa_dict.json"

# Set of organism names from the Taxonomy file that could not be matched
# against the names in teh Prokaryotic database
def UNMATCHED_TAXA_SET():
    return config.WORK_FILES_DIR() + "unmatched_taxa_set.json"


