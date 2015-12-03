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

# Dump of ProkCogs
def PROK_COGS_SET():
    return config.WORK_FILES_DIR() + "cogs_set.json"

# Dump of ProkDna.key() -> Taxa
def PROK_TAXA_DICT():
    return config.WORK_FILES_DIR() + "taxa_dict.json"

# Set of organism names from the Taxonomy file that could not be matched
# against the names in teh Prokaryotic database
def UNMATCHED_TAXA_SET():
    return config.WORK_FILES_DIR() + "unmatched_taxa_set.json"


