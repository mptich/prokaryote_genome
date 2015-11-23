# Global definitions and utilities for the scripts

# Directory containing subdirecories with organisms' genomes
def PROKARYOTS_DIR():
    return ROOT_DIR() + "BACTERIA/"

# Directory for the files created by utilities
def WORK_FILES_DIR():
    return ROOT_DIR() + "WorkFiles/"

# Input on taxonomy
def TAXONOMY_FILE():
    return ROOT_DIR() + "Taxonomy/prokaryot_families.csv"

# List of all directories with Prokaryot genomes
def PROKARYOT_DIRS_FILE():
    return ROOT_DIR() + "BACTERIA/genome_dirs.txt"

# Dump of the dictionary directory -> ProkGenome
def PROK_GENOME_DICT():
    return WORK_FILES_DIR() + "genome_dict.json"

# Returns the absolute path to the root directory
def ROOT_DIR():
    return "/Users/morel/biology/"

# Shared utilities directory
def SHARED_PROG_DIR():
    return ROOT_DIR() + "devel/shared/"




