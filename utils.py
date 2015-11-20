# Global definitions and utilities for the scripts

# Directory containing subdirecories with organisms' genomes
def ORGANISMS():
    return "../BACTERIA"

# Input on taxonomy
def TAXONOMY():
    return "../Taxonomy/prokaryot_families.csv"

# List of all directories with genomes
def GENOME_DIRS_FILE():
    return "./genome_dirs.txt"

# Dump of the dictionary directory -> ProkGenome
def GENOME_DICT_JSON():
    return "./genome_dict.json"

# Make our exception more specific
class GenomeError(Exception):

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
