#!/usr/bin/env python
"""
This script will take:

    * combined table of gene counts and annotations 
    * clade definitions YAML

it will produce:

    * table of annotated counts for each  requested clade

Output table files will be named with the input file plus the clade name
unless a naming prefix is given iwth -o <prefix>

Usage:
  split_table_mult.py [options] <annot_counts> <tax_yaml>
  split_table_mult.py -h | --help
  split_table_mult.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.
  -o <out_base> Name output files with this prefix

"""
import docopt
import yaml
import pandas

def main():
    """
    The body of this script

    The docopt arguments will look like this
    {
      "--help": false, 
      "--version": false, 
      "<annot_counts>": "table.tsv", 
      "<tax_yaml>": "clades.yml"
      "-o": None
    }
    """

    # parse the command line options
    arguments = docopt.docopt(__doc__, version='0.0.1')
    hit_counts_annot_file = arguments['<annot_counts>']
    taxon_definitions_file = arguments['<tax_yaml>']
    output_prefix = arguments['-o'] if arguments['-o'] \
                                            else hit_counts_annot_file

    # load the BIG table
    hit_counts_annot = pandas.read_csv(hit_counts_annot_file,
                                       sep='\t',
                                       index_col=0,
                                       header=0)

    # loop over clade definitions
    for clade_name, keep_taxa, drop_taxa in get_taxon_definitions(taxon_definitions_file):

        # filter raw data table on taxonomic designations
        hit_counts_filt_annot = filter_by_taxon(hit_counts_annot,
                                                keep_taxa,
                                                drop_taxa)

        # write out filtered file
        hit_counts_filt_annot_file = f"{output_prefix}.{clade_name}"
        hit_counts_filt_annot.to_csv(hit_counts_filt_annot_file, sep='\t')


def get_taxon_definitions(taxon_definitions_file):
    """
    read in a YAML file defining clades. the structure should look like:

    Crocosphaera:
        keep:
            genus:
                - Crocoshaera
    Prochlorococcus:
        keep:
            family:
                - Cyanobiaceae
        drop:
            genus:
                - Synechococcus_C
                - Vulcanococcus
                - Unknown
                - Synechococcus_D
                - PCC7001
                - RCC307
                - BAIKAL-G1
                - Cyanobium
                - Synechococcus_B
                - WH-5701
    """
    with open(taxon_definitions_file) as yaml_in:
        tax_defs = yaml.load(yaml_in, Loader=yaml.SafeLoader)

    for clade_name, clade_defs in tax_defs.items():
        yield clade_name, clade_defs['keep'], clade_defs.get('drop',{})


def filter_by_taxon(hit_counts_annot, taxa_to_keep, taxa_to_drop):
    """
    return only rows that match taxa_To_keep and not_taxa_to_Drop

    do this using DAtaFrame.query()

    eg: 
        taxa_to_keep = {'genus': ['Crocosphaera']}
    becomes:
        hit_counts_annot.query('genus == "Crocosphaera"')
    """
    # make sure single taxa are in lists
    taxa_to_keep = _check_tax_def_format(taxa_to_keep)
    taxa_to_drop = _check_tax_def_format(taxa_to_drop)

    keep_query_string = " and ".join(f'{rank} == "{name}"'
                                     for rank in taxa_to_keep
                                     for name in taxa_to_keep[rank]
                                    )
    drop_query_string = " and ".join(f'{rank} != "{name}"'
                                     for rank in taxa_to_drop
                                     for name in taxa_to_drop[rank]
                                    )

    # let's assume the keep list is nonzero, but check the drop list
    if len(drop_query_string) > 0:
        query_string = keep_query_string + " and " + drop_query_string
    else:
        query_string = keep_query_string

    # run the query and return the filtered result
    return hit_counts_annot.query(query_string)

def _check_tax_def_format(taxa_def_dict):
    return {rank: ([names,] if isinstance(names, str) else names) 
            for rank, names in taxa_def_dict.items()}

if __name__ == "__main__":
    main()
