#!/usr/bin/env python
"""
Usage:
  plot_counts.py <tax_name> <count_file> <pdf_file>
  plot_counts.py -h | --help
  plot_counts.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.

"""
import docopt
import re

import pandas
import numpy as np
import sklearn
import xlrd

import plotnine as p9
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from scipy.stats import kruskal
from statsmodels.sandbox.stats.multicomp import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from adjustText import adjust_text

# constants
CruiseSamp_metaD_file = '/data/urisheyn/Mesoscope_CruiseData/170801_MESOSCOPE_DigitizedLogbookFinal20190508.txt'

# map group names to regular expression patterns for finding samples
sample_group_patterns = {
    'SURF_L1': r'SL1.+015$',
    'SURF_L2': r'SL2.+015$',
    'DCM_L1': r'SL1.+\d\d(?<!015)$',
    'DCM_L2': 'SL2.+\d\d(?<!015)$',
}

#comparisons and the interesting ones:
# lists of sample for each comparison
SURFsamples = ['SURF_L2', 'SURF_L1']
DCMsamples  = ['DCM_L2', 'DCM_L1']
L1samples   = ['SURF_L1', 'DCM_L1']
L2samples   = ['SURF_L2','DCM_L2']


# Dictionary for comparisons samples
compar = {'SURF_L2vsL1': SURFsamples,
          'DCM_L2vsL1':DCMsamples,
          'L1_SURFvsDCM':L1samples,
          'L2_SURFvsDCM':L2samples
         }

compar_keys = ['SURF_L2vsL1', 'DCM_L2vsL1', 'L1_SURFvsDCM', 'L2_SURFvsDCM']
import_comp = ['SURF_L2vsL1', 'DCM_L2vsL1']

ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum']
other_annots = ['function', 'COG', 'EGGNOG', 'KEGG', 'PFAM', 'TIGRFAM']
annot_for_analysis = ranks + other_annots


def main():

    # parse the command line options
    arguments = docopt.docopt(__doc__, version='0.0.1')
    # plot_counts.py <tax_name> <count_file> <pdf_file>
    pdf_plots(arguments['<tax_name>'],
              arguments['<count_file>'],
              arguments['<pdf_file>'])

def pdf_plots(tax_name, tax_counts_file, pdf_file, rank='genus'):

    # loading the filt hit count data
    hit_counts_filt_annot = pandas.read_csv(tax_counts_file, 
                                            sep='\t', 
                                            index_col=0, 
                                            header=0)

    # separate data from annots
    #data columns:
    data_columns = [c for c in hit_counts_filt_annot.columns if c.startswith("MS")]
    #non data columns:
    non_data_columns = [c for c in hit_counts_filt_annot.columns if not c.startswith("MS")]

    # use re to build sample lists
    time_courses = {
        name:[c for c in data_columns if re.search(pattern, c)]
        for name, pattern in sample_group_patterns.items()
    }

    #substituting Nan with unkown 
    hit_counts_filt_annot.fillna('Unknown', inplace=True)

    # making list of unique value for each annotation in a dic 
    #  and giving it a number to map for aq unique color for down stream ploting
    annot_color_indexes = {}
    for annot in annot_for_analysis:
        unique_annots = np.unique(hit_counts_filt_annot[annot])
        annot_color_indexes[annot] = pandas.Series(
            {label:i for i, label in enumerate(unique_annots)},
            name=f"{annot}_color_#"
        )

    #putting the unique number per annot columns into the annot dataframe under a new copy.
    annotations_color = hit_counts_filt_annot[annot_for_analysis]

    # merging the annot dataframe with the unique number per value dataframe 
    for annot in annot_for_analysis:
        # on the annot column for each iteration
        annotations_color = annotations_color.join(annot_color_indexes[annot],
                                                   on=annot,
                                                   how='left')

    rank_sums = aggregate_on_column(hit_counts_filt_annot,
                                                data_columns,
                                                rank) 

    #DataFrame with sum per genus values in the shape of the expression dataframe
    genusDevValues = hit_counts_filt_annot[[rank,]].join(rank_sums, on=rank)[data_columns]

    # division of every taxon expression value by the sum of its genus in the coresponding sample: 
    #nan are coming from the devision of 0 by sum = 0 only! (genes of genus that are not detected at all for a specific location)
    genus_norm_data = pandas.DataFrame(hit_counts_filt_annot[data_columns]/genusDevValues) 
    genus_norm_data.fillna(0, inplace=True) #rplacing nan by 0 

    #merging the expression data (genus norm and not) with the new colored annotation dataframe 
    Norm_filter_taxa_data_annot = hit_counts_filt_annot[data_columns].join(annotations_color, 
                                                                           how='left')
    #substituting Nan with unkown 
    Norm_filter_taxa_data_annot.fillna('Unknown', inplace=True)
    genus_norm_exp_data_annot = genus_norm_data.join(annotations_color, how='left')
    #substituting Nan with unkown 
    genus_norm_exp_data_annot.fillna('Unknown', inplace=True)


    #sum per annot dic genus norm data (only if more than one unique value)
    SumPerAnnot_genus_norm_data = {annot:aggregate_on_column(genus_norm_exp_data_annot,
                                                                         data_columns,
                                                                         annot) 
                                   for annot in annot_for_analysis
                                   if len(annot_color_indexes[annot]) > 1
                                  } 

    #sum per annot dic norm data
    SumPerAnnot_norm_data = {annot:aggregate_on_column(Norm_filter_taxa_data_annot,
                                                                   data_columns,
                                                                   annot) 
                             for annot in annot_for_analysis
                             if len(annot_color_indexes[annot]) > 1
                            } 

    # open the PDF file for plotting
    pdf = PdfPages(pdf_file)

    # plot the normalized totals grouped y taxon for various ranks
    for annot in ranks:
        if annot not in SumPerAnnot_genus_norm_data:
            # skip anything with only one value
            continue
        annot_table = SumPerAnnot_norm_data[annot]
        f = barplot_per_annot(annot, annot_table, time_courses, data_columns,
                                          figsize=[20,16])
        pdf.savefig(bbox_inches='tight')
        plt.close()

    # plot the normalized totals: one quartet of bar plots per annotation column
    for annot in other_annots:
        annot_table = SumPerAnnot_genus_norm_data[annot]
        f = barplot_per_annot(annot, annot_table, time_courses, data_columns,
                                          figsize=[20,16])
        pdf.savefig(bbox_inches='tight')
        plt.close()

    # line plots of the normalized and genus-normalized data 
    f = line_plots(Norm_filter_taxa_data_annot, genus_norm_exp_data_annot, 
                               time_courses, rank, name=tax_name, figsize=[16,20])
    pdf.savefig(bbox_inches='tight')
    plt.close()

    pdf.close() 

#########
# local functions (consider moving to the background environment)

#Kruskal wraper function for apply
def WrapKruskal (row, splits=[18]):
    sets = np.split(row.values, splits)
    try:
        Pval = kruskal(*sets).pvalue 
    except ValueError as e:
        if e.args == ('All numbers are identical in kruskal',):
            # rows are the same!
            Pval = None
        else:
            raise
    return Pval

#Making a comparison table between the different mesoscope locations (maight need correction)
def make_comparison_tables(Compar, time_courses, norm_data_per_taxa):
    
    matrix = np.array(norm_data_per_taxa[[c for c in norm_data_per_taxa.columns if c.startswith("MS")]])
    ExpressionValues = matrix.reshape(matrix.shape[0] * matrix.shape[1])
    ExpMinVal = ExpressionValues[ExpressionValues.nonzero()].min()
    ExpMaxVal = ExpressionValues[ExpressionValues.nonzero()].max()
    AddedVal = ExpMinVal/ExpMaxVal
    
    Norm_Comparisons = {}
    for c in Compar:
        coursnames = Compar[c]
        cols = time_courses[coursnames[0]] + time_courses[coursnames[1]]
        SampData = norm_data_per_taxa[cols]
        set_1 =  norm_data_per_taxa[time_courses[coursnames[0]]]
        set_2 =  norm_data_per_taxa[time_courses[coursnames[1]]]
        FoldChange = (set_1.mean(axis = 1)+AddedVal) /(set_2.mean(axis = 1)+AddedVal) 
        new_norm_data_per_taxa=norm_data_per_taxa.copy()
        new_norm_data_per_taxa[coursnames[0]+"_mean"] = set_1.mean(axis = 1)
        new_norm_data_per_taxa[coursnames[1]+"_mean"] = set_2.mean(axis = 1)
        new_norm_data_per_taxa[c+"_FoldChange"] = FoldChange
        new_norm_data_per_taxa[c+"_pval"] = SampData.apply(WrapKruskal, axis= 1)    
        pval_series = new_norm_data_per_taxa[c + "_pval"].dropna()
        pval_index = pval_series.index
        passes, corr_pvalues, sidak, bonf = multipletests(pval_series.values, method='fdr_bh')
        corr_pval_series = pandas.Series(data=corr_pvalues, index=pval_index, name=c+"_corr_pvalues")
        new_norm_data_per_taxa = new_norm_data_per_taxa.join(corr_pval_series)
        Norm_Comparisons[c] = new_norm_data_per_taxa[[c + "_FoldChange",]].join(corr_pval_series, how='right')
    return Norm_Comparisons


#Aggregating expression level in eache sample (column) by list of annotations
def aggregate_on_column(data_table, data_columns, annot_column):
    return data_table[data_columns + [annot_column,]].groupby(annot_column).agg(sum)


def load_meta_data(CruiseSamp_metaD_file):

    #laoding the sample metadata cruise info
    #cruise metadata including the volume filtered for each sample:
    CruiseSamp_metaD = pandas.read_csv(CruiseSamp_metaD_file, sep='\t', index_col=0, header=0)
    CruiseSamp_metaD['Mesoscope_Trans_sample_name'] = CruiseSamp_metaD.index

    # color for each location
    Loc_key_color = pandas.DataFrame(CruiseSamp_metaD.Location.unique())
    Loc_key_color['color'] = ("red", "green", "blue", "black")
    Loc_key_color.columns = ("Location","color")

    CruiseSamp_metaD = CruiseSamp_metaD.merge(Loc_key_color,
                                              on='Location',
                                              left_index=True,
                                              right_index=False)
    CruiseSamp_metaD.set_index('Mesoscope_Trans_sample_name',
                               inplace=True)

    #Changing date time format for the cruise metadata
    def read_date(date):
        return xlrd.xldate.xldate_as_datetime(date, 0) # wrapper function for apply:

    # joining the date and time columns
    CruiseSamp_metaD['Date_Time'] = CruiseSamp_metaD.Date + CruiseSamp_metaD.Time_Start
    # apply on read_date to get new colmun with new format
    CruiseSamp_metaD['newdateTime'] = \
        pandas.to_datetime(CruiseSamp_metaD['Date_Time'].apply(read_date), 
                           errors='coerce')
    # extracting time of day rom the date new time 
    CruiseSamp_metaD['TimeOfDay'] = CruiseSamp_metaD.newdateTime.dt.time 
    CruiseSamp_metaD['day.night'] = np.where(CruiseSamp_metaD['DayTimePoints'] > 4,
                                             'night', 
                                             'day')

    return CruiseSamp_metaD

def barplot_per_annots(SumPerAnnot_norm_data, 
                      time_courses,
                      data_columns):
    """
    Creates multiple figures, one for each annoation column. Best for use in Jupyter to generate multiple figures quickly
    """
    for annot in SumPerAnnot_norm_data:
        annot_expression_data =  SumPerAnnot_norm_data[annot].copy()
        barplot_per_annot(annot, annot_expression_data,
                          time_courses,
                          data_columns)

def barplot_per_annot(annot, annot_expression_data,
                      time_courses, data_columns, figsize=[30,20], ):
    """
    Generates a single barplot for the given annotation. Sutiable for creating a single DF page.
    """

    annot_expression_data['sum'] = annot_expression_data[data_columns].sum(1)

    annot_expression_data.sort_values('sum', axis=0, ascending=False, inplace=True, kind='quicksort', na_position='last')

    top_annot_expression_data = annot_expression_data.iloc[0:30,:][annot_expression_data.iloc[0:30,:].index!='Unknown']

    #top_annot_expression_data.loc['other'] = annot_expression_data.iloc[30:,:].sum() 

    time_points = list(range(0,72,4))

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=figsize, sharex=True, sharey=True)
    for c,ax in zip(time_courses,axes.flatten()):
        data = top_annot_expression_data[time_courses[c]].transpose().apply(lambda x: x*100/sum(x), axis=1)
        data['time_point (h)'] = time_points
        data.set_index('time_point (h)', inplace=True)
        data.plot(kind="bar", title=c, stacked=True, ax=ax, legend=False, fontsize=18 )
    plt.figlegend(data.columns, loc='right', prop={'size': 17})
    fig.suptitle(annot, fontsize=60)

def line_plots(Norm_filter_taxa_data_annot,
               genus_norm_exp_data_annot,
               time_courses, rank, figsize=[8,10],
               name='Raw counts'):
    """ PLOT 1: sum in time for each time course """
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=figsize, sharex=True)

    # raw counts (normalized only by spiked standards)
    tot_sum_dic = {location:Norm_filter_taxa_data_annot[time_courses[location]].sum().values for location in time_courses}
    tot_sum_dic_df = pandas.DataFrame(tot_sum_dic, index= range(0, 72, 4))
    tot_sum_dic_df.plot(title=name, ax=ax1)

    # normalize by location to see the variation
    tot_sum_dic = {location:(Norm_filter_taxa_data_annot[time_courses[location]].sum().values
                             / Norm_filter_taxa_data_annot[time_courses[location]].sum().sum())
                   for location in time_courses}
    tot_sum_dic_df = pandas.DataFrame(tot_sum_dic, index= range(0, 72, 4))
    tot_sum_dic_df.plot(title='Normalized by location', ax=ax2)

    # normalize by rank total
    tot_sum_dic = {location:genus_norm_exp_data_annot[time_courses[location]].sum().values for location in time_courses}
    tot_sum_dic_df = pandas.DataFrame(tot_sum_dic, index= range(0, 72, 4))
    tot_sum_dic_df.plot(title=f'Normalized by {rank}', ax=ax3)

    return fig


if __name__ == "__main__":
    main()
