from __future__ import division
import os
import shutil
import pandas as pd
import numpy as np
import utils
import re
import pdb
from glob import glob

IS_LOADED = False

def load(metadata_file, group, use_sample_data):
    '''
    Loads metadata from metadata_file for use in all other pipeline functions. Must be called before any other functions are called. 
    '''

    # Stores metadata as dictionary of dataframes. 
    global METADATA

    # Stores the METADATA filename
    global FILENAME
    FILENAME = metadata_file
    global USE_SAMPLE_DATA
    USE_SAMPLE_DATA = use_sample_data

    # Stores information about the file naming standard
    global FILETYPE_TO_DIR
    global FILETYPE_TO_REGEX
    global FILETYPE_TO_TEMPLATE
    global ELEMENT_TO_VALUES
    global FILETYPE_TO_GLOB
    global COLUMN_TO_DATATYPE
    global IS_LOADED

    assert os.path.isfile(metadata_file), \
        'Error: Cant find metadata file %s.'%metadata_file

    # Clear existing metadata
    METADATA = {}
    FILETYPE_TO_DIR = {}
    FILETYPE_TO_REGEX = {}
    FILETYPE_TO_TEMPLATE = {}
    FILETYPE_TO_GLOB = {}
    ELEMENT_TO_VALUES = {}

    # Load excel file
    excel_file = pd.ExcelFile(metadata_file)

    # Check for 'column_datatypes' sheet. If found, use to set data types
    # of specified columns
    sheet = 'column_datatypes'
    COLUMN_TO_DATATYPE = {}
    if sheet in excel_file.sheet_names:
        df = excel_file.parse(sheet,skiprows=2).fillna(False)
        for i, row in df.iterrows():
            COLUMN_TO_DATATYPE[row['column']] = eval(row['data_type'])

    # Load sheets into METADATA dict
    for sheet_name in excel_file.sheet_names:

        # Excel often adds this sheet without consent. Ignore it.
        if sheet_name=='Sheet2':
            continue

        # Get sheet as df
        df = excel_file.parse(sheet_name,skiprows=2).fillna(False)

        # Convert columns to datatype if listed in 'column_datatypes' sheet. 
        for col in df.columns:
            if col in COLUMN_TO_DATATYPE:
                df.loc[:,col] = df[col].astype(COLUMN_TO_DATATYPE[col])
        
        # Only use valid rows of df
        if 'use' in df.columns:
            df = df[df['use']]

        # Save dataframe in METADATA dict
        METADATA[sheet_name] = df


    ######
    # Enforce constraints and set attributes
    global GROUP
    global GROUPS
    global SAMPLES
    global LIBRARIES
    global BC_SAMPLES
    global JCT_SAMPLES
    global LIDS
    global AMPLICONS
    global READS

    def filter_sheet(sheet,column,values):
        ''' 
        Removes rows of sheet for which the entries in a 
        specified column do not appare in a provides set of values.

        sheet = metadata sheet name
        column = metadata column name
        values = list of acceptable values
        returns nothing. 
        '''

        # Make sure df has the specified column
        assert sheet in METADATA, \
            'Error: sheet %s not in METADATA.'
        df = METADATA[sheet]

        # Make sure all values are in column of df
        assert set(values) <= set(df[column])

        # Remove rows of df that if df[column] is not in values
        indices = [v in values for v in df[column]]
        tmp_df = df[indices]

        # Make sure that the df was not entirely annihilated
        assert len(tmp_df) > 0, \
            'Error: no filtered rows found in sheet=%s, column=%s'%\
            (sheet, column, values) 
        METADATA[sheet] = tmp_df


    def column_union(sheets,columns):
        '''
        Returns the set of values found in a set of columns across a set of
        sheets. 
        '''

        # If a  string was provided for sheets or columns, make a list instead
        if type(sheets)==str:
            sheets = [sheets]
        if type(columns)==str:
            columns = [columns]

        # Tally values in all specified columns of all specified sheets
        values = []
        for sheet in sheets:

            # Make sure sheet exists
            assert sheet in METADATA, \
                'Error: sheet %s does not exist.'%sheet

            for column in columns:
                # Make sure sheet has column
                assert column in METADATA[sheet].columns, \
                    'Error: column %s is not in sheet %s.'%(column, sheet)

                # Record unique values in column
                values += list(set(METADATA[sheet][column].values))

        # Remove False values from values
        values = [v for v in values if not v is False]

        # Make sure all values are unique
        assert len(values)==len(set(values)), \
            'Error: extracted values are not unique.'

        # Return the set of accumulated values. 
        return values

    # If group is 'all', then use all experiments
    GROUP = group
    if group == 'all':
        GROUPS = list(set(METADATA['experiments']['group'].values))
    else:
        GROUPS = [group]

    # Filter experiments sheet
    filter_sheet(sheet='experiments',column='group',values=GROUPS)
    
    # Set attributes from experiments sheet
    SAMPLES = column_union(
        sheets='experiments',
        columns=['library','total_bc_sample','exon_bc_sample','jct_sample'])
    LIBRARIES = column_union(
        sheets='experiments',
        columns='library')
    BC_SAMPLES = column_union(
        sheets='experiments',
        columns=['total_bc_sample','exon_bc_sample']) 
    JCT_SAMPLES = column_union(
        sheets='experiments',
        columns='jct_sample')

    # Filter next layer of sheets
    filter_sheet(sheet='samples',column='sample',values=SAMPLES)
    filter_sheet(sheet='libraries',column='sample',values=LIBRARIES)
    filter_sheet(sheet='bc_samples',column='sample',values=BC_SAMPLES)
    filter_sheet(sheet='jct_samples',column='sample',values=JCT_SAMPLES)

    # Set next layer of attributes
    LIDS = column_union(
        sheets='samples',
        columns='LID')
    AMPLICONS = column_union(
        sheets=['libraries','bc_samples','jct_samples'],
        columns='amplicon')

    # Filter next layer of sheets
    filter_sheet(sheet='illumina_runs',column='LID',values=LIDS)
    filter_sheet(sheet='amplicons',column='amplicon',values=AMPLICONS)

    # Set last layer of attributes
    READS = column_union(
        sheets='amplicons',
        columns=['read1','read2'])

    # Fitler last layer of sheets
    filter_sheet(sheet='reads',column='read',values=READS)

    ######

    # If using sample data
    if use_sample_data:

        # Alter file names in illumina_runs sheet
        df = METADATA['illumina_runs']
        for index, row in df.iterrows():
            df.loc[index,'read1_file'] = 'sample.'+row['read1_file']
            df.loc[index,'read2_file'] = 'sample.'+row['read2_file']

        # Also alter associated patterns in pipeline_filetypes
        df = METADATA['pipeline_filetypes']
        for index, row in df.iterrows():
            if row['file_type'] == 'illumina_run':
                df.loc[index,'glob'] = 'sample.'+row['glob']
                df.loc[index,'regex'] = 'sample\.'+row['regex']

    # Set global variables containing valid file information
    for i, row in METADATA['pipeline_filetypes'].iterrows():
        file_type = str(row['file_type'])
        
        FILETYPE_TO_DIR[file_type] = str(row['directory'])

        FILETYPE_TO_REGEX[file_type] = str(row['regex'])

        FILETYPE_TO_GLOB[file_type] = str(row['glob'])

        # Add to FILETYPE_TO_TEMPLATE only if there is a template
        template = str(row['template'])
        if template:
            elements = re.split('\W+',str(row['elements']))
            FILETYPE_TO_TEMPLATE[file_type] = (template,elements)

    # Load valid values for columns
    ELEMENT_TO_VALUES = {}
    df = METADATA['pipeline_filename_elements']
    for i, row in df.iterrows():
        key = row['element']
        if row['unrestricted']:
            continue
        elif row['use_values']:
            values_str = re.split('\W+',row['values'])
            indices = (METADATA['column_datatypes']['column']==row['column'])
            assert sum(indices)==1,\
                'Error: cant identify unique row in column_datatypes sheet for key %s'%key
            values_type = eval(\
                METADATA['column_datatypes']['data_type'][indices].values[0]
                )
            values = [values_type(v) for v in values_str]
        elif row['use_sheet']:
            sheet_name = row['sheet']
            column  = row['column']
            values = METADATA[sheet_name][column].values
        ELEMENT_TO_VALUES[key] = values

    # Get dictionary mapping samples to read1 and samples to read2
    global SAMPLE_TO_READ1
    global SAMPLE_TO_READ2
    global READ_TO_REGEX
    SAMPLE_TO_READ1 = {}
    SAMPLE_TO_READ2 = {}
    READ_TO_REGEX = {}
    samples_df = METADATA['samples'].set_index('sample')
    amplicons_df = METADATA['amplicons'].set_index('amplicon')
    reads_df = METADATA['reads'].set_index('read')
    for sample in samples_df.index.values:

        # Get amplicon corresponding to sample
        amplicon = samples_df.loc[sample,'amplicon']

        # Get reads corresponding to sample
        SAMPLE_TO_READ1[sample] = amplicons_df.loc[amplicon,'read1']
        SAMPLE_TO_READ2[sample] = amplicons_df.loc[amplicon,'read2']

    for read in reads_df.index.values:
        READ_TO_REGEX[read] = re.compile(reads_df.loc[read,'regex'])

    # Register metadata as loaded
    IS_LOADED = True

def parse_read(read_name,sequence):
    assert IS_LOADED, 'Error: metadata not yet loaded.'

    regex = READ_TO_REGEX[read_name]

    m = re.match(regex,sequence)

    # If parse failed, return False
    if not m:
        return False

    # Otherwise, return out_dict with all possible elements
    else:
        raw_dict = m.groupdict()
        out_dict = {}
        for key in ['ss','bc','jct']:
            key_rc = key+'_rc'
            if key in raw_dict:
                out_dict[key]=str(raw_dict[key])
            elif key_rc in raw_dict:
                out_dict[key]=utils.rc(str(raw_dict[key_rc]))
            else:
                out_dict[key]=None
        return out_dict


def parse_read_given_sample(sample,read_num,sequence):
    ''' Parses a read, returns a dict containing extracted elements. '''
    assert IS_LOADED, 'Error: metadata not yet loaded.'
    assert read_num in [1,2],\
        'Error: invalid read_num %d'%read_num

    assert sample in SAMPLE_TO_READ1,\
        'Error: sample %s not in SAMPLE_TO_READ1'%sample

    assert sample in SAMPLE_TO_READ2,\
        'Error: sample %s not in SAMPLE_TO_READ2'%sample

    # Get read name
    read = SAMPLE_TO_READ1[sample] if read_num==1 else SAMPLE_TO_READ2[sample]

    # Return parsed read
    return parse_read(read,sequence)

def get_splits_from_LID(LID):
    assert IS_LOADED, 'Error: metadata not yet loaded.'
    ''' Returns splits (in list) corresponding to a given LID'''
    read_files = get_existing_filenames('split_fastq',{'LID':LID})
    splits = list({decode_filename(f,'split_fastq')['split'] \
        for f in read_files})
    return splits

def get_filename_and_group():
    assert IS_LOADED, 'Error: metadata not yet loaded.'
    return FILENAME, GROUP

def get_use_sample_data():
    assert IS_LOADED, 'Error: metadata not yet loaded.'
    return USE_SAMPLE_DATA

def get_all_element_values(element):
    ''' Returns a list of all valid values for a specified element '''
    assert IS_LOADED, 'Error: metadata not yet loaded.'

    assert element in ELEMENT_TO_VALUES,\
        'Error: element %s is not found in ELEMENT_TO_VALUES.'%element
    return ELEMENT_TO_VALUES[element]

def get_sheet(sheet_name):
    ''' Returns specified sheet as a pandas dataframe. '''
    assert IS_LOADED, 'Error: metadata not yet loaded.'

    assert sheet_name in METADATA,\
        'Error: sheet_name %s not found in METADATA.'%sheet_name
    return METADATA[sheet_name].copy()

def filename_to_report_filename(file_name):
    ''' Returns the name of the report file corresponding to file_name'''
    assert IS_LOADED, 'Error: metadata not yet loaded.'

    directory, basename = os.path.split(file_name)
    report_file = '%s/report.%s'%(directory,basename)
    return report_file

def get_filetype_directory(file_type):
    ''' Returns a directory from a given file_type '''
    assert IS_LOADED, 'Error: metadata not yet loaded.'

    assert file_type in FILETYPE_TO_DIR,\
        'Error: unrecognized file_type %s.'%file_type
    return FILETYPE_TO_DIR[file_type]

def get_existing_filenames(file_type, elements_dict={}):
    ''' Returns a list of all existing files for specified file_type.
    Only returns files whose names have specific elements'''
    assert IS_LOADED, 'Error: metadata not yet loaded.'

    directory = FILETYPE_TO_DIR[file_type]
    g = FILETYPE_TO_GLOB[file_type]
    all_file_names = glob('%s/%s'%(directory,g))

    # Only keep file names that have elements matching in elements_dict
    keep_file_names = []
    for file_name in all_file_names:
        keep = True
        try:
            file_dict = decode_filename(file_name,file_type)
            for key, value in elements_dict.iteritems():
                if file_dict[key]!=value:
                    keep = False
        except:
            keep = False

        if keep:
            keep_file_names.append(file_name)

    return keep_file_names

def validate_element(element,value):
    '''
    Checks to make sure that an element parsed from a sequence is valid. 
    Raises exception if not. 
    '''
    assert IS_LOADED, 'Error: metadata not yet loaded.'

    df = METADATA['pipeline_filename_elements']
    assert element in df['element'].values,\
        'Error: unrecognized element %s'%element

    row = df[df['element'] == element].squeeze()

    if row['use_values']:
        values = re.split('\W+',row['values'])

    elif row['use_sheet']:
        values = METADATA[row['sheet']][row['column']].values

    # If necessary, change datatype of values extracted from METADATA
    if element in COLUMN_TO_DATATYPE:
        datatype = COLUMN_TO_DATATYPE[element]
        values = [datatype(v) for v in values]

    if not row['unrestricted']:
        assert value in values,\
            'Error: invalid value %s for element %s.'%(element,value)

def validate_filename(file_name, file_type, get_elements=False):
    ''' 
    Validates whether a file_name adheres to conventions.  
    Raises exception if not. 
    Able to return elements in a given file name if requested.
    '''
    assert IS_LOADED, 'Error: metadata not yet loaded.'

    directory, basename = os.path.split(file_name)

    # check directory
    assert file_type in FILETYPE_TO_DIR,\
        'Error: file_type %s not found in FILETYPE_TO_DIR.'%file_type
    
    assert directory == FILETYPE_TO_DIR[file_type],\
        'Error: directory %s is invalid for file_type %s'%\
        (directory, file_type)

    # check basename
    assert file_type in FILETYPE_TO_REGEX,\
        'Error: file_type %s not found in FILETYPE_TO_REGEX.'%file_type
    
    # Extract elements
    pattern = FILETYPE_TO_REGEX[file_type]
    p = re.compile(pattern)
    m = re.match(pattern,basename)
    assert bool(m),\
        'Error: file basename %s is invalid for file_type %s'%\
        (basename, file_type) 
    elements_dict = m.groupdict()

    # Chage element datatype if necessary
    for key in elements_dict:
        if key in COLUMN_TO_DATATYPE:
            datatype = COLUMN_TO_DATATYPE[key]
            elements_dict[key] = datatype(elements_dict[key])

    # Verify elements_dict contains all necessary keys
    if file_type in FILETYPE_TO_TEMPLATE:
        xxx, elements = FILETYPE_TO_TEMPLATE[file_type]
        assert set(elements)==set(elements_dict.keys()),\
            'Error: elements %s do not match elements_dict %s'%\
            (str(elements),str(elements_dict))

     # Validate all elements
    for element, value in elements_dict.iteritems():
        validate_element(element,value)

    # Return elements_dict if requested
    if get_elements:
        return elements_dict

def decode_filename(file_name, file_type):
    '''
    Parses informative elements out of a filename. 
    Note: just a wrapper for validate_filename
    '''
    assert IS_LOADED, 'Error: metadata not yet loaded.'

    # Return dictionary, e.g. {'LID': '296963', 'read_num': '2'}
    return validate_filename(file_name, file_type, get_elements=True)


def encode_filename(elements_dict, file_type, report=False):
    '''
    Creates the valid file name containing the information in elements_dict.
    '''
    assert IS_LOADED, 'Error: metadata not yet loaded.'

    # Get directory
    assert type(file_type) != dict,\
        'Error: looks like the order of arguments to encode_filename is wrong.'
    assert file_type in FILETYPE_TO_DIR,\
        'Error: file_type %s is not in FILETYPE_TO_DIR.'%file_type
    directory = FILETYPE_TO_DIR[file_type]

    # Get template and elements
    assert file_type in FILETYPE_TO_TEMPLATE,\
        'ERROR: file_type %s is not in FILETYPE_TO_TEMPLATE.'%file_type
    template, elements = FILETYPE_TO_TEMPLATE[file_type]

    # Make sure user has supplied correct info in elements_dict
    assert set(elements)==set(elements_dict.keys()),\
        'Error: elements %s do not match elements_dict %s'%\
        (str(elements),str(elements_dict))

    # Create file basename
    basename = template%tuple([elements_dict[e] for e in elements])

    # Create full file name
    file_name = '%s/%s'%(directory,basename)

    # Verify validity of file name
    validate_filename(file_name,file_type)

    # Create report file_name
    if report:
        file_name = '%s/report.%s'%(directory,basename)

    return file_name