#!/usr/bin/env python

import pandas as pd
import sys
import pdb

from code import metadata

def test_validate_filename():
    ''' Tests the function metadata.validate_filename. '''
    sheet_name = 'validate_filename'
    print 'Running test in sheet %s.'%sheet_name
    df = pd.read_excel(tests_file,sheetname=sheet_name,skiprows=2).fillna(False)

    # Iterate through 
    for test_num, row in df.iterrows():
        file_name = row['file_name']
        file_type = row['file_type']
        expected_abort = bool(row['abort'])
        expected_return_value = eval(row['return_value'])

        # Describe test
        print 'test_num: %d'%test_num
        print 'file_type: %s'%file_type
        print 'file_name: %s'%file_name
        print 'expected_abort: %s'%expected_abort
        print 'expected_return_value: %s'%expected_return_value

        # Test abort and return value. If either tests fail, abort. 
        if expected_abort == False:
            return_value = metadata.validate_filename(file_name, file_type)
            assert return_value == expected_return_value,\
                'Return value %s does not match %s.'%\
                (return_value,expected_return_value)     
        elif expected_abort == True:
            try:
                metadata.validate_filename(file_name, file_type)
                abort = False
            except:
                abort = True
            assert abort == True,\
                'Error: Abort expected but did not happend.'

        # If both tests succed, declare success.
        print 'Success!\n'

    # return number of successful tests
    return len(df)

def test_decode_filename():
    ''' Tests the function metadata.decode_filename. '''
    sheet_name = 'decode_filename'
    print 'Running test in sheet %s.'%sheet_name
    df = pd.read_excel(tests_file,sheetname=sheet_name,skiprows=2).fillna(False)

    # Iterate through 
    for test_num, row in df.iterrows():
        file_name = row['file_name']
        file_type = row['file_type']
        expected_abort = bool(row['abort'])
        expected_return_value = eval(row['return_value'])

        # Describe test
        print 'test_num: %d'%test_num
        print 'file_type: %s'%file_type
        print 'file_name: %s'%file_name
        print 'expected_abort: %s'%expected_abort
        print 'expected_return_value: %s'%expected_return_value

        # Test abort and return value. If either tests fail, abort. 
        if expected_abort == False:
            return_value = metadata.decode_filename(file_name, file_type)
            assert return_value == expected_return_value,\
                'Return value %s does not match %s.'%(v,return_value)     
        elif expected_abort == True:
            try:
                metadata.validate_filename(file_name, file_type)
                abort = False
            except:
                abort = True
            assert abort == True,\
                'Error: Abort expected but did not happend.'

        # If both tests succed, declare success.
        print 'Success!\n'

    # return number of successful tests
    return len(df)

if __name__ == '__main__':

    # Get excel files from user
    metadata_file = sys.argv[1]
    group = sys.argv[2]
    tests_file = sys.argv[3]

    # Load metadata
    metadata.load(metadata_file, group=group, use_sample_data=False)

    # Run tests
    num_tests = 0
    num_tests += test_validate_filename()
    num_tests += test_decode_filename()
    print 'All %d tests succeeded!'%num_tests
