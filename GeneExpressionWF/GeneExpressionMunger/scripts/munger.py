"""
usage: munger.py [-h] [-f FILES [FILES ...]] [-r INDEX_FILE] [-c] [-o OUTPUT_FILENAME]

optional arguments:
  -h, --help            show this help message and exit
  -f FILES [FILES ...], --files FILES [FILES ...]
                        TCGA Gene Expression TXT files
  -c, --csv
  -r FILE_INDEX         A file with input files, one each line. Merged with -f files.
  -o OUTPUT_FILENAME, --output_filename OUTPUT_FILENAME
"""

import pandas as pd
import argparse
import sys

# To read TXT files:
# df = pd.read_table(filename)
# To read MAF files:
# df = pd.read_table(filename, skiprows=1) # to skip the version row

def get_dataframe_list(args, data_fields=('gene', 'raw_counts')):

    # get a list of dataframes
    dfs, files = [], args['files'] or []
    # create an index using the filenames
    # this will prevent having an overlong command line for 100's or 1000's of files
    if args['file_index']:
        with open(args['file_index']) as fp:
            files.extend(fp.readlines())
    files = sorted(filter(None, set([f.strip() for f in files])))
    # now iterate over the files and get the looooong list of dataframes
    for f in files:
        # Get only specific columns with usecols
        df = pd.read_table(f, usecols=data_fields)
        dfs.append(df)
    return dfs, files # a list of dataframes and the files index

def get_metadata_tag(filename):
    """ Gets a filename (without extension) from a provided path """
    UNCID = filename.split('/')[-1].split('.')[0]
    TCGA = filename.split('/')[-1].split('.')[1]
    return TCGA

def merge_texts(args):
    # get the list of dataframes
    dfs, filenames = get_dataframe_list(args)
    # rename the columns of the first df
    df = dfs[0].rename(columns={'raw_counts': 'raw_counts_' + get_metadata_tag(filenames[0])})
    # enumerate over the list, merge, and rename columns
    for i, frame in enumerate(dfs[1:], 2):
        df = df.merge(frame, on='gene').rename(columns={'raw_counts':'raw_counts_' + get_metadata_tag(filenames[i-1])})
    return df

def get_csv(args, df, filename='GEX_dataframe.csv', header_opt=False, index_opt=False):
    # if csv is true and an output filename is given, rename
    # there is a default filename, so it should pass if --csv is True
    if args['csv'] and args['output_filename']:
        return df.to_csv(path_or_buf=filename, header=header_opt, index=index_opt)

def get_transpose(df):
    df_transpose = df.transpose()
    df_transpose = df_transpose.rename(index = {'gene':'case'})
    return df_transpose

def main(args):
    df = merge_texts(args)
    get_csv(args, df, filename=str(args['output_filename']) + '_by_gene.csv', header_opt=True)
    if args['transpose']:
        get_csv(args, get_transpose(df), filename=str(args['output_filename']) + '_by_case.csv', header_opt=False, index_opt=True)
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", help="TCGA Gene Expression TXT files", nargs="+")
    parser.add_argument("-c", "--csv", action="store_true", default=False)
    parser.add_argument("-t", "--transpose", action="store_true", default=False)
    parser.add_argument("-o", "--output_filename", type=str, default="GEX_dataframe")
    parser.add_argument("-r", "--file_index", type=str, default=None)

    args = parser.parse_args()
    args = vars(args)
    if not args['files'] and not args['file_index']:
        parser.print_help()
        sys.exit(0)

    df = main(args)