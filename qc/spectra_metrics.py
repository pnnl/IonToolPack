import h5py
import hdf5plugin
import pandas as pd

def ExtractSpectraMetadataMetrics(mzaFile):
    with h5py.File(mzaFile, 'r') as mza:
        metadata = mza["Metadata"]
        # Convert metadata to a DataFrame
        metadata = pd.DataFrame(metadata[:]) 

    # Calculate metrics spectra summary statistics
    myColumns = ['COUNT', 'MEANTIC', 'MEDIANTIC', 'MAXTIC']
    df = pd.DataFrame()
    # Group by 'MSLevel' column and extract 'TIC' column
    grouped_tic = metadata.groupby('MSLevel')['TIC']

    # Iterate over each group and make calculations
    for mslevel, tic_group in grouped_tic:
        count = tic_group.count()
        mean_tic = tic_group.mean()
        median_tic = tic_group.median()
        max_tic = tic_group.max()
        df['MS' + str(mslevel) + myColumns[0]] = [count]
        df['MS' + str(mslevel) + myColumns[1]] = [mean_tic]
        df['MS' + str(mslevel) + myColumns[2]] = [median_tic]
        df['MS' + str(mslevel) + myColumns[3]] = [max_tic]

    return df

# testdf = ExtractSpectraMetadataMetrics("E:/QCdecoder/code/test_data/LC-IM-MS-Agilent/DataMza/Agile_Pput_CJ019_WT_IMS_5ul_R1_01_20V_redo_12Sep20_Kristin_MA-d3-c3-Min20.mza")
# print(testdf)
