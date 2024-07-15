import pandas as pd
import os

RawFormats = [".mza", ".raw", ".d", ".mzML"]

def FormatDataframeSamples(df, basePathCsvRuns=""):

    # format column names: LABELSAMPLEGROUP, MSRUNID, MSRUN, MSRUNFORMAT, MSRUNPATH
    df.columns = df.columns.str.upper()
    df.columns = df.columns.str.replace(' ', '')
    df.columns = df.columns.str.replace('_', '')
    df.columns = df.columns.str.replace('DATASETFOLDERPATH', 'MSRUNPATH')
    df.columns = df.columns.str.replace('DATASET', 'MSRUN')

    # check if acq_start column exists, otherwise generate it
    if "ACQSTART" not in df.columns:
        df["ACQSTART"] = range(1, len(df)+1)

    df["ACQSTART"] = pd.to_datetime(df['ACQSTART'], errors='ignore')
    # sort data by ACQSTART
    df = df.sort_values(by=["ACQSTART"], ascending=[True], ignore_index=True)

    # check if LABELSAMPLEGROUP column exists, otherwise generate it
    if "LABELSAMPLEGROUP" not in df.columns:
        df["LABELSAMPLEGROUP"] = ""
        for i, row in df.iterrows():
            if "blank" in str(row.MSRUN).lower():
                df.at[i,"LABELSAMPLEGROUP"] = "Blank"
            elif "qc" in str(row.MSRUN).lower():
                df.at[i,"LABELSAMPLEGROUP"] = "QC"
            else:
                df.at[i,"LABELSAMPLEGROUP"] = "Sample"

    # check if MSRUNID column exists and has unique values, otherwise generate it
    if "MSRUNID" in df.columns and len(set(df["MSRUNID"])) != len(df):
        df["index"] = range(1, len(df)+1)
        df["MSRUNID"] = df["MSRUNID"].astype(str) + '_' + df["index"].astype(str)
        df.drop(columns=["index"])
    if "MSRUNID" not in df.columns:
        df["MSRUNID"] = range(1, len(df)+1)

    # Extract MS data format as separate column:
    if "MSRUNFORMAT" not in df.columns:
        df["MSRUNFORMAT"] = ""
    for i, row in df.iterrows():
        runx = row.MSRUN
        for x in RawFormats:
            if runx.endswith(x): # if provided in the raw file name
                runx = runx.removesuffix(x)
            # check if file exists with path provided as absolute:    
            xpath = os.path.join(row.MSRUNPATH, runx + x)
            if os.path.exists(xpath):
                df.at[i,"MSRUN"] = runx
                df.at[i,"MSRUNFORMAT"] = x
                break
            # check if file exists with path provided as relative
            xpath = os.path.join(basePathCsvRuns, row.MSRUNPATH, runx + x)
            if os.path.exists(xpath):
                df.at[i,"MSRUN"] = runx
                df.at[i,"MSRUNFORMAT"] = x
                df.at[i,"MSRUNPATH"] = os.path.join(basePathCsvRuns, row.MSRUNPATH)
                break
    
    notFound = df[(df["MSRUNFORMAT"] == "")]
    if len(notFound) > 0:
        print("MS run not found:")
        for i, row in notFound.iterrows():
            print("     " + os.path.join(row.MSRUNPATH, row.MSRUN))

    # Keep only MS runs found:
    df = df[(df["MSRUNFORMAT"] != "")]
    if len(df) == 0:
        raise Exception("No valid raw MS data files found.")
        
    df = df[["LABELSAMPLEGROUP", "MSRUNID", "MSRUN", "MSRUNPATH", "MSRUNFORMAT", "ACQSTART"]]
    df["MSRUNPATH"] = [os.path.normpath(x) for x in df["MSRUNPATH"]]

    return df

