import numpy as np
import pandas as pd
from qc.ion_batch import GetHighResCoordinates

def DetectTopmostIons(dfions, dfruns, topIons=4, minMzDistDetectCentroidMS=0.005):
    dfions = dfions.copy()
    dfions["rtRegion"] = 0
    dfions["MZ"] = dfions["MZ"].astype(int)
    dftopmost = pd.DataFrame()
    Nrt = int(topIons) # take more blocks across RT because the mass is usually very stable and good
    selectedMzs = [] # to keep only the first MZ, each ion will be extracted from all runs
    for label in dfions['LABELSAMPLEGROUP'].unique():
        df = dfions[(dfions["LABELSAMPLEGROUP"] == label) & (~dfions["MZ"].isin(selectedMzs))]
        if len(df) == 0:
            continue
        # Step 1: Partition the space
        rtbins = pd.cut(df['RT'], bins=Nrt, labels=False)
        # Step 2: Select rows with maximum frequency and intensity in each zone
        rtRegion = 1
        for rt_bin in range(Nrt):
            zone_rows = df[(rtbins == rt_bin) & (~df["MZ"].isin(selectedMzs))]
            if not zone_rows.empty:
                zone_rows = zone_rows[zone_rows['FREQ'] == max(zone_rows['FREQ'])]
                #max_row = zone_rows.loc[(zone_rows['FREQ'] * zone_rows['INTENSITY']).idxmax()]
                max_row = zone_rows.loc[(zone_rows['INTENSITY']).idxmax()]
                df.at[max_row.name, 'rtRegion'] = rtRegion
                rtRegion = rtRegion + 1
                selectedMzs.append(df.loc[max_row.name, "MZ"])
        dftopmost = pd.concat([dftopmost, df[df["rtRegion"] > 0]])
    
    dftopmost.reset_index(inplace=True, drop=True)
    # keep only the first MS run for each LABELSAMPLEGROUP
    dfruns = dfruns.copy()
    dfruns.drop_duplicates(['LABELSAMPLEGROUP'], inplace=True, keep='first')
    dfruns.reset_index(inplace=True, drop=True)
    for k in dftopmost.index:
        ionmz = dftopmost["MZ"][k]
        ionrt = dftopmost["RT"][k]/10 # <- ToDo: correct scaling in PCA.py
        mzaFile = list(dfruns["MZAPATH"][dfruns["LABELSAMPLEGROUP"] == (dftopmost["LABELSAMPLEGROUP"][k])])[0] + ".mza"
        [ionmz,ionrt] = GetHighResCoordinates(mzaFile, ionmz, ionrt, rtrange=1, mzrange=1, minMzDistCentroid=minMzDistDetectCentroidMS) # these tolerances must be kept at unit resolution
        dftopmost.loc[k,"MZ"] = ionmz
        dftopmost.loc[k,"RT"] = ionrt 
    
    dftopmost = dftopmost[(dftopmost["MZ"] > 0) & (dftopmost["RT"] > 0)]
    dftopmost.sort_values(by=["LABELSAMPLEGROUP", "rtRegion", "FREQ", "INTENSITY"], ascending=[True, True, False, False], inplace=True)
    dftopmost.drop_duplicates(subset=["LABELSAMPLEGROUP", "rtRegion"], inplace=True, keep='first')
    dftopmost.drop(columns=["rtRegion"], inplace=True)
    dftopmost.reset_index(inplace=True, drop=True)
    dftopmost["MOLECULE"] = ["Ion" + str(k+1) for k in dftopmost.index]
    #dftopmost["MOLECULE"] += "-MZ" + str(round(dftopmost["MZ"], ndigits=2))+ "-RT" + str(round(dftopmost["RT"], ndigits=1))
    return dftopmost
