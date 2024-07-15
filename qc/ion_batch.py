import os
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mza.mza import GetExtractedIonRetention, GetClosestSpectrum, GetExtractedIonArrival
from scipy.signal import find_peaks

def GetHighResCoordinates(mzaFile, mz, rt, rtrange=0.5, mzrange=0.5, minMzDistCentroid = 0.005):
    df = GetExtractedIonRetention(mzaFileName=mzaFile, mz=mz, msLevel=1, startRT=rt - rtrange, endRT=rt + rtrange, mztolhalfwidth=mzrange/2)
    mz_array = []
    if len(df) > 0:
        df.sort_values("rt", inplace=True, ignore_index=True)
        # Get spectrum at apex:
        index = find_peaks(df["intensity"], width=3)[0]
        if len(index) > 0:
            #index = index[np.array(abs(df["RT"][index] - rt)).argmin()]
            index = index[np.array(df["intensity"][index]).argmax()] # keep the most intense
        else:
            return [0, 0]
        rtapex = df.loc[index,"rt"]

        [mz_array, intensity_array] = GetClosestSpectrum(mzaFileName=mzaFile, msLevel=1, rt=rtapex)
        indexes = np.argwhere((mz_array >= mz - mzrange) & (mz_array <= mz + mzrange))
        mz_array = mz_array[indexes].flatten() # keep only the selected indexes
        intensity_array = intensity_array[indexes].flatten()

    if len(mz_array) == 0:
        return [0, 0]
    
    # Check if spectrum is profile or centroid:
    apexmz = np.array(intensity_array).argmax()
    if apexmz == 0:
        apexmz+=1
    if apexmz == len(intensity_array) - 1:
        apexmz-=1
    minMzDist = min(mz_array[apexmz+1] - mz_array[apexmz], mz_array[apexmz] - mz_array[apexmz-1])
    if minMzDist < minMzDistCentroid: # profile
        apexmz = find_peaks(intensity_array, width=3)[0]
        if len(apexmz) > 0:
            #apexmz = apexmz[np.array(abs(mz_array[apexmz] - mz)).argmin()]
            apexmz = apexmz[np.array(intensity_array[apexmz]).argmax()] # keep the most intense
        else:
            return [0, 0]
    
    return [mz_array[apexmz], rtapex]


def GenerateImageIonBatch(dfruns, outputFolder, mz, rt, molecule, suffixImage, mzHalfWindowXIC=0.01, rtrange=0.3, mzrange=0.075, at=0, atrange=1.5):
    # Check and flag if ion mobility data:
    isIMdata = False
    if at > 0: # at = arrival time
        # check if first mza file has ion mobility values:
        mzaFile = dfruns.loc[0,"MZAPATH"] + ".mza" 
        with h5py.File(mzaFile, 'r') as mza:
            metadata = mza["Metadata"]
            if len(metadata["IonMobilityBin"] > 0):
                isIMdata = True
    mzvals = []
    flagIsMzCentroid = []
    intensvals = []
    rtvals = []
    xicvals = []
    atvals = []
    atxicvals = []
    
    for mzaFile in dfruns["MZAPATH"]:
        mzaFile = mzaFile + ".mza"
        df = GetExtractedIonRetention(mzaFileName=mzaFile, mz=mz, msLevel=1, startRT=rt-rtrange, endRT=rt+rtrange, mztolhalfwidth=mzHalfWindowXIC)
        if len(df) > 0:
            df.sort_values("rt", inplace=True, ignore_index=True)
            rtvals.append(df["rt"])
            xicvals.append(df["intensity"])
            # Get spectrum at apex:
            index = find_peaks(df["intensity"], width=3)[0]
            if len(index) > 0:
                index = index[np.array(df["intensity"][index]).argmax()] # keep the most intense
            else:
                # Get spectrum at apex:
                index = np.array(df["intensity"]).argmax()
            [mz_array, intensity_array] = GetClosestSpectrum(mzaFileName=mzaFile, msLevel=1, rt=df.loc[index,"rt"])
            indexes = np.argwhere((mz_array >= mz - mzrange) & (mz_array <= mz + mzrange))
            mz_array = mz_array[indexes].flatten() # keep only the selected indexes
            intensity_array = intensity_array[indexes].flatten()
            mzvals.append(mz_array)
            intensvals.append(intensity_array)
        else:
            rtvals.append([])
            xicvals.append([])
            mzvals.append([])
            intensvals.append([])
        
        if isIMdata:
            df = GetExtractedIonArrival(mzaFileName=mzaFile, mz=mz, msLevel=1, rt=rt, startAT = at-atrange, endAT = at+atrange, mztolhalfwidth = mzHalfWindowXIC)
            if len(df) > 0:
                df.sort_values("at", inplace=True, ignore_index=True)
                atvals.append(df["at"])
                atxicvals.append(df["intensity"])
            else:
                atvals.append([])
                atxicvals.append([])
    
    # Generate overlaid ion figures: ----------------------------
    # Create a figure with subplots
    npanels = 2
    if isIMdata:
        npanels = 3
    for label in dfruns["LABELSAMPLEGROUP"].unique():
        indexes = dfruns[dfruns["LABELSAMPLEGROUP"] == label].index
        suffixImageNew = suffixImage.replace("MZ", label + "-MZ")
        legendText = dfruns.loc[indexes, "legend"]
        if len(legendText) > 15:
            legendText = legendText[0:14]

        figheight = 4 * npanels
        fig, axes = plt.subplots(nrows=npanels, ncols=1, figsize=(6, figheight), gridspec_kw={'hspace': 0.5})
        
        build_overlaid_plot([mzvals[i] for i in indexes.values],
                            [intensvals[i] for i in indexes.values], 
                            mz, mzrange, "$\it{m/z}$", "Intensity", axes[0])
        
        build_overlaid_plot([rtvals[i] for i in indexes.values], 
                            [xicvals[i] for i in indexes.values], 
                            rt, rtrange, "Retention time", "Intensity", axes[1])
        
        if isIMdata:
            build_overlaid_plot([atvals[i] for i in indexes.values], 
                                [atxicvals[i] for i in indexes.values], 
                                at, atrange, "Arrival time", "Intensity", axes[2])

        plt.legend(legendText, loc ="right", fontsize=9)
        fig.suptitle(suffixImageNew)
        plt.savefig(os.path.join(outputFolder, suffixImageNew + ".jpg"), bbox_inches="tight")
        plt.savefig(os.path.join(outputFolder, suffixImageNew + ".pdf"), bbox_inches="tight")
        plt.close()

    # Compute errors: ----------------------------
    df = pd.DataFrame({"MSRUN": dfruns["MZAPATH"], 
                       "LABELSAMPLEGROUP": dfruns["LABELSAMPLEGROUP"],
                       "MSRUNID": dfruns["MSRUNID"],
                       "MOLECULE": molecule,
                       "MZTARGET": mz, 
                       "RTTARGET": rt, 
                       "MZ": -1, 
                       "RT": -1, 
                       "ABUNDANCE": 0})
    if isIMdata:
        df["ATTARGET"] = at
        df["AT"] = -1
    df["MSRUN"] = [os.path.basename(x).replace(".mza", "") for x in df["MSRUN"]]
    for k in range(0, len(mzvals)):
        # get local maxima per array
        if len(intensvals[k]) > 0:
            apexmz = find_peaks(intensvals[k], width=3, rel_height=0.8)[0]
            if len(apexmz) > 0:
                apexmz = apexmz[np.array(abs(np.array(mzvals[k])[apexmz] - mz)).argmin()] # keep the smallest error apex
                df.loc[k, "MZ"] = mzvals[k][apexmz]

        if len(xicvals[k]) > 0:
            # ToDo: save peak metrics to integrate only fwhm as abundance
            apexrt = find_peaks(xicvals[k], width=3, rel_height=0.8)[0]
            if len(apexrt) > 0:
                # apexrt = apexrt[np.array(abs(np.array(rtvals[k])[apexrt] - rt)).argmin()] # keep the smallest error apex
                apexrt = apexrt[np.array(xicvals[k][apexrt]).argmax()] # keep the most intense apex
                df.loc[k, "RT"] = rtvals[k][apexrt]
                df.loc[k, "ABUNDANCE"] = sum([x/1000 for x in xicvals[k]]) # scale intensity to avoid overflow
        
        if isIMdata and len(atxicvals[k]) > 0:
            apexat = find_peaks(atxicvals[k], width=3, rel_height=0.8)[0]
            if len(apexat) > 0:
                #apexat = apexat[np.array(abs(np.array(atvals[k])[apexat] - at)).argmin()] # keep the smallest error apex
                apexat = apexat[np.array(atxicvals[k][apexat]).argmax()] # keep the most intense apex
                df.loc[k, "AT"] = atvals[k][apexat]

    if np.argwhere(df["MZ"] == -1).size > 0:
        df.loc[df["MZ"] == -1, "MZ"] = np.nan
    if np.argwhere(df["RT"] == -1).size > 0:
        df.loc[df["RT"] == -1, "RT"] = np.nan

    df["MZERROR"] = df["MZ"] - mz
    df["MZERRORPPM"] = (df["MZ"] - mz) / mz * 1e6
    df["RTERROR"] = df["RT"] - rt

    if isIMdata:
        if np.argwhere(df["AT"] == -1).size > 0:
            df.loc[df["AT"] == -1, "AT"] = np.nan
        df["ATERROR"] = df["AT"] - at

    df = df[~(np.isnan(df["MZ"]) | np.isnan(df["RT"]))]
    return df


def build_overlaid_plot(x,y, xcenter, xrange, xlabel1, ylabel1, ax):
    for k in range(0,len(x)):
        ax.plot(x[k], y[k], '-')
    # format the x tick labels
    #ticks_loc = ax.get_xticks()
    #ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    #ax.set_xticklabels(['{:.3f}'.format(x) for x in ticks_loc])
    ax.tick_params(axis='x', labelrotation=10) 
    ax.set_xlim([xcenter-xrange, xcenter+xrange*1.2]) # Add a some room on the right side for the legend
    # Plot vertical line at expected on the x-axis
    ax.ticklabel_format(axis='x', style='plain')
    ax.axvline(xcenter, color='red', linestyle='--', linewidth=0.8)
    ax.set_xlabel(xlabel1)
    ax.set_ylabel(ylabel1)

