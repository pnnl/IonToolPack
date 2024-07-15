import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import h5py
from mza.mza import GetClosestSpectrum, GetExtractedIonRetention, GetExtractedIonArrival

def GenerateMS2plot(mzaFile, outputFilename, molecule, precMz, mztolhalfwidth, rt, fragsMz, fragsIntensity, at=0, mzHalfWindowXIC=0.01, rtrange=0.3, mzrange=0.07, atrange=1.5, minMzDistCentroid = 0.005): # mzrange=0.1
    # Check and flag if ion mobility data:
    isDIAdata = False
    isIMdata = False
    if at > 0: # at = arrival time
        # check if first mza file has ion mobility values:
        with h5py.File(mzaFile, 'r') as mza:
            metadata = mza["Metadata"]
            metadata = metadata[(metadata["MSLevel"] == 2) & (metadata["IsolationWindowTargetMz"] == 0) | (metadata["IsolationWindowLowerOffset"] > 1) & (metadata["IsolationWindowUpperOffset"] > 1)]
            if len(metadata):
                isDIAdata = True
            if len(metadata["IonMobilityBin"] > 0):
                isIMdata = True
    rtvals = []
    xicvals = []
    atvals = []
    atxicvals = []
    msLevels = [1] # To add trace for precursor in first position
    for x in fragsMz: msLevels.append(2)
    legendText = ["p" + str(round(precMz, ndigits=2))]
    for x in fragsMz: legendText.append(str(round(x, ndigits=2)))
    lineStyles = ['--'] # To plot precursor trace dashed
    for x in fragsMz: lineStyles.append('-')
    fragstemp = [precMz]
    for x in fragsMz: fragstemp.append(x)
    for k in range(len(fragstemp)):
        mz = fragstemp[k]
        if isDIAdata:
            df = GetExtractedIonRetention(mzaFileName=mzaFile, mz=mz, msLevel=msLevels[k], startRT=rt-rtrange, endRT=rt+rtrange, mztolhalfwidth=mzHalfWindowXIC)
            if len(df) > 0:
                df.sort_values("rt", inplace=True, ignore_index=True)
                rtvals.append(df["rt"])
                xicvals.append(df["intensity"])
            else:
                rtvals.append([])
                xicvals.append([])
        if isIMdata:
            df = GetExtractedIonArrival(mzaFileName=mzaFile, mz=mz, msLevel=msLevels[k], rt=rt, startAT = at-atrange, endAT = at+atrange, mztolhalfwidth = mzHalfWindowXIC)
            if len(df) > 0:
                df.sort_values("at", inplace=True, ignore_index=True)
                atvals.append(df["at"])
                atxicvals.append(df["intensity"])
            else:
                atvals.append([])
                atxicvals.append([])

    # Scale precursor intensity:
    preci = 0
    if len(xicvals) > 1:
        maxFragIntensity = max(map(lambda x: max(x, default=0), xicvals[1:]), default=0)
        maxPrecIntensity = max(xicvals[preci], default=0)
        # scale precursor intensity to 5% above the mas fragment intensity
        xicvals[preci] = [(x/maxPrecIntensity) * maxFragIntensity for x in xicvals[preci]]
    if len(atxicvals) > 1:
        maxFragIntensity = max(map(lambda x: max(x, default=0), atxicvals[1:]), default=0)
        maxPrecIntensity = max(atxicvals[preci], default=0)
        # scale precursor intensity to 5% above the mas fragment intensity
        atxicvals[preci] = [(x/maxPrecIntensity) * maxFragIntensity for x in atxicvals[preci]]

    # Generate overlaid ion figures: ----------------------------
    # Create a figure with subplots
    npanels = 1
    if isDIAdata:
        npanels = 2
    if isIMdata:
        npanels = 3
    # Spectrum plot:
    minmz = np.min(fragsMz)
    maxmz = np.max(fragsMz)
    [mz_array, intensity_array] = GetClosestSpectrum(mzaFileName=mzaFile, msLevel=2, rt=rt, at=at, precursorMz=precMz, mztolhalfwidth=mztolhalfwidth)
    if len(mz_array) < 2:
        return
    # Normalize intensity:
    maxIntensityExperimental = max(intensity_array)
    intensity_array = [x/maxIntensityExperimental for x in intensity_array]
    maxIntensityReference = max(fragsIntensity)
    fragsIntensity = [x/maxIntensityReference for x in fragsIntensity]
    
    # Generate subplot for the experimental spectrum
    #fig = plt.figure(figsize=(8, 6))
    figheight = 4 * npanels
    fig, axes = plt.subplots(nrows=npanels, ncols=1, figsize=(6, figheight), gridspec_kw={'hspace': 0.5})
    mspanel = axes
    if npanels > 1:
        mspanel = axes[0]

    # Check if spectrum is profile or centroid:
    apexmz = np.array(intensity_array).argmax()
    if apexmz == 0 and len(intensity_array) > 1:
        apexmz+=1
    if apexmz == len(intensity_array) - 1 and len(intensity_array) > 0:
        apexmz-=1
    minMzDist = min(mz_array[apexmz+1] - mz_array[apexmz], mz_array[apexmz] - mz_array[apexmz-1])
    if minMzDist > minMzDistCentroid: # plot centroid
        mspanel.bar(mz_array, intensity_array, width=0.5, label='Experimental')
    else: # plot profile
        mspanel.plot(mz_array, intensity_array, '-', label='Experimental')
    mspanel.set_ylabel("Intensity")
    mspanel.set_title("Intensity/Max: Exp=" + str(round(maxIntensityExperimental)) + ", Ref=" + str(round(maxIntensityReference)))
    fig.suptitle(f"Fragment ions for {molecule}")
    mspanel.legend()
    # Generate subplot for reference spectrum in mirror format
    fragsIntensity = np.array(fragsIntensity) * -1
    mspanel.bar(fragsMz, fragsIntensity, color='red', width=0.5, label='Reference')
    mspanel.legend()
    # Plot horizontal line at 0 on the x-axis
    mspanel.axhline(0, color='black', linestyle='--', linewidth=1)
    # Adjust minimum and maximum x-axis limits for the current molecule
    mspanel.set_xlim(minmz - 50, maxmz + 50)
    mspanel.set_xlabel("m/z")

    if isDIAdata:
        build_overlaid_plot(rtvals, 
                            xicvals, 
                            lineStyles,
                            rt, rtrange, "Retention time", "Intensity", axes[1])
        axes[1].legend(legendText, loc ="right", fontsize=9)

    if isIMdata:
        build_overlaid_plot(atvals, 
                            atxicvals, 
                            lineStyles,
                            at, atrange, "Arrival time", "Intensity", axes[2])
    
    plt.savefig(outputFilename + ".jpg")
    plt.savefig(outputFilename + ".pdf")
    plt.close(fig)


def build_overlaid_plot(x,y, lineStyles, xcenter, xrange, xlabel1, ylabel1, ax):
    for k in range(0,len(x)):
       ax.plot(x[k], y[k], linestyle=lineStyles[k])

    ax.tick_params(axis='x', labelrotation=10) 
    ax.set_xlim([xcenter-xrange, xcenter+xrange*1.2])
    # Plot vertical line at expected on the x-axis
    ax.axvline(xcenter, color='red', linestyle='--', linewidth=0.8)
    ax.set_xlabel(xlabel1)
    ax.set_ylabel(ylabel1)
