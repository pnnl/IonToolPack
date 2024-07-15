import h5py
import hdf5plugin
import numpy as np
from PIL import Image 
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# GenerateImageTimeVsMz
# m/z dimension is rounded and summed, used at unit resolution
# retention time is rounded to 1 decimal places and summed

def GenerateImageTimeVsMz(mzaFile, outputFullPath, LcmsImageMinIntensityPercentage=10, LcmsImageMaxIntensityCeilingPercentage=70):

    mza = h5py.File(mzaFile, 'r')
    # Reading Metadata table:
    metadata = mza["Metadata"]
    # Keep only MS1, and TFS for ion mobility data
    metadata = metadata[(metadata["MSLevel"] == 1) & (metadata["IonMobilityBin"] == 0)]
    metadata = metadata[np.argsort(metadata["RetentionTime"])]

    lcms = {} # Create a dictionary with keys as RT and values as mz dictionaries
    maxMz = 0
    maxIntensity = 0
    for k in range(0, metadata.size):
        scan = metadata["Scan"][k]
        mzapath = str(metadata["MzaPath"][k], 'utf-8')

        mz_array = None
        if "Full_mz_array" in mza:
            # Get array of m/z values (common for all spectra in the file) and map mzbins to m/z:
            full_mz = mza["Full_mz_array"][:]
            mzbins = mza["Arrays_mzbin" + mzapath + "/" + str(scan)][:]
            mz_array = np.array([full_mz[i] for i in mzbins])
        else:
            mz_array = np.array(mza["Arrays_mz" + mzapath + "/" + str(scan)][:])
        intensities_array = mza["Arrays_intensity" + mzapath + "/" + str(scan)][:]
        intensities_array = intensities_array/1000 # scale intensity to avoid overflow 

        spectrum = {} # Create a dictionary with keys as unit mz and values as intensities
        for i in range(0, mz_array.size):
            mzx = int(np.floor(mz_array[i]))

            if mzx > maxMz:
                maxMz = mzx

            if mzx in spectrum:
                spectrum[mzx] += intensities_array[i]
            else:
                spectrum[mzx] = intensities_array[i]

        rtx = np.round(metadata["RetentionTime"][k], decimals=1)
        if rtx in lcms:
            for mzx in spectrum:
                if mzx in lcms[rtx]:
                    lcms[rtx][mzx] += spectrum[mzx]
                else:
                    lcms[rtx][mzx] = spectrum[mzx]
        else:
            lcms[rtx] = spectrum

        maxIntx = np.max(list(lcms[rtx].values()))
        if maxIntx > maxIntensity:
            maxIntensity = maxIntx
    mza.close()
    
    # Create image from 0 coordinate center (bottom left in image) to be comparable across runs:
    # Create a new RGB image with the same size as the input image
    width = int(np.max(list(lcms.keys())) * 10) + 1
    height = maxMz + 1

    img = Image.new('RGB', (width, height))

    # Get the pixel access object for the image
    pixels = img.load()

    # Define the dark blue to yellow color gradient
    colors = [(0, 0, 0.5), (1, 1, 0)]
    # Create a colormap with the defined gradient
    cmap = LinearSegmentedColormap.from_list('mycmap', colors)

    for i in lcms.keys(): # width or x, size[0]
        for j in lcms[i].keys(): # height or y, size[1]

            x = int(lcms[i][j])
            if x < (maxIntensity * (LcmsImageMinIntensityPercentage/100)):
                continue
            if x > (maxIntensity * (LcmsImageMaxIntensityCeilingPercentage/100)):
                x = maxIntensity * (LcmsImageMaxIntensityCeilingPercentage/100)
            # scale intensity value to 255 and max:
            x = int(np.log10(lcms[i][j]))
            x /= np.log10(maxIntensity * (LcmsImageMaxIntensityCeilingPercentage/100))
            x = np.array(cmap(x)[0:3]) # ignore the 4th channel (alpha)
            x *= 255
            # Set the pixel in the RGB image
            pixels[int(i*10), maxMz-j] = tuple(map(int,tuple(x)))
    
    #pixels = np.transpose(pixels)
    img.save(outputFullPath + ".jpg")
    img.close()
