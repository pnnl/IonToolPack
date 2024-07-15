import os
import numpy as np
from PIL import Image
from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as pltcolors
import matplotlib.backends.backend_pdf

def PerformPCA(images, msRuns, labelGroups, runIds, outputFolder, mzaFiles, display = False, minIntensityPresencePercentage=80):

    # create an empty list to store the image data
    images1D = []

    minW = 1000000
    minH = minW
    # loop through all the images in the folder and check their sizes to keep the smallest dimension:
    for filename in images:
        if filename.endswith('.jpg'):
            # load the image and flatten its data
            img = Image.open(filename).convert('L')
            width, height = img.size
            if width < minW:
                minW = width
            if height < minH:
                minH = height
    
    # loop images to crop them, add their flattened data to the lists
    for filename in images:
        if filename.endswith('.jpg'):
            img = Image.open(filename).convert('L')
            width, height = img.size
            # box for cropping is a 4-tuple defining the left, upper, right, and lower pixel coordinate. 
            #   Origin is top left of image. MZ x RT origin is bottom left.
            img = img.crop((0, height - minH, minW, height))
            img = np.array(img).flatten()
            images1D.append(img)
            
    # convert the image data and labelGroups to numpy arrays
    images1D = np.array(images1D)
    labelGroups = np.array(labelGroups)
    runIds = np.array(runIds)

    # perform PCA on the image data
    pca = PCA(n_components=2)
    pca.fit(images1D)
    transformed_data = pca.transform(images1D)

    # Create data frame for saving formated output:
    df = pd.DataFrame({"MSRUN": msRuns, "LABELSAMPLEGROUP": labelGroups, "MSRUNID": runIds, "PC1": transformed_data[:, 0], "PC2": transformed_data[:, 1], "MZAPATH": mzaFiles})
    df["MSRUN"] = [os.path.basename(x).replace(".jpg", "") for x in df["MSRUN"]]
    df["PC1"] = round(df["PC1"], ndigits=3)
    df["PC2"] = round(df["PC2"], ndigits=3)
    df.to_csv(outputFolder + "/PCA.csv", index=False)

    myColors = list(pltcolors.TABLEAU_COLORS)
    while len(myColors) < len(set(df["LABELSAMPLEGROUP"])):
        myColors.extend(list(pltcolors.TABLEAU_COLORS)) # Add more colors until the number of groups if covered

    color_dict = {}
    # create a dictionary to map labelGroups to colors:
    for k,value in enumerate(set(df["LABELSAMPLEGROUP"])):
        color_dict[value] = myColors[k]

    # create an array of colors based on the labelGroups
    colors = [color_dict[label] for label in df["LABELSAMPLEGROUP"]]

    plt.close('all')
    # plot the transformed data using matplotlib and color the points based on the labelGroups
    #plt.scatter(df["PC1"], df["PC2"], c=colors, label=df["LABELSAMPLEGROUP"])
    # Iterate over unique label groups
    for label_group in df["LABELSAMPLEGROUP"].unique():
        group_df = df[df["LABELSAMPLEGROUP"] == label_group]
        plt.scatter(group_df["PC1"], group_df["PC2"], label=label_group)
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', fontsize=8)
    plt.subplots_adjust(right=0.8)  # Adjust the space on the right to make room for legend

    # Label runs:
    for i in range(df.shape[0]):
        plt.text(x=df["PC1"][i]+0.3, y=df["PC2"][i]+0.3, s=(str(df["MSRUNID"][i])), 
            rotation=40, rotation_mode="anchor",     
            fontdict=dict(color="black",size=8))
    plt.xlabel("PC1, explained variance " + "{:.3f}".format(pca.explained_variance_ratio_[0]))
    plt.ylabel("PC2, explained variance " + "{:.3f}".format(pca.explained_variance_ratio_[1]))

    #plt.savefig(outputFolder + "/PCA.jpg")
    plt.savefig(outputFolder + "/PCA.pdf")
    if display:
        plt.show(block=False)

    # Find common ions per LABELSAMPLEGROUP: --------------
    images1DFreq = np.copy(images1D)
    # normalize pixel values to 0 and 1:
    for k in range(images1DFreq.shape[0]):
            threshold = np.max(images1DFreq[k]) * (minIntensityPresencePercentage/100)
            images1DFreq[k][images1DFreq[k] < threshold] = 0
            images1DFreq[k][images1DFreq[k] > 0] = 1

    # build a data frame with detected ions and frequencies:
    dfIons = pd.DataFrame()
    df["IndexOriginal"] = range(1, df.shape[0]+1)
    for label in df['LABELSAMPLEGROUP'].unique():
        if "blank" in label.lower(): # skip blanks, no needed to check RT and m/z shift
            continue
        runIndexes = np.array(df[df['LABELSAMPLEGROUP'] == label]["IndexOriginal"]) - 1
        dfimg = pd.DataFrame(np.transpose(images1DFreq[runIndexes]))
        dfimg["FREQ"] = dfimg.sum(axis=1)
        dfimg["indexPixel1D"] = range(0, dfimg.shape[0])
        dfimg["MZ"] = minH - np.floor(dfimg["indexPixel1D"] / minW) - 1 # Substract minH to calculate m/z: bottom left in image, but 1D array starts at top left (pixel flattened image 1D)
        dfimg["RT"] = dfimg["indexPixel1D"] % minW
        # Take the INTENSITY from the first run in group for this label:
        dfimg["INTENSITY"] = np.transpose(images1D[runIndexes[0]])
        dfimg = dfimg[(dfimg["FREQ"] > 0) & (dfimg["INTENSITY"] > 0)]
        dfimg["LABELSAMPLEGROUP"] = label
        dfIons = pd.concat([dfIons, dfimg[["LABELSAMPLEGROUP", "MZ", "RT", "FREQ", "INTENSITY"]]])

    return dfIons
        #dfimg.to_csv(outputFolder + "/Ions-" + label + ".csv", index=False)

