import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mza.mza import CreateMZA, Extract2DIonIntensityFrame

def GenerateXISurfacePlot(mzaFile, outputFilename, molecule, precMz, rt, at, mzHalfWindowXIC=0.01, rtrange=0.3, atrange=1.5):
    # Check if ion mobility data:
    with CreateMZA(mzaFile) as mza:
        metadata = mza["Metadata"]
        metadata = metadata[(metadata["MSLevel"] == 1)]
        if len(metadata["IonMobilityBin"] > 0) == 0:
            return # data has no ion mobility separation

        df = Extract2DIonIntensityFrame(mza, precMz, msLevel=1, startRT = rt-rtrange, endRT = rt+rtrange, startAT = at-atrange, endAT = at+atrange, mztolhalfwidth=mzHalfWindowXIC)
        df = df[df["intensity"] >= 1]
        # TODO: need to adjust figure size (w x h) and marker size (s) based on sampling frequency
        figWidth = 5 #len(np.unique(df["rtbin"]))/10 # 6.2 
        figHeight = 5 #len(np.unique(df["atbin"]))/10 #6
        #print(molecule)
        #print("w h = " + str(figWidth) + " " + str(figHeight))
        fig, ax = plt.subplots()
        ax.patch.set_facecolor('black')  # Set background color to black
        ax.set_xlim(rt-rtrange, rt+rtrange)
        ax.set_ylim(at-atrange, at+atrange)
        fig.set_figwidth(figWidth)
        fig.set_figheight(figHeight) # 12 for proteomics and 6 for metabolomics
        scatter = ax.scatter(x=df["rt"], 
                    y=df["at"], 
                    c=np.log10(df["intensity"]), 
                    cmap='viridis',
                    marker ='s',
                    s=30,
                    #linewidths = 0.9,
                    edgecolors = 'face')
        
        ax.set_xlabel("Retention time")
        ax.set_ylabel("Arrival time")
        # Add color bar
        cbar = plt.colorbar(scatter)
        cbar.set_label('Log10(Intensity)')
        fig.suptitle(f"Extracted ion surface for {molecule}")
        plt.savefig(outputFilename + ".jpg")
        plt.savefig(outputFilename + ".pdf")
        plt.close(fig)
