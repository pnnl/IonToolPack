import os
import pandas as pd
import multiprocessing
from multiprocessing.pool import Pool
import math
from qc.utils import FormatDataframeSamples
from qc.image_time_vs_mz import GenerateImageTimeVsMz
from qc.pca import PerformPCA
from qc.auto_ion_tracking import DetectTopmostIons
from qc.ion_batch import GenerateImageIonBatch
from qc.ms2 import GenerateMS2plot
from qc.xis import GenerateXISurfacePlot
from qc.spectra_metrics import ExtractSpectraMetadataMetrics
from qc.detailed_anomaly_detection import detect_outliers, plot_heatmap, detect_outsidetolerances
import string
import subprocess
import traceback 
import warnings
import glob
import time
import tomllib

def mza_conversion(myArgs, execPath):
    subprocess.run(execPath + myArgs, shell=True)

def qc_pipeline(dfruns, outputPath, configFile=""):
    if configFile == "":
        configFile = "config.toml" # get default config
    with open(configFile, "rb") as f:
        config = tomllib.load(f)

    mzaPath = os.path.normpath(os.path.join(outputPath, "DataMza"))
    resultsPath = os.path.normpath(os.path.join(outputPath, "ResultsQC"))

    # Record the start time
    start_time = time.time()

    if dfruns is None or len(dfruns) == 0:
        print("Error: No MS runs selected!")
        return
    if outputPath is None or outputPath == "" or not os.path.exists(outputPath):
        print("Error: Invalid output folder!")
        return
    
    dfruns = FormatDataframeSamples(dfruns)
                    
    # check and create folders if it does not exist
    if not os.path.exists(resultsPath):
        os.makedirs(resultsPath)

    # Use parallel processing
    nProcesses = math.floor(multiprocessing.cpu_count() * 0.6) # Use n = 70% of computer cores
    if nProcesses < 1:
        nProcesses = 1

    # ---------------------------------------------------------------
    # 1) Iterate list of LC-MS runs and check if mza format exist, otherwise convert it:
    print("Converting raw files to mza format...")

    minIntensityMza = config["MinIntensityMza"]
    dfruns["MZAPATH"] = ""
    myFiles = []
    useMultithreading = True
    for i, row in dfruns.iterrows():
        if row["MSRUNFORMAT"] != ".mza":
            if not os.path.exists(os.path.join(mzaPath, row["MSRUN"] + ".mza")):
                # mza file does not exist in project path, convert raw MS file to mza
                xpath = os.path.join(row["MSRUNPATH"], row["MSRUN"] + row["MSRUNFORMAT"])
                myFiles.append(' -file "' + xpath + '" -out "' + mzaPath + '" -intensityThreshold ' + str(minIntensityMza))
                if row["MSRUNFORMAT"] == ".d" and os.path.exists(os.path.join(xpath, "AcqData")): # don't use multithreading here for Agilent .d because it will be used per file by mza.exe
                    useMultithreading = False
            dfruns.loc[i,"MZAPATH"] = os.path.join(mzaPath, row["MSRUN"])        
        else:
            # mza file exists in initial path provided
            dfruns.loc[i,"MZAPATH"] = os.path.join(row["MSRUNPATH"], row["MSRUN"])

    if len(myFiles) > 0 and not os.path.exists(mzaPath):
        os.makedirs(mzaPath)
    execPath = '"' + os.path.join(os.getcwd(), 'mza', '"mza.exe')
    if useMultithreading and len(myFiles) > 0:
        nProcesses = min(nProcesses, len(myFiles))
        # create and configure the process pool
        with Pool(nProcesses) as pool:
            # issue tasks to the process pool
            for i in range(0,len(myFiles)):
                pool.apply_async(mza_conversion, args=(myFiles[i], execPath))
            pool.close() # close the process pool
            pool.join() # wait for all tasks to complete
    else:
        # Without parallel processing:
        for i in range(0,len(myFiles)):
           mza_conversion(myFiles[i], execPath)

    # ---------------------------------------------------------------
    # 2) Image time-vs-mz: Generate an image for each MS run (only most intense peaks)
    print("Generating time-vs-m/z images...")

    # check and create folder if it does not exist
    if not os.path.exists(os.path.join(resultsPath, "images-time-vs-mz")):
        os.makedirs(os.path.join(resultsPath, "images-time-vs-mz"))

    # create a list to generate images only if not found
    myFiles = []
    myOutputs = []
    for row in dfruns.itertuples():
        msfile = os.path.join(resultsPath, "images-time-vs-mz", 
                                str(row.LABELSAMPLEGROUP) + "_" + str(row.MSRUNID) + "_" + row.MSRUN)
        if not os.path.exists(msfile + ".jpg"):
            myFiles.append(row.MZAPATH + '.mza')
            myOutputs.append(msfile)
                
    # create and configure the process pool
    if len(myFiles) > 0:
        nProcesses = min(nProcesses, len(myFiles))
        with Pool(nProcesses) as pool:
            # issue tasks to the process pool
            for i in range(0,len(myFiles)):
                pool.apply_async(GenerateImageTimeVsMz, args=(myFiles[i], myOutputs[i], config["TimeVsMzImageMinIntensityPercentage"], config["TimeVsMzImageMaxIntensityCeilingPercentage"]))
            pool.close() # close the process pool
            pool.join() # wait for all tasks to complete

    # ---------------------------------------------------------------
    # 3) PCA and common ions: Perform PCA based on the LC-MS images and detect common ions
    print("Performing PCA analysis...")

    if not os.path.exists(os.path.join(resultsPath, "PCA.csv")) or not os.path.exists(os.path.join(resultsPath, "AutoTracked-Ions.csv")):
        myFiles = []
        myMzaPaths = []
        myRuns = []
        myGroups = []
        myRunIds = []
        for row in dfruns.itertuples():
            outputfile = os.path.join(resultsPath, "images-time-vs-mz", 
                                      str(row.LABELSAMPLEGROUP) + "_" + str(row.MSRUNID) + "_" + row.MSRUN + ".jpg")
            if os.path.exists(outputfile):
                myFiles.append(outputfile)
                myRuns.append(row.MSRUN)
                myGroups.append(row.LABELSAMPLEGROUP)
                myRunIds.append(row.MSRUNID)
                myMzaPaths.append(row.MZAPATH)
        if(len(myFiles) == 0):
            print("Error: No time-vs-mz images found!")
            return

        dfIons = PerformPCA(myFiles, myRuns, myGroups, myRunIds, outputFolder = resultsPath, mzaFiles = myMzaPaths, display = True, minIntensityPresencePercentage=config["MinIntensityPresencePercentage"])

        dfautoIons = DetectTopmostIons(dfIons, dfruns, config["AutoTrackedIonsTopN"], config["MinMzDistDetectCentroidMS"])
        pd.DataFrame.to_csv(dfautoIons, os.path.join(resultsPath, "AutoTracked-Ions.csv"), index=False)

    # ---------------------------------------------------------------
    # 4) SpectraMetrics: Generate a data frame with metrics of spectra for each mza file
    print("Extracting metrics spectra summary statistics...")
    spectraMetricsFile = os.path.join(resultsPath, "Metrics_Spectra.csv")
    nProcesses = min(nProcesses, len(dfruns))
    result = None
    with Pool(nProcesses) as pool:
        result = pool.starmap(ExtractSpectraMetadataMetrics, zip(dfruns["MZAPATH"] + '.mza'))

    result = pd.concat(result, ignore_index=True)
    result = pd.concat([dfruns[["MSRUN", "LABELSAMPLEGROUP", "MSRUNID"]],result], axis=1)
    result.to_csv(spectraMetricsFile, index=False)

    # ---------------------------------------------------------------    
    # 5) ImageIonBatch: Generate for each sample group (user specified: Blank, QC, Other)
    ionsFiles = []
    ionsMetricsFiles = []
    messages = []
    outputFolders = []
    userIonsMS2OutputFolder = "user-ions-ms2"
    userIonsXISOutputFolder = "user-ions-xis"
    if not os.path.exists(os.path.join(resultsPath, "Metrics_AutoTracked-Ions.csv")):
        ionsFiles.append(os.path.join(resultsPath, "AutoTracked-Ions.csv"))
        ionsMetricsFiles.append(os.path.join(resultsPath, "Metrics_AutoTracked-Ions"))
        messages.append("Generating overlaid ion images for auto tracked ions...")
        outputFolders.append("overlaid-images-ions")
    if not os.path.exists(os.path.join(resultsPath, "Metrics_User-Ions.csv")):
        ionsFiles.append(os.path.join(outputPath, "User-Ions.csv"))
        ionsMetricsFiles.append(os.path.join(resultsPath, "Metrics_User-Ions"))
        messages.append("Generating overlaid ion images for user-specified ions (theoretical or reference values)...")
        outputFolders.append("overlaid-images-user-ions")
    
    dfruns["legend"] = dfruns["LABELSAMPLEGROUP"].astype(str) + "_" + dfruns["MSRUNID"].astype(str)
    for i in range(0, len(ionsFiles)):
        try:
            ionsfile = ionsFiles[i]
            if not os.path.exists(ionsfile):
                continue

            print(messages[i])
            dfions = pd.read_csv(ionsfile)
            dfions.columns = dfions.columns.str.upper()

            # check and create folder if it does not exist
            outputfolder = os.path.join(resultsPath, outputFolders[i])
            if not os.path.exists(outputfolder):
                os.makedirs(outputfolder)
            userIonsMS2 = False
            if "User-Ions.csv" in ionsfile and "FRAGSMZ" in dfions.columns:
                userIonsMS2 = True
                userIonsMS2OutputFolder = os.path.join(resultsPath, userIonsMS2OutputFolder)
                os.makedirs(userIonsMS2OutputFolder)
            userIonsXIS = False # Generate Extracted Ion Surface images if both RT and AT are provided
            if "User-Ions.csv" in ionsfile and "RT" in dfions.columns and "AT" in dfions.columns:
                userIonsXIS = True
                userIonsXISOutputFolder = os.path.join(resultsPath, userIonsXISOutputFolder)
                os.makedirs(userIonsXISOutputFolder)
        
            # check columns ions file:
            if "MZVIEWHALFWINDOW" not in dfions.columns:
                dfions['MZVIEWHALFWINDOW'] = config["MZVIEWHALFWINDOW"]
            if "RTVIEWHALFWINDOW" not in dfions.columns:
                dfions['RTVIEWHALFWINDOW'] = config["RTVIEWHALFWINDOW"]
            if "ATVIEWHALFWINDOW" not in dfions.columns:
                dfions['ATVIEWHALFWINDOW'] = config["ATVIEWHALFWINDOW"]
            if "MZXICHALFWINDOW" not in dfions.columns:
                dfions['MZXICHALFWINDOW'] = config["MZXICHALFWINDOW"]
            # TODO: don't use inf in limit by default
            #dfions['RTVIEWHALFWINDOW'] = dfions['RTVIEWHALFWINDOW'].fillna(float('inf')) 
            #dfions['RTVIEWHALFWINDOW'] = dfions['RTVIEWHALFWINDOW'].replace("all", float('inf'))
            
            if "AT" not in dfions.columns: # AT = arrival time
                dfions["AT"] = 0

            dferrors = pd.DataFrame()
            for k in range(0, dfions.shape[0]):
                ionmz = dfions["MZ"][k]
                ionrt = dfions["RT"][k]
                ionat = dfions["AT"][k]
                molecule = dfions["MOLECULE"][k]
                mzViewHalfWindow = dfions["MZVIEWHALFWINDOW"][k]
                rtViewHalfWindow = dfions["RTVIEWHALFWINDOW"][k]
                mzHalfWindowXIC = dfions["MZXICHALFWINDOW"][k]
                atViewHalfWindow = dfions["ATVIEWHALFWINDOW"][k]

                # format MOLECULE name as valid file name
                valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
                molecule = ''.join(c for c in molecule if c in valid_chars)
                molecule = molecule.replace(' ','-')  
                suffixImage = molecule + "-MZ" + str(round(ionmz, ndigits=2)) + "-RT" + str(round(ionrt, ndigits=1))
                if ionat > 0:
                    suffixImage = suffixImage + "-AT" + str(round(ionat, ndigits=1))

                print("     analyzing: " + suffixImage)
                dfx = GenerateImageIonBatch(dfruns,
                                outputFolder=outputfolder,
                                mz=ionmz, 
                                rt=ionrt,
                                molecule=molecule,
                                suffixImage=suffixImage,
                                mzHalfWindowXIC=mzHalfWindowXIC,
                                rtrange=rtViewHalfWindow,
                                mzrange=mzViewHalfWindow, 
                                at=ionat,
                                atrange=atViewHalfWindow)
                if len(dfx) > 0:
                    dferrors = pd.concat([dferrors, dfx], ignore_index=True)

                # 6) ImageIonBatch MS/MS: Generate for each mza file 
                if userIonsMS2:
                    fragsMz = [float(value) for value in dfions["FRAGSMZ"][k].split(';')]
                    if "FRAGSINTENSITY" in dfions.columns:
                        fragsIntensity = [float(value) for value in dfions["FRAGSINTENSITY"][k].split(';')]
                    else:
                        fragsIntensity = [1 for value in fragsMz]
                    for j in range(len(dfruns)):
                        mzaFile = dfruns["MZAPATH"][j] + ".mza"
                        outputFilename = os.path.join(userIonsMS2OutputFolder, molecule + "-" + dfruns["legend"][j])
                        GenerateMS2plot(mzaFile, 
                                        outputFilename, 
                                        molecule + "-" + dfruns["legend"][j], 
                                        ionmz, 
                                        mzHalfWindowXIC, 
                                        ionrt, 
                                        fragsMz, 
                                        fragsIntensity, 
                                        ionat,
                                        config["MinMzDistDetectCentroidMS"])
            
                # 7) Extracted Ion Surface (XIS): Generate XIS for each mza file 
                if userIonsXIS:
                    nProcesses = min(nProcesses, len(dfruns))
                    # create and configure the process pool
                    with Pool(nProcesses) as pool:
                        # issue tasks to the process pool
                        for j in range(0,len(dfruns)):
                            mzaFile = dfruns["MZAPATH"][j] + ".mza"
                            outputFilename = os.path.join(userIonsXISOutputFolder, molecule + "-" + dfruns["legend"][j])
                            pool.apply_async(GenerateXISurfacePlot, args=(mzaFile, 
                                                    outputFilename, 
                                                    molecule + "-" + dfruns["legend"][j], 
                                                    ionmz, 
                                                    ionrt, 
                                                    ionat, 
                                                    mzHalfWindowXIC, 
                                                    max(rtViewHalfWindow,1),
                                                    max(atViewHalfWindow,3)))
                        pool.close() # close the process pool
                        pool.join() # wait for all tasks to complete
            if len(dferrors) > 0:
                pd.DataFrame.to_csv(dferrors, ionsMetricsFiles[i] + ".csv", index=False)
            else:
                print("     no ions found.")
        except:
            # printing stack trace 
            traceback.print_exc()
    
    # 8) Detailed anomaly detection:
    filesMetrics = glob.glob(resultsPath + "/Metrics*Ions.csv")
    if len(filesMetrics) > 0:
        print("Detailed anomaly detection: outliers and QC ions outside tolerances...")
        warnings.filterwarnings("ignore", category=UserWarning) # Ignore the UserWarning
        for f in filesMetrics:
            df = pd.read_csv(f)
            df = df[~df['LABELSAMPLEGROUP'].str.contains('blank', case=False)] # ignore blanks
            outliers = detect_outliers(df, 
                                       mzerrorppmAbsThreshold=config["MZERRORPPM"], 
                                       rterrorAbsThreshold=config["RTERROR"], 
                                       aterrorAbsThreshold=config["ATERROR"])
            if len(outliers) > 0:
                plot_heatmap(outliers, f.replace(".csv", "_Outliers"))
                pd.DataFrame.to_csv(outliers, f.replace(".csv", "_Outliers.csv"), index=False)

            # generate plots for outside tolerances, only QCs
            df = pd.read_csv(f)
            df = df[~df['LABELSAMPLEGROUP'].str.contains('blank', case=False)] # ignore blanks
            df = df[df['LABELSAMPLEGROUP'].str.contains('qc', case=False)]
            outsideTolerance = detect_outsidetolerances(df, 
                                    mzerrorppmAbsThreshold=config["MZERRORPPM"], 
                                    rterrorAbsThreshold=config["RTERROR"], 
                                    aterrorAbsThreshold=config["ATERROR"],
                                    abundanceerrorAbsThreshold=config["ABUNDANCEERROR"])
            if len(outsideTolerance) > 0:
                plot_heatmap(outsideTolerance, f.replace(".csv", "_Outside-tolerances"))
                pd.DataFrame.to_csv(outsideTolerance, f.replace(".csv", "_Outside-tolerances.csv"), index=False)

        warnings.resetwarnings() # Reset warnings to their default behavior
    
    print("Done!")
    end_time = time.time()
    # Calculate the total running time in minutes
    total_running_time = (end_time - start_time)/60
    strTime = f"Total Execution Time: {total_running_time:.2f} minutes"
    strTime += "\nNumber of MS runs: " + str(len(dfruns))
    print(strTime)
    with open(os.path.join(resultsPath, 'xlog_' + time.strftime("%Y-%m-%d-%H-%M-%S", time.gmtime(end_time)) + '.txt'),"w") as txtfile:
        txtfile.write(strTime)

