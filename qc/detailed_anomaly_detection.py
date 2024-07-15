import pandas as pd
from sklearn.ensemble import IsolationForest
from scipy.stats import zscore
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


 # rterrorAbsThreshold: IM DIA 0.2, Metab DDA 0.03. Prot DDA 0.5
def detect_outliers(data, mzerrorppmAbsThreshold = 15, rterrorAbsThreshold = 0.3, aterrorAbsThreshold = 0.1):

    # Group by MOLECULE if it exists and create columns for each molecule and metric
    if 'MOLECULE' in data.columns:
        data.drop(columns=["MZERROR"], inplace=True) # to keep only error in ppm
        # keep the first 4 columns plus all columns named with the substring "ERROR"
        error_columns = [col for col in data.columns if 'ERROR' in col]
        filtered_columns = data.columns[:4].tolist() + error_columns
        data = data[filtered_columns]
        # transform data frame to wide format, where each error column is named with a suffix of the MOLECULE value and the content is the corresponding error value
        data = data.pivot_table(index=['MSRUN', 'LABELSAMPLEGROUP', 'MSRUNID'],
                         columns='MOLECULE',
                         values=error_columns)
        # Flatten multi-level column index
        data.columns = [f"{col[0]}_{col[1]}" for col in data.columns]
        # Reset index to make MSRUN, LABELSAMPLEGROUP, and MSRUNID as columns
        data = data.reset_index()
    
    outliers_df = pd.DataFrame()
    for group_key, group in data.groupby('LABELSAMPLEGROUP'):

        # Perform outlier detection for each metric column 
        for column in data.columns[4:].tolist():
            # Select specified column for outlier detection
            X = group[[column]]
            X = X.dropna()
            if X.shape[0] < 3:
                continue  # Skip if there are less than 3 samples in the group
    
            # Adjust contamination based on variance
            if X.var().values[0] < 0.01:  # Adjust threshold as needed
                cont = 0.0001  # A very small value
            else:
                cont = 0.01  # A small value, but higher than for low-variance data

            clf = IsolationForest(contamination=cont, random_state=0)
            outliers = clf.fit_predict(X)
            # Filter dataframe to get rows with outliers
            outlier_indices = X.index[outliers == -1]

            # Apply thresholds to filter instacens with small errors:
            if 'MZERRORPPM' in column and len(outlier_indices) > 0:
                indexes = X[X[column].abs() > mzerrorppmAbsThreshold].index
                outlier_indices = list(set(indexes).intersection(outlier_indices))
            elif 'RTERROR' in column and len(outlier_indices) > 0:
                indexes = X[X[column].abs() > rterrorAbsThreshold].index
                outlier_indices = list(set(indexes).intersection(outlier_indices))
            elif 'ATERROR' in column and len(outlier_indices) > 0:
                indexes = X[X[column].abs() > aterrorAbsThreshold].index
                outlier_indices = list(set(indexes).intersection(outlier_indices))

            if len(outlier_indices) > 0:
                outlier_scores = clf.decision_function(X.loc[outlier_indices])
                # Create DataFrame for outlier rows
                outlier_rows = pd.DataFrame({
                    'LABELSAMPLEGROUP': [group_key] * len(outlier_indices),  # Store the sample group
                    'MSRUNID': group.loc[outlier_indices, 'MSRUNID'].values,
                    'Metric': [column] * len(outlier_indices),  # Store the name of the metric
                    'MetricValue': group.loc[outlier_indices, column].values,
                    'OutlierScore': outlier_scores,  # Store the outlier scores,
                    'Zscore': zscore(X[column])[outlier_indices] # Store Z-score
                })
                # Concatenate outlier rows with outliers_df
                outliers_df = pd.concat([outliers_df, outlier_rows])
    return outliers_df


def calculate_percentage_error(group):
    mean_abundance = group['ABUNDANCE'].mean()
    group['ABUNDANCEERROR'] = (abs(group['ABUNDANCE'] - mean_abundance) / mean_abundance) * 100
    return group


def detect_outsidetolerances(data, mzerrorppmAbsThreshold = 15, rterrorAbsThreshold = 0.3, aterrorAbsThreshold = 0.1, abundanceerrorAbsThreshold = 30):

    # Calculate abundance percentage error:
    data = data.groupby(['LABELSAMPLEGROUP', 'MOLECULE']).apply(calculate_percentage_error).reset_index(drop=True)

    # Group by MOLECULE if it exists and create columns for each molecule and metric
    if 'MOLECULE' in data.columns:
        data.drop(columns=["MZERROR"], inplace=True) # to keep only error in ppm
        # keep the first 4 columns plus all columns named with the substring "ERROR"
        error_columns = [col for col in data.columns if 'ERROR' in col]
        filtered_columns = data.columns[:4].tolist() + error_columns
        data = data[filtered_columns]
        # transform data frame to wide format, where each error column is named with a suffix of the MOLECULE value and the content is the corresponding error value
        data = data.pivot_table(index=['MSRUN', 'LABELSAMPLEGROUP', 'MSRUNID'],
                         columns='MOLECULE',
                         values=error_columns)
        # Flatten multi-level column index
        data.columns = [f"{col[0]}_{col[1]}" for col in data.columns]
        # Reset index to make MSRUN, LABELSAMPLEGROUP, and MSRUNID as columns
        data = data.reset_index()
    
    outsidetol_df = pd.DataFrame()
    for group_key, group in data.groupby('LABELSAMPLEGROUP'):
        # Perform outlier detection for each metric column 
        for column in data.columns[4:].tolist():
            # Select specified column for z score calculation
            X = group[[column]]
            X = X.dropna()
            if X.shape[0] < 3:
                continue 
            indexes = []
            # Apply thresholds to filter instacens with small errors:
            if 'MZERRORPPM' in column:
                indexes = X[X[column].abs() > mzerrorppmAbsThreshold].index
            elif 'RTERROR' in column:
                indexes = X[X[column].abs() > rterrorAbsThreshold].index
            elif 'ATERROR' in column :
                indexes = X[X[column].abs() > aterrorAbsThreshold].index
            elif 'ABUNDANCEERROR' in column :
                indexes = X[X[column].abs() > abundanceerrorAbsThreshold].index

            if len(indexes) > 0:
                new_rows = pd.DataFrame({
                    'LABELSAMPLEGROUP': [group_key] * len(indexes),  # Store the sample group
                    'MSRUNID': group.loc[indexes, 'MSRUNID'].values,
                    'Metric': [column] * len(indexes),  # Store the name of the metric
                    'MetricValue': group.loc[indexes, column].values,
                    'Zscore': zscore(X[column])[indexes] # Store Z-score
                })
                # Concatenate 
                outsidetol_df = pd.concat([outsidetol_df, new_rows])
    #outsidetol_df.dropna(inplace=True)
    return outsidetol_df


def plot_heatmap(data, outputPath):
    if len(data) < 1:
        return
    
    data['MSRUN'] = data['LABELSAMPLEGROUP'] + '_' + [str(x) for x in data['MSRUNID']]
    dfwide = data.pivot_table(index=['MSRUN'],
                            columns='Metric',
                            values='Zscore')
    dfannot = data.pivot_table(index=['MSRUN'],
                            columns='Metric',
                            values='MetricValue')
    
    vmin = min(-1, min(data['Zscore']))
    vmax = max(max(data['Zscore']), 1)
    num_rows, num_columns = data.shape

    # Set the minimum width for each column
    min_widths = [1] * num_columns  # List of minimum widths for each column
    total_width = sum(min_widths) # Calculate the total width of the heatmap
    # Set the aspect ratio based on the total width
    aspect_ratio = total_width / (num_rows)
    total_height = sum([0.2] * num_rows)
    total_height = max(total_height, 1) # to make height 1 inch minimum
    # Create the heatmap with adjusted aspect ratio
    plt.figure(figsize=(total_width, total_height))
    #plt.figure(figsize=(num_rows*0.5, num_columns*0.5))
    myfmt = ".2f"
    if np.nanmax(abs(dfannot.values)) > 9999.99:
        myfmt = ".2g" # set format to scientific notation if numbers are big:
    ax = sns.heatmap(dfwide, cmap='coolwarm', annot=dfannot, fmt=myfmt, cbar_kws={'label': 'Zscore'}, center=0, vmin=vmin, vmax=vmax, linewidths=0.5)
    ax.set(xlabel="", ylabel="")
    ax.xaxis.tick_top()
    for text in ax.texts: # Rotate annotations
        #text.set_rotation(90)
        text.set_fontsize(7)
        text.set_color('black')
    plt.xticks(rotation=90)
    #plt.yticks(rotation=-45)
    plt.subplots_adjust(right=0.8) # Adjust the space on the right to make room for legend
    #plt.show(block=False)
    plt.savefig(outputPath + ".pdf", format="pdf", bbox_inches="tight")


# import glob
# # # List of file paths
# filesMetrics = glob.glob("E:/IonToolPack/PeakQC/manuscript/Metab_IM-DIA-Agilent-HILIC-neg_mass-error/ResultsQC/Metrics*.csv")
# # filesMetrics = glob.glob("E:/IonToolPack/PeakQC/manuscript/Metab_DDA-lumos-HILIC-neg_LC-issues/ResultsQC/Metrics*.csv")
# # #filesMetrics = glob.glob("E:/IonToolPack/PeakQC/manuscript/**/Metrics*.csv", recursive=True)
# # Ignore the UserWarning
# import warnings
# warnings.filterwarnings("ignore", category=UserWarning)

# for f in filesMetrics:
#     if f.endswith('Outliers.csv'):
#         continue
#     print(f)
#     df = pd.read_csv(f)
#     outliers = detect_outliers(df)
#     plot_outlier_heatmap(outliers, f.replace(".csv", "_Outliers"))

#     if 'Spectra' not in f: # generate plots for outside tolerances
#         df = pd.read_csv(f)
#         df = df[df['LABELSAMPLEGROUP'].str.contains('qc', case=False)]
#         outsideTolerance = detect_outsidetolerance(df)
#         plot_outlier_heatmap(outsideTolerance, f.replace(".csv", "_Outside-tolerance"))


# # Reset warnings to their default behavior
# warnings.resetwarnings()
# input("Press the Enter key to continue: ") 
