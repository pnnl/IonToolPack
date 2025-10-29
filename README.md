
# IonToolPack

<div align="center">
<mark><strong>-> Download the latest IonToolPack.exe version from the <a href="https://github.com/pnnl/IonToolPack/releases">Releases</a>.</strong></mark>
</div>

IonToolPack is a software suite housing tools for mass spectrometry data. It reads data from multiple instrument formats, requires no installation and provides omics agnostic functionalities (metabolomics, lipidomics, proteomics, etc.) through a simplified and intuitive GUI.

**Available Features (New GUI with tabs per tool):**
- **Mirador**: Interactive data visualization including extracted ion chromatograms (XIC), extracted ion mobility (XIM) heatmaps, and MS/MS mirror plots with customizable m/z, RT, and arrival time ranges and tolerances
- **PeakQC**: Automated quality control pipeline with PCA analysis, outlier detection and comprehensive metrics extraction for either user specified ion targets or auto-tracked ions  
- **TandemMatch**: MS/MS spectral library matching with support for MSP and CSV library formats
- **PeakQuant**: Targeted MS1 peak abundance extraction for quantitation
- **CompareFeatures**: Cross-platform feature comparison tool for harmonizing and analyzing results (CSV files) from different acquisition methods or processing software


## Usage
1. Download the latest version (Release section, right panel) and decompress it
2. Double click IonToolPack.exe
3. Import raw MS files and click “Process”
4. See example input and output files in <a href="https://github.com/pnnl/IonToolPack/tree/master/test_data">test_data</a>.

## MS data supported
Supported formats include Agilent 'd', Thermo '.raw', Bruker 'd', and mzML, and for different types of MS acquisition methods:
* LC-MS
* LC-IMS-MS
* With/without fragmentation spectra in DDA or DIA mode
* Direct infusion 

## Contact
aivett.bilbao@pnnl.gov

## Contributors
* Aivett Bilbao, aivett.bilbao@pnnl.gov
* Andrea Harrison, andrea.harrison@pnnl.gov

## References

If you use this tool or any portions of this code please cite: 
* Harrison et al. "PeakQC: A Software Tool for Omics-Agnostic Automated Quality Control of Mass Spectrometry Data". Journal of the American Society for Mass Spectrometry 2024 https://doi.org/10.1021/jasms.4c00146.
* Bilbao et al. "MZA: A Data Conversion Tool to Facilitate Software Development and Artificial Intelligence Research in Multidimensional Mass Spectrometry". Journal of Proteome Research 2023 https://doi.org/10.1021/acs.jproteome.2c00313.
