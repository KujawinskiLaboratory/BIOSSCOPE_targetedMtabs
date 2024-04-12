# BIOSSCOPE_targetedMtabs
12 April 2024\
Krista Longnecker\
Woods Hole Oceanographic Institution

The most important file here is the Excel file labeled  ```KujawinskiWHOI_targetedMetabolites.2024.04.12.xlsx```. The details are as follows:
* The data has New_ID as the leftmost column, and will have the same number of rows as the discrete data file that is available on the BIOS-SCOPE Google Drive as of January 5, 2024
* The remaining columns are the metabolites, one column per metabolite
* A matrix-matched standard curve was used to convert from peak areas to concentration
* The units for metabolites are pM
*	The data are corrected for extraction efficiency, with the following exceptions
    * in the 2016 to 2019 there is no extraction efficiency information for:
        * 2-deoxyguanosine
        * hypoxanthine
    * 2021 data has no extraction efficiency for these metabolites:
        * 2-deoxyguanosine
        * hydroxypyruvic acid
        * hypoxanthine
        * myo-inositol
    * 2023 has no extraction efficiency for these metabolites:
        * cystine
        * glycine
        * hydroxocobalamin
        * hypoxanthine
        * lumichrome
        * nicotinamide
        * pantoic acid
        * beta-alanine (isom. alanine)
        * sarcosine
        * trigonelline
        * uric acid
     
The 2016 to 2019 data are already in the MetaboLights repository as [MTBLS2356](https://www.ebi.ac.uk/metabolights/editor/MTBLS2356/descriptors); these data are covered in the Longnecker et al. (2024) publication available [here](https://dx.doi.org/10.1002/lno.12497).

The data files here have metabolites that were detected during this project. There is a longer list of metabolites that were measured. Also note that these data are not corrected for limit of detection (LOD)/limit of quantification (LOQ). The data were collected on two different mass spectrometers, and the original TSQ showed changes in LOD over time. I decided to be as generous as possible and leave all possible values in here.

### mfiles folder
There are three m-files in here because this part of the computer processing was done in stages.
* ```matchKujTargetedMetabolites_2016to2019.m``` --> used to convert the samples processed up to July 2019 into a format that matches the BIOS-SCOPE discrete file
* ```matchKujTargetedMetabolites_2021.m``` --> used to process the samples from July 2019 up through 2021; these are the final TSQ samples
* ```matchKujTargetedMetabolites_2023.m``` --> used to process the samples collected up through July 2023; these are the first Altis samples

### dataFiles folder
This folder holds the MATLAB and Excel files that are either needed to merge in new targeted mtabolomics data, or are the output from prior merges.

### Altis_processingFiles folder
These are the files used by Krista to read in and process El-MAVEN files\
You should just need to edit ```riMAVEN.m``` to process one batch of samples <i>after</i> the El-MAVEN processing has been done. Hence this step of the analysis requires knowledge about the individual metabolites, how to use El-MAVEN, and how to use MATLAB.
