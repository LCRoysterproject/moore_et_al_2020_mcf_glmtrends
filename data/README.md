This is the `data` folder structure outline. The `data` folder will include all data relevant to this specific repository project. This folder will include the **master dataset** that is unedited from the the `master_data` repository (in the `archive` folder). The dataset used for project analysis will preferably be a single data file in two formats (.csv, .xslx) located in this folder (`data`). This file can be named "data.csv". This dataset does not need additional naming parameters, unless it is subsetted.      

### Naming Standards  
Every dataset file is required to be in lowercase, no uppercase (camelcase) nor all caps, all names with separate words need to include a underscore ( _ ) with **no spaces**, no dates in the names unless it helps with the descriptions of the content (i.e., discharge_1941_2018.csv, landings_1986_2018.txt).  

Mandatory folders in `data` are:  
- archive  
  
 The `archive` folder consists of the **master dataset** and any raw datasets used to update the project repository dataset. 

Allowable optional folders in `data` could be:  
- water  
- heights  
- counts 
- transect  
- shell_biomass 
  
It is recommended only to add additional folders if there are more than 10 different datasets (this is **NOT** normal for one individual repository project). 
  
### General Data Do’s and Don’ts  
For all data collected, the methods, instruments, and protocol summaries used for data collection must be included in the README.md.    
  
### Color and format  
Excel files should **NOT** contain any color formatting, comments, formulas, graphs, signatures, or designs. These datasets will be returned to the submitter to clear the formatting. These features may not migrate into new file formats and may get lost.  
  
### Blanks  
Please do not leave cells blank. Blank cells leave a lot up for interpretation-- are the data missing? Was it not collected? Did something go wrong? Therefore, if there is no information available, use “N_A” for character values and “-999” for numerical values.  
  
### Headers  
All headers, variables, and columns must be defined, and units should be provided, no uppercase or spaces between words. Define any acronyms or abbreviations used in the dataset in the README.md. An outside user should be able to understand your dataset using this information.
