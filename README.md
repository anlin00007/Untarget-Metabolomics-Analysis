Untarget-Metabolomics-Analysis
==============================

Environment: Windows/Linux
Software: Python, R, MySQL
Package: Python (MySQLdb, csv), R (XCMS, CAMERA, Multtest)

Useage: 1. Import hmdbchem1.sql to local mysql database
        2. Set up raw data path in Metaanalysis.R
        3. Set up parameters for analysis
        4. Run metaanalysis.R
        5. Filter peak results by fold change and p-value
        6. Save filtered results in a new csv file
        7. Run annotation.py with the filtered results
        8. Retrieve annotation data from Metlin database
        9. Integrate annotation data
