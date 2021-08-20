# WikiNetworks
WikiNetworks is a user-friendly Python package to process pathways from the WikiPathways database into machine-readable network representations.



## Summary
WikiPathways is a database of 2979 manually created biologically pathways across 31 species. These pathways are primarily visual references and are not directly usable for network-based analyses. The WikiNetworks software reconstructs networks from these pathways by combining geometric information and manual annotations from the hand-drawn pathways. WikiNetworks allows users to (a) download manually created gpml files from WikiPathways by specifying pathway id and (b) convert them to machine-interpretable networks suitable for use with pathway/network analysis tools and (c) for import into pathway visualization tools.  
## Availability and Implementation: 
WikiNetworks is written entirely in Python and is available on github.com/Thakar-Lab and on pip. Tested for use with Python 3.7 in Windows and Linux environments. 
## Contact: 
Juilee_Thakar@URMC.Rochester.edu; Mukta_Palshikar@URMC.Rochester.edu
## Dependencies (Python packages):
 - networkx
 - re
 - urllib
 - csv
 - itertools
 - sys
 - bs4
 - random
 - requests
 - binascii
 - bioservices
 - numpy
 - shapely
 - scipy
 - pandas
 - fiona
 - geopandas
 - matplotlib
 - seaborn
 - glob
 - string
