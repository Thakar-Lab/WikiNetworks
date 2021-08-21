# WikiNetworks
WikiNetworks is a user-friendly Python package to process pathways from the WikiPathways database into machine-readable network representations.
## Summary
Summary: WikiPathways is a database of 2979 biological pathways across 31 species created using the drawing software PathVisio. These pathways are not directly usable for network-based topological analyses due to differences in curation styles and drawings. We developed the WikiNetworks package to standardize and construct directed networks from these pathway representations. 
The WikiNetworks software reconstructs networks from these pathways by combining geometric information and manual annotations from the hand-drawn pathways. WikiNetworks performs significantly better than existing tools. This enables the use of high-quality WikiPathways resources for network-based topological analysis of high-throughput data.  WikiNetworks allows users to (a) download manually created gpml files from WikiPathways by specifying pathway id and (b) convert them to machine-interpretable networks suitable for use with pathway/network analysis tools and (c) for import into pathway visualization tools.  
## Availability and Implementation: 
WikiNetworks is written entirely in Python and is available on github.com/Thakar-Lab and on pip. Tested for use with Python 3.7 in Windows and Linux environments. 
## Contact: 
Juilee_Thakar@URMC.Rochester.edu; Mukta_Palshikar@URMC.Rochester.edu
## Bug Tracker:
https://github.com/Thakar-Lab/WikiNetworks/issues
## Tutorials:
https://github.com/Thakar-Lab/WikiNetworks/tree/main/Tutorials