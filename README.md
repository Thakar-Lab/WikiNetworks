![PyPI - Python Version](https://img.shields.io/pypi/pyversions/wikinetworks?style=for-the-badge)

# WikiNetworks
WikiNetworks is a user-friendly Python package to process pathways from the WikiPathways database into machine-readable network representations.
## Summary
Summary: WikiPathways is a database of 2979 biological pathways across 31 species created using the drawing software PathVisio. These pathways are not directly usable for network-based topological analyses due to differences in curation styles and drawings. We developed the WikiNetworks package to standardize and construct directed networks from these pathway representations. 
The WikiNetworks software reconstructs networks from these pathways by combining geometric information and manual annotations from the hand-drawn pathways. WikiNetworks performs significantly better than existing tools. This enables the use of high-quality WikiPathways resources for network-based topological analysis of high-throughput data.  WikiNetworks allows users to (a) download manually created gpml files from WikiPathways by specifying pathway id and (b) convert them to machine-interpretable networks suitable for use with pathway/network analysis tools and (c) for import into pathway visualization tools.  
## Availability and Implementation: 
WikiNetworks is written entirely in Python and is available on github.com/Thakar-Lab and on pip. Tested for use with Python 3.7 in Windows and Linux environments.
## Citation:
Palshikar MG, Hilchey SP, Zand MS, Thakar J. WikiNetworks: translating manually created biological pathways for topological analysis. Bioinformatics. 2021 Oct 12:btab699. doi: 10.1093/bioinformatics/btab699. Epub ahead of print. PMID: [34636843](https://pubmed.ncbi.nlm.nih.gov/34636843/) 

## Installation Instructions

**WikiNetworks requires Python 3.6 or newer.**

With pip:

    pip install --upgrade wikinetworks

## Contact: 
Juilee_Thakar@URMC.Rochester.edu; Mukta_Palshikar@URMC.Rochester.edu
## Bug Tracker:
https://github.com/Thakar-Lab/WikiNetworks/issues
## Tutorials:
https://github.com/Thakar-Lab/WikiNetworks/tree/main/Tutorials
## Basic Usage:
### Import all functions from wikinetworks


```python
from wikinetworks import *
```

### Create an instance of the WikiPathways database using the bioservices package


```python
__all__ = ["WikiPathways"]
s = WikiPathways()
```

    WARNING [bioservices:WikiPathways]:  URL of the services contains a double //.Check your URL and remove trailing /


### Specify the pathways from WikiPathways that you wish to convert into networks.
Wikipathways entries are specified by Wikipathways ID - you can find this ID in the URL for the pathway of interest.

Here, we select the Sandbox pathway test https://www.wikipathways.org/index.php/Pathway:WP4 to demonstrate WikiNetworks' functionality.

Note that the pathway ID is the first element in a list, and is *not* specified as a string. The user can extend this list to include other pathways to download and process.


```python
pathID = "WP4"
```

### Download and process pathways into networks

This is the driver function for the WikiNetworks processing algorithm. This function downloads GPML files for the requested pathway(s) from WikiPathways, processes the pathway into a network and attempts to correct drawing errors, and finally outputs the network as:
1. a NetworkX digraph object that can be used for downstream programmatic applications, 
1. a graphml file that can be imported into other simulation, processing or visualization tools such as CasQ, CellNOpt and Cytoscape, and
1. a simple interaction format (SIF) file that can be used similarly to the graphml file and is also easily human-readable.


```python
graph = runParsePathway(s, pathID)
```
