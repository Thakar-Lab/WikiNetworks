# <font color='red'>WikiNetworks tutorial</font> 

This tutorial demonstrates the use of the WikiNetworks package to download and process pathway representations from WikiPathways database into machine-readable networks that can be used with visualization software or in downstream programming applications

*Author: Mukta Palshikar, mukta_palshikar@urmc.rochester.edu*

### <font color='red'>Install WikiNetworks</font>  

**WikiNetworks requires Python 3.6 or newer.**

With conda:
    
    conda wikinetworks create -f wikinetworks.yml
    conda activate wikinetworks
    conda install wikinetworks
    
With pip:

    pip install --upgrade wikinetworks

### <font color='red'>Import all functions from wikinetworks</font>   


```python
from wikinetworks import *
```


### <font color='red'> Create an instance of the WikiPathways database using the bioservices package</font>  


```python
__all__ = ["WikiPathways"]
s = WikiPathways()
```

### <font color='red'>Specify the pathways from WikiPathways that you wish to convert into networks. </font>  
Wikipathways entries are specified by Wikipathways ID - you can find this ID in the URL for the pathway of interest.

For example, the URL for the Apoptosis (*Homo sapiens*) pathway is https://www.wikipathways.org/index.php/Pathway:WP254 and the Wikipathways ID for this pathway is WP254. The ID begins with 'WP' and is followed by a string of digits.


Here, we select the Sandbox pathway test https://www.wikipathways.org/index.php/Pathway:WP4 to demonstrate WikiNetworks' functionality.

Note that the pathway ID is the first element in a list, and is *not* specified as a string. The user can extend this list to include other pathways to download and process.


```python
pathID = ["WP4"]
```

### <font color='red'>Download and view the curation tags for the pathway of interest</font>   

#### <font color='red'> Download a dictionary of curation tags for the selected pathway using the WikiPathways API </font>   


```python
curationTagsDict = getCurationTags(s, pathID) 
```

#### <font color='red'>Print pathway and curation information</font>

<font color = 'red'>*Caveat*: </font>We advise users to carefully examine the curation tags. Curation tags usually indicate whether a manual curator found issues in the pathway entry on WikiPathways, and that these errors have not yet been fixed by the pathway creator. While WikiNetworks attempts to automatically find and fix these common errors, we encourage the users of WikiNetworks to go over the output network and double-check interactions. Despite our best efforts, it's possible that the WikiNetworks algorithm does not find or incompletely fixes the pathway issues. This is due to the widely varying drawing styles and quality of the pathways in Wikipathways.


```python
processCurationTags(curationTagsDict)
```

### <font color='red'>Download and process pathways into networks</font> 

This is the driver function for the WikiNetworks processing algorithm. This function downloads GPML files for the requested pathway(s) from WikiPathways, processes the pathway into a network and attempts to correct drawing errors, and finally outputs the network as:
1. a NetworkX digraph object that can be used for downstream programmatic applications, 
1. a graphml file that can be imported into other simulation, processing or visualization tools such as CasQ, CellNOpt and Cytoscape, and
1. a simple interaction format (SIF) file that can be used similarly to the graphml file and is also easily human-readable.


```python
graph = runParsePathway(s, pathID)
```

#### <font color='red'>Print basic network statistics</font> 

Print the number of nodes and edges in the network.


```python
print("Nodes:", len(graph.nodes())) #print number of nodes in the pathway
print("Edges:", len(graph.edges())) #print number of edges in the pathway
```

#### <font color='red'>Print the names of the output files </font> 


```python
print("SIF graph:", '_'.join([pathID, "edgeList.csv"])) #location of the edgelist file. SIF = simple interaction format
print("Graphml: ", ''.join([pathID, "_graph.graphml"])) #location of the graphml file
```

### <font color='red'>Download and process multiple pathways into networks</font>

To download and process multiple pathways from WikiPathways, it is straightforward to specify a list of pathways of interest and process them in a loop, as shown below.



```python
testPathways = ["WP4", "WP2727","WP2456"] # Specify the test pathways using the WikiPathways IDs
for pathID in testPathways:
    curationTagsDict = getCurationTags(s, pathID) # download a dictionary of curation tags for the selected pathway
    processCurationTags(curationTagsDict) # print pathway and curation information
    graph = runParsePathway(s, pathID) # download and process pathways, generate a NetworkX digraph object
    print("Nodes:", len(graph.nodes())) # print number of nodes in the pathway
    print("Edges:", len(graph.edges())) # print number of edges in the pathway
    print("SIF graph:", '_'.join([pathID, "edgeList.csv"])) # location of the edgelist file. SIF = simple interaction format
    print("Graphml: ", ''.join([pathID, "_graph.graphml"])) # location of the graphml file
    print("**********\n")

```
