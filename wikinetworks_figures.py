#Thie file reproduces all figures in the manuscript WikiNetworks: translating manually created biological pathways for topological analysis
# Author: Mukta Palshikar, mukta_palshikar@Urmc.rochester.edu


# Import all functions from wikinetworks
from wikinetworks import *

if __name__ == '__main__':
    __all__ = ["WikiPathways"]
    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.max_rows', 10000)

    s=WikiPathways()
    testPathways = ['WP4482', 'WP3935', 'WP4868', 'WP49', 'WP127', 'WP195', 'WP205', 'WP231', 'WP286', 'WP23'] #['WP254', 'WP23', "WP4495", "WP437",  "WP4495", "WP195", "WP49", "WP286", "WP395", "WP127", "WP364", "WP205", "WP231", "WP366"] # "WP1835", "WP1840",  "WP1919", "WP1836",
    for pathID in testPathways:
        graph = runParsePathway(s, pathID)
        print("Nodes:", len(graph.nodes()))
        print("Edges:", len(graph.edges()))
        print("SIF graph:", '_'.join([pathID, "edgeList.csv"]))
        print("Graphml: ", ''.join([pathID, "_graph.graphml"]))
        print("***")
    
    #unitTest("WP4482_112969.gpml", "WP4482")
    #testToyPathways()

    evaluateProgramOutput()
    evaluateCytoscapeOutput()
    #evaluateCytoscapeOutput_fairplay()    
    makeComparisonPlots()
    #compareCytoscapeWithProgramOutput()
