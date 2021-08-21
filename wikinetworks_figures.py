#Thie file reproduces all figures in the manuscript WikiNetworks: translating manually created biological pathways for topological analysis
# Author: Mukta Palshikar, mukta_palshikar@Urmc.rochester.edu


# Import all functions from wikinetworks
from wikinetworks import *

# Define functions required for plotting/analysis

def makeMatchedDistHist(featureDFs, pathwayID):
    df = featureDFs['interactDF']
    #print(df.head(10))
    plotName = ''.join([pathwayID, "_disthist.svg"])
    #sns.distplot(df['matched_Node1_dist'], hist = True, kde = False, label='Node 1')
    #sns.distplot(df['matched_Node2_dist'], hist = True, kde = False, label='Node 2')
    sns.displot(data = df, x = 'matched_Node1_dist', kind = "hist", label='Node 1')
    sns.displot(data = df, x = 'matched_Node2_dist', kind = "hist", label='Node 2')

    # Plot formatting
    plt.legend(prop={'size': 12})
    plt.title(' '.join(['Distribution of distance between endpoint and matched node: ', pathwayID]))
    plt.xlabel('Matched distance')
    plt.ylabel('Frequency')

    sns.despine()
    plt.tight_layout()
    plt.savefig(plotName, bbox_inches='tight', transparent=True)
    plt.clf()

def makeCombinedEndPointBarPlot(interactDFs=None):
    if interactDFs is None:
        interactDFs = glob.glob("*interactDF.csv")
    #print(interactDFs)
    concatInteractDF = pd.DataFrame()
    for df in interactDFs:
        pathwayID = re.sub("_interactDF.csv",  "", df)
        #print(pathwayID)
        df = pd.read_csv(df, engine = "python", sep = None, escapechar = "\\")
        df = df.assign(PathwayID = str(pathwayID))
        if len(concatInteractDF) == 0:
            concatInteractDF = df
        else:
            concatInteractDF = concatInteractDF.append(df)
            #concatInteractDF = concatInteractDF.reset_index()
    #print(concatInteractDF)
    #pd.write_csv(concatInteractDF,"concatInteractDF.csv", sep="\t")
    
    concatInteractDF = concatInteractDF.assign(index = concatInteractDF.index)

    makeCombinedInteractionBarPlot(concatInteractDF)

    node_count = pd.DataFrame(concatInteractDF.groupby(by = "PathwayID")['matched_Node1_cat'].value_counts(normalize=True))
    node_count['matched_Node2_cat'] = concatInteractDF.groupby(by = "PathwayID")['matched_Node2_cat'].value_counts(normalize=True)
    node_count = node_count.melt(ignore_index=False)
    node_count.index.names = ["PathwayID", "Category"]
    node_count = node_count.reset_index()

    #Make plot
    sns.catplot(x="value", y = "Category", data=node_count, height=6,  kind = "bar", legend=False, color = "grey", dodge = True, estimator = np.median, ci="sd").set(xlim=(0,1))
    #Label plot
    #plt.legend(prop={'size': 12}, title = 'Endpoint of interaction')
    plt.title(' '.join(['Endpoint Categories']))
    plt.xlabel('Proportion', fontsize=12)
    plt.ylabel('Category', fontsize=12)
    
    #Finish plot
    sns.despine()
    plt.tight_layout()
    plt.savefig("combined_endpoint_categories_1.png", bbox_inches='tight', transparent=True)
    plt.clf()

    #Make a plot to show whether the matched node ids are the same as the explicitly described node ids, ie, 'information added'
    temp = ["Explicit match"  if list(concatInteractDF.node1_graphref)[i] == list(concatInteractDF.matched_Node1)[i] else "Inferred match" for i in range(len(concatInteractDF)) ]
    #temp = 
    concatInteractDF = concatInteractDF.assign(Node1_match = temp)
    concatInteractDF = concatInteractDF.assign(Node2_match= np.where(concatInteractDF.node2_graphref == concatInteractDF.matched_Node2, "Explicit match", "Inferred match"))
    explicitMatch = pd.DataFrame(concatInteractDF.groupby(by = "PathwayID")['Node1_match'].value_counts(normalize=True))
    explicitMatch['Node2_match'] = concatInteractDF.groupby(by = "PathwayID")['Node2_match'].value_counts(normalize=True)
    explicitMatch.index.names = ["PathwayID", "ExplicitMatch"]
    explicitMatch = explicitMatch.melt(ignore_index=False)
    
    explicitMatch = explicitMatch.reset_index()

    sns.catplot(x="value", y = "ExplicitMatch", data=explicitMatch,  height=6,  kind = "bar", legend=False, color = "grey", dodge = True, estimator = np.median, ci="sd").set(xlim=(0,1)) #hue="variable",
    #Label plot
    #plt.legend(prop={'size': 12}, title = 'Endpoint of interaction')
    plt.title(' '.join(['Endpoint Categories']))
    plt.xlabel('Proportion', fontsize=12)
    plt.ylabel('Category', fontsize=12)
    
    #Finish plot
    sns.despine()
    plt.tight_layout()
    plt.savefig("combined_endpoint_categories_2.png", bbox_inches='tight', transparent=True)
    plt.clf()


def makeEndPointBarPlot(featureDFs, pathwayID):
    # Plotting a bar graph of the categories of nodes to which endpoints are mapped

    plotName = ''.join([pathwayID, "_endpoint_categories.svg"])

    #Endpoint 1
    node_count1  = pd.DataFrame(featureDFs['interactDF']['matched_Node1_cat'].value_counts())
    node_count1 = node_count1.rename(columns={'matched_Node1_cat':"Endpoint 1"})

    #Endpoint 2
    node_count2  = pd.DataFrame(featureDFs['interactDF']['matched_Node2_cat'].value_counts())
    node_count2 = node_count2.rename(columns={'matched_Node2_cat':"Endpoint 2"})

    #Prepare dataframe for plotting
    node_count = pd.concat([node_count1, node_count2], axis=1, sort=True)
    node_count = pd.melt(node_count.reset_index(), id_vars='index', value_vars = ['Endpoint 1', 'Endpoint 2'])

    #Make plot
    sns.catplot(x="index", y="value", hue="variable", data=node_count, height=6, kind="bar", palette="muted", legend=False)

    #Label plot
    plt.legend(prop={'size': 12}, title = 'Endpoint of interaction')
    plt.title(' '.join(['Endpoint Categories: ', pathwayID]))
    plt.ylabel('Frequencies', fontsize=12)
    plt.xlabel('Category', fontsize=12)

    #Finish plot
    sns.despine()
    plt.tight_layout()
    plt.savefig(plotName, bbox_inches='tight', transparent=True)
    plt.clf()


def makeCombinedInteractionBarPlot(concatInteractDF):
    concatInteractDF  = concatInteractDF.assign(OriginalEdge = [True if len(str(edge)) == 5 else False for edge in list(concatInteractDF.edgeID)])
    edge_count = pd.DataFrame(concatInteractDF.groupby(by = "PathwayID")['OriginalEdge'].value_counts(normalize=True))
    #edge_count.index.names = ["PathwayID", "OriginalEdge"]
    edge_count = edge_count.melt(ignore_index=False)
    edge_count = edge_count.reset_index()
    edge_count.OriginalEdge = [str("Explicit Edge") if temp == True else str("Inferred Edge") for temp in edge_count.OriginalEdge]
    sns.catplot(x="value", y = "OriginalEdge", data=edge_count,  kind = "bar", legend=False, color = "grey", dodge = True, estimator = np.median, ci="sd").set(xlim=(0,1)) #hue="variable",
    
    #Label plot
    plt.title(' '.join(['Edge Categories']))
    plt.xlabel('Median Proportion Across All Test Networks', fontsize=12)
    plt.ylabel('', fontsize=12) #Is Explicit Edge

    #Finish plot
    sns.despine()
    plt.tight_layout()
    plt.savefig("combined_edge_categories.png", bbox_inches='tight', transparent=True)
    plt.clf()

    #Break up inferred edges into categories
    

    return(None)
	
def unitTest(name, pathwayID):
    gpml = open(str(name), 'r')
    gpml = gpml.read()
    gpml = str(gpml)

    featureList = makeFeatureLists(gpml, isFromFile=True)
    featureDFs = getFeatureDFs(featureList)
    #pipeline
    featureDFs['interactDF'] = mapEndPoints(featureDFs)
    featureDFs['interactDF'] = mapPointsToInteractions(featureDFs, featureList)
    featureDFs['interactDF'] = processInteractDF(featureDFs, featureList)
    
    #make diagnostic plots
    makeMatchedDistHist(featureDFs, pathwayID)
    makeEndPointBarPlot(featureDFs, pathwayID)
    #write out graphml
    graph = makeGraph(featureDFs, featureList)
    #write out edgelist
    featureDFs['interactDF'].to_csv('_'.join([pathwayID, "interactDF.csv"]))
    featureDFs['datanodeDF'].to_csv('_'.join([pathwayID, "datanodeDF.csv"]))
    edgeListToCSV(graph, pathwayID)
    #Write out graphml
    #for source, target, edgeData in graph.edges(data=True):
    #    print([source, target, edgeData])
    nx.write_graphml_lxml(graph, ''.join([pathwayID, "_graph.graphml"]))
    return graph
def evaluateProgramOutput():
    programOuts = ["WP49_edgeList.csv", "WP127_edgeList.csv", "WP195_edgeList.csv", "WP205_edgeList.csv", "WP231_edgeList.csv", "WP286_edgeList.csv", "WP23_edgeList.csv", "WP3935_edgeList.csv", "WP4482_edgeList.csv", "WP4868_edgeList.csv", 'WP4462_edgeList.csv', 'WP2038_edgeList.csv', 'WP453_edgeList.csv','WP545_edgeList.csv', 'WP306_edgeList.csv'] #, "test_1_edgeList.csv"
    manualOuts = ["WP49_edgeList_edited.csv", "WP127_edgeList_edited.csv", "WP195_edgeList_edited.csv", "WP205_edgeList_edited.csv",  "WP231_edgeList_edited.csv",  "WP286_edgeList_edited.csv", "WP23_edgeList_edited.csv", "WP3935_edgelist_edited.csv", "WP4482_edgelist_edited.csv", "WP4868_edgelist_edited.csv", 'WP4462_edgelist_edited.csv', 'WP2038_edgelist_edited.csv', 'WP453_edgelist_edited.csv','WP545_edgelist_edited.csv', 'WP306_edgelist_edited.csv']
    fileNames = pd.DataFrame(list(zip(programOuts, manualOuts)), columns = ["programOuts", "manualOuts"])
    fileList = dict(zip(programOuts, manualOuts))

    netDiffResult = pd.DataFrame()

    for programOut in fileList.keys():
        print(programOut)
        manualOut = fileList[programOut]
        netDiff = getNetworkDiffs(manualOut, programOut)
        if netDiffResult.shape == (0, 0):
            netDiffResult = pd.DataFrame(netDiff).T
        else:
            netDiffResult = netDiffResult.append(pd.DataFrame(netDiff).T)
        #print(netDiffResult)

    netDiffResult.index = programOuts
    netDiffResult.columns = ["number_of_manualNodes","number_of_programNodes","manualNodes","programNodes","extra_manualNodes", "extra_programNodes","number_of_extra_manualEdges","number_of_extra_programEdges","extraEdgesInManual", "extraEdgesInProgram", "propTrueNodes", "propTrueEdges", "truePositiveEdges", "recall", "precision","truePositiveRate","falsePositiveRate", "f1Score"]
    #print(netDiffResult)
    netDiffResult.to_csv("netDiffResult_wikinetworks.tsv", sep= "\t")
    makeEvaluationPlots(netDiffResult)
    makeCombinedEndPointBarPlot(interactDFs=None)

def compareCytoscapeWithProgramOutput():
    programOuts = ["WP49_edgeList.csv", "WP127_edgeList.csv", "WP195_edgeList.csv", "WP205_edgeList.csv", "WP231_edgeList.csv", "WP286_edgeList.csv", "WP23_edgeList.csv", "WP3935_edgeList.csv", "WP4868_edgeList.csv", "WP4482_edgeList.csv", 'WP4462_edgeList.csv', 'WP2038_edgeList.csv', 'WP453_edgeList.csv','WP545_edgeList.csv', 'WP306_edgeList.csv']
    manualOuts = ["WP49_cytoscape.csv", "WP127_cytoscape.csv", "WP195_cytoscape.csv", "WP205_cytoscape.csv", "WP231_cytoscape.csv", "WP286_cytoscape.csv", "WP23_cytoscape.csv", "WP3935_cytoscape.csv", "WP4868_cytoscape.csv", "WP4482_cytoscape.csv", 'WP4462_cytoscape.csv', 'WP2038_cytoscape.csv', 'WP453_cytoscape.csv','WP545_edgelist_edited.csv', 'WP306_edgelist_edited.csv'] 
    fileNames = pd.DataFrame(list(zip(programOuts, manualOuts)), columns = ["programOuts", "cytoscapeOuts"])
    fileList = dict(zip(programOuts, manualOuts))

    netDiffResult = pd.DataFrame()

    for programOut in fileList.keys():
        #print(programOut)
        manualOut = fileList[programOut]
        netDiff = getNetworkDiffs(manualOut, programOut)
        if netDiffResult.shape == (0, 0):
            netDiffResult = pd.DataFrame(netDiff).T
        else:
            netDiffResult = netDiffResult.append(pd.DataFrame(netDiff).T)
        #print(netDiffResult)

    netDiffResult.index = programOuts
    netDiffResult.columns = ["number_of_cytoscapeNodes","number_of_programNodes","cytoscapeNodes","programNodes","extra_cytoscapeNodes", "extra_programNodes","number_of_extra_cytoscapeEdges","number_of_extra_programEdges","extraEdgesIncytoscape", "extraEdgesInProgram", "propTrueNodes", "propTrueEdges", "truePositiveEdges", "recall", "precision","truePositiveRate","falsePositiveRate", "f1Score"]
    #print(netDiffResult)
    netDiffResult.to_csv("netDiffResult_wikinetworksCytoscape.tsv", sep= "\t")

def createCombinedEdgelists():
    """Combine cytoscape and wikinetworks outputs to create a consensus that can be manually improved for evaluation of both programs"""
    cytoscapeEdgelists = glob.glob("*_cytoscape.csv")
    wikinetworksEdgelists = glob.glob("*_edgeList.csv")
    for cel in cytoscapeEdgelists:
        consensusFileName = str(cel)[:-14]+"_consensusEdgeList.csv"
        print(consensusFileName)
        wel = str(cel)[:-14]+"_edgeList.csv"
        wel = pd.read_csv(wel, escapechar = "\\")
        wel = wel[["Source", "Target", "signal"]]
        wel["Program"] = "Wikinetworks"
        cel = pd.read_csv(cel, escapechar = "\\")
        cel = cel[["Source", "Target", "signal"]]
        cel["Program"] = "Cytoscape"
        consensus_edgelist = wel.append(cel, ignore_index=True)
        consensus_edgelist = consensus_edgelist[~consensus_edgelist.Source.str.contains("group", na=False)]
        consensus_edgelist = consensus_edgelist[~consensus_edgelist.Target.str.contains("group", na=False)]
        consensus_edgelist = consensus_edgelist.drop_duplicates(subset = ["Source", "Target", "signal"])
        consensus_edgelist.to_csv(consensusFileName, quoting=csv.QUOTE_ALL, index=False)

def evaluateCytoscapeOutput():
    programOuts = ["WP49_cytoscape.graphml", "WP127_cytoscape.graphml", "WP195_cytoscape.graphml", "WP205_cytoscape.graphml", "WP231_cytoscape.graphml", "WP286_cytoscape.graphml", "WP23_cytoscape.graphml", "WP3935_cytoscape.graphml", "WP4868_cytoscape.graphml", "WP4482_cytoscape.graphml", 'WP4462_cytoscape.graphml', 'WP2038_cytoscape.graphml', 'WP453_cytoscape.graphml','WP545_cytoscape.graphml', 'WP306_cytoscape.graphml']
    #Generate edgelists from cytoscape graphmls
    for i in range(0, len(programOuts)):
        graph = programOuts[i]
        programOuts[i] = str(graph)[:-8]+".csv"
        graph = nx.read_graphml(graph)
        nodeNames = nx.get_node_attributes(graph,"name") 
        graphIDs = nx.get_node_attributes(graph,"GraphID")
        for temp in list(graph.nodes):
            if not temp in nodeNames.keys():
                nodeNames[temp] = graphIDs[temp]
            if nodeNames[temp] == "":
                nodeNames[temp] = "unnamed"
        nx.relabel_nodes(graph, mapping = nodeNames, copy=False)
        targetArrows = nx.get_edge_attributes(graph, "Target Arrow Shape")
        #print(targetArrows)
        for temp in list(graph.edges):
            if not temp in targetArrows.keys():
                targetArrows[temp] = "Unknown"
        nx.set_edge_attributes(graph, targetArrows, "targetArrows")
        nx.write_edgelist(graph, programOuts[i], delimiter=",", data = ["targetArrows"])
        print(programOuts[i])
        #tempDF = pd.read_csv(programOuts[i], header=None, sep=None, escapechar = "\\", engine = "python")
        tempDF = nx.to_pandas_edgelist(graph, source = 'Source', target = 'Target')
        tempDF = tempDF[['Source', 'Target', 'targetArrows']]
        print(tempDF.head())
        tempDF.columns = ["Source", "Target", "signal"]
        tempDF.signal = [checkArrow(arrow)[1] for arrow in tempDF["signal"]]
        tempDF.to_csv(programOuts[i], quoting=csv.QUOTE_ALL)
        print(tempDF)
    
    manualOuts = ["WP49_edgeList_edited.csv", "WP127_edgeList_edited.csv", "WP195_edgeList_edited.csv", "WP205_edgeList_edited.csv",  "WP231_edgeList_edited.csv",  "WP286_edgeList_edited.csv", "WP23_edgeList_edited.csv", "WP3935_edgelist_edited.csv", "WP4868_edgelist_edited.csv", "WP4482_edgelist_edited.csv", 'WP4462_edgelist_edited.csv', 'WP2038_edgelist_edited.csv', 'WP453_edgelist_edited.csv', 'WP545_edgelist_edited.csv', 'WP306_edgelist_edited.csv']
    fileNames = pd.DataFrame(list(zip(programOuts, manualOuts)), columns = ["programOuts", "manualOuts"])
    fileList = dict(zip(programOuts, manualOuts))
    netDiffResult = pd.DataFrame()
    for programOut in fileList.keys():
        print("###")
        print(programOut)
        manualOut = fileList[programOut]
        print(manualOut)
        netDiff = getNetworkDiffs(manualOut, programOut)
        if netDiffResult.shape == (0, 0):
            netDiffResult = pd.DataFrame(netDiff).T
        else:
            netDiffResult = netDiffResult.append(pd.DataFrame(netDiff).T)
        #print(netDiffResult)

    netDiffResult.index = programOuts
    netDiffResult.columns = ["number_of_manualNodes","number_of_programNodes","manualNodes","programNodes","extra_manualNodes", "extra_programNodes","number_of_extra_manualEdges","number_of_extra_programEdges","extraEdgesInManual", "extraEdgesInProgram", "propTrueNodes", "propTrueEdges", "truePositiveEdges", "recall", "precision","truePositiveRate","falsePositiveRate", "f1Score"]
    #print(netDiffResult)
    netDiffResult.to_csv("netDiffResult_cytoscape.tsv", sep= "\t")
    makeEvaluationPlots(netDiffResult, softwareName="cytoscape")

def evaluateCytoscapeOutput_fairplay():
    programOuts = ["WP49_cytoscape.graphml", "WP127_cytoscape.graphml", "WP195_cytoscape.graphml", "WP205_cytoscape.graphml", "WP231_cytoscape.graphml", "WP286_cytoscape.graphml", "WP23_cytoscape.graphml", 'WP4462_cytoscape.graphml', 'WP2038_cytoscape.graphml', 'WP453_cytoscape.graphml', 'WP545_cytoscape.graphml', 'WP306_cytoscape.graphml']
    #Generate edgelists from cytoscape graphmls
    for i in range(len(programOuts)):
        graph = programOuts[i]
        programOuts[i] = "fairplay_"+str(graph)[:-8]+".csv"
        graph = nx.read_graphml(graph)
        nodeNames = nx.get_node_attributes(graph,"name") 
        graphIDs = nx.get_node_attributes(graph,"GraphID")
        for temp in list(graph.nodes):
            if not temp in nodeNames.keys():
                nodeNames[temp] = graphIDs[temp]
        nx.relabel_nodes(graph, mapping = nodeNames, copy=False)
        targetArrows = nx.get_edge_attributes(graph, "Target Arrow Shape")
        #print(targetArrows)
        for temp in list(graph.edges):
            if not temp in targetArrows.keys():
                targetArrows[temp] = "Unknown"
        nx.set_edge_attributes(graph, targetArrows, "targetArrows")
        for temp in list(graph.nodes):
            if "group" in temp:
                #print(temp)
                graph = passThroughUnlabeled(temp, graph)
        nx.write_edgelist(graph, programOuts[i], delimiter="\t", data = ["targetArrows"]) # 
        #print(programOuts[i])
        tempDF = pd.read_csv(programOuts[i], header=None, sep="\t", escapechar = "\\")
        tempDF.columns = ["Source", "Target", "signal"]
        tempDF.signal = [checkArrow(arrow)[1] for arrow in tempDF["signal"]]
        tempDF.to_csv(programOuts[i], quoting=csv.QUOTE_ALL)
        #print(tempDF)
    
    manualOuts = ["WP49_edgeList_edited.csv", "WP127_edgeList_edited.csv", "WP195_edgeList_edited.csv", "WP205_edgeList_edited.csv",  "WP231_edgeList_edited.csv",  "WP286_edgeList_edited.csv", "WP23_edgeList_edited.csv", 'WP4462_edgeList_edited.csv', 'WP2038_edgeList_edited.csv', 'WP453_edgeList_edited.csv', 'WP545_edgelist_edited.csv', 'WP306_edgelist_edited.csv']
    fileNames = pd.DataFrame(list(zip(programOuts, manualOuts)), columns = ["programOuts", "manualOuts"])
    fileList = dict(zip(programOuts, manualOuts))
    netDiffResult = pd.DataFrame()
    for programOut in fileList.keys():
        #print(programOut)
        manualOut = fileList[programOut]
        netDiff = getNetworkDiffs(manualOut, programOut)
        if netDiffResult.shape == (0, 0):
            netDiffResult = pd.DataFrame(netDiff).T
        else:
            netDiffResult = netDiffResult.append(pd.DataFrame(netDiff).T)
        #print(netDiffResult)

    netDiffResult.index = programOuts
    netDiffResult.columns = ["number_of_manualNodes","number_of_programNodes","manualNodes","programNodes","extra_manualNodes", "extra_programNodes","number_of_extra_manualEdges","number_of_extra_programEdges","extraEdgesInManual", "extraEdgesInProgram", "propTrueNodes", "propTrueEdges", "truePositiveEdges", "recall", "precision","truePositiveRate","falsePositiveRate", "f1Score"]
    #print(netDiffResult)
    netDiffResult.to_csv("netDiffResult_cytoscape_fairplay.tsv", sep= "\t")
    makeEvaluationPlots(netDiffResult, softwareName="cytoscape_fairplay")

def makeEvaluationPlots(netDiffResult, softwareName = "wikinetworks"): 
    #sns.set_theme(context = "paper", style="whitegrid")
    #Proportion of true edges = number(intersection of manual and program edges)/number(program edges)
    netDiffResult["Pathway"] = [re.sub("_edgeList.csv", "", temp) for temp in netDiffResult.index]
    ax = sns.barplot(x = "propTrueEdges", y = "Pathway", data = netDiffResult, color="#6495ED")
    plt.xlim(0, 1)
    ax.figure.savefig("propTrueEdges"+"_"+str(softwareName)+".png")
    plt.clf()
    #Proportion of true nodes = number(intersection of manual and program nodes)/number(program nodes)
    ax = sns.barplot(x = "propTrueNodes", y = "Pathway", data = netDiffResult, color="#6495ED")
    ax.figure.savefig("propTrueNodes"+"_"+str(softwareName)+".png")
    plt.xlim(0, 1)
    plt.clf()
    #Number of edges left out
    #netDiffResult["number_of_extra_programEdges"]
    ax = sns.barplot(y = "number_of_extra_programEdges", x = "Pathway", data = netDiffResult, color="#6495ED")
    ax.figure.savefig("number_of_extra_programEdges"+"_"+str(softwareName)+".png")
    plt.clf()
    #Number of nodes left out
    #netDiffResult["extra_programNodes"]
    ax = sns.barplot(y = "extra_programNodes", x = "Pathway", data = netDiffResult, color="#6495ED")
    ax.figure.savefig("extra_programNodes"+"_"+str(softwareName)+".png")
    plt.clf()
    #Number of edges added
    #netDiffResult["number_of_extra_manualEdges"]
    ax = sns.barplot(y = "number_of_extra_manualEdges", x = "Pathway", data = netDiffResult, color="#6495ED")
    ax.figure.savefig("number_of_extra_manualEdges"+"_"+str(softwareName)+".png")
    plt.clf()
    #Number of nodes added
    #netDiffResult["extra_manualNodes"]
    ax = sns.barplot(y = "extra_manualNodes", x = "Pathway", data = netDiffResult, color="#6495ED")
    ax.figure.savefig("extra_manualNodes"+"_"+str(softwareName)+".png")
    plt.clf()
    #Number of true positive edges
    ax = sns.barplot(y = "truePositiveEdges", x = "Pathway", data = netDiffResult, color="#6495ED")
    ax.figure.savefig("truePositiveEdges"+"_"+str(softwareName)+".png")
    plt.clf()
    #ROC curve: x = false positive rate, y = true positive rate
    #false positive = edge where there is no edge, ie, the edge is in program but not in manual, ie, extra edges in program = "number_of_extra_programEdges"
    #true positive = edge where there is edge in both program and in manual
    #false negative = missing edge in program, ie, extra edges in manual, = "number_of_extra_manualEdges"
    #true negative = correctly predicted no edge when there is no edge, (n-1)*n/2 where n = number of nodes
    #True positive rate = true positives/(true positives + false negatives)
    #False positive rate = false positives/(true negatives + false positives)
    netDiffResult.truePositiveRate = netDiffResult.truePositiveRate.astype(float)
    netDiffResult.falsePositiveRate = netDiffResult.falsePositiveRate.astype(float)
    fig, ax = plt.subplots()
    sns.scatterplot(x = "falsePositiveRate", y = "truePositiveRate", data = netDiffResult, color="blue", ax = ax)
    sns.lineplot(x = "falsePositiveRate", y = "truePositiveRate", data = netDiffResult, color="red", ax = ax)
    ax.figure.savefig("ROC_curve"+"_"+str(softwareName)+".png")
    plt.show()
    plt.clf()    

    #Precision-Recall curve: x = recall, y = precision
    #Recall = True Positives / (True Positives + False Negatives)
    #Precision = True Positives / (True Positives + False Positives)

    #recall = intersection between edges in program and manual/(intersection between edges in program and manual + "number_of_extra_manualEdges")
    #precision = intersection between edges in program and manual/(intersection between edges in program and manual + "number_of_extra_programEdges")
    netDiffResult.precision = netDiffResult.precision.astype(float)
    netDiffResult.recall = netDiffResult.recall.astype(float)
    fig, ax = plt.subplots()

    sns.scatterplot(x = "recall", y = "precision", data = netDiffResult, color="blue", ax = ax)
    #sns.lineplot(x = "recall", y = "precision", data = netDiffResult, color="red", ax = ax)
    for i in range(netDiffResult.shape[0]):
        plt.text(x=netDiffResult.recall[i],y=netDiffResult.precision[i],s=netDiffResult.Pathway[i], fontdict=dict(color='red',size=6), bbox=dict(facecolor='white',alpha=0.5))
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    ax.figure.savefig("precision_recall_curve"+"_"+str(softwareName)+".png")
    plt.clf()
    return None

def testToyPathways():
    
    testFiles = glob.glob("test*.gpml")
    
    #print(testFiles)
    
    for tf in testFiles:
        #print(tf)
        pathwayID = tf[:-5]
        #print(pathwayID)
        unitTest(tf, pathwayID)
    return None

def makeComparisonPlots():
    cytoscape = pd.read_csv("netDiffResult_cytoscape.tsv", sep= "\t", index_col=0)
    cytoscape["Method"] = "cytoscape"
    #cytoscapeFairplay = pd.read_csv("netDiffResult_cytoscape_fairplay.tsv", sep= "\t", index_col=0, escapechar = "\\")
    #cytoscapeFairplay["Method"] = "cytoscape_fairplay"
    wikinetworks = pd.read_csv("netDiffResult_wikinetworks.tsv", sep= "\t", index_col=0, escapechar = "\\")
    wikinetworks["Method"] = "wikinetworks"
    combinedDF =  pd.concat([wikinetworks, cytoscape]) #, cytoscapeFairplay])
    #combinedDF["propTrueEdges"] = [100*temp for temp in combinedDF["propTrueEdges"]]
    #print(combinedDF)
    combinedDF["Pathway"] = [re.sub("_edgeList.csv", "", temp) for temp in combinedDF.index]
    combinedDF["Pathway"] = [re.sub("_cytoscape.csv", "", temp) for temp in combinedDF["Pathway"]]
    #combinedDF["Pathway"] = [re.sub("fairplay_", "", temp) for temp in combinedDF["Pathway"]]
    sns.set_context("paper")#, rc={"font.size":10,"axes.titlesize":12,"axes.labelsize":12})
    #plt.figure(figsize=(8,5))
    ax = sns.barplot(x = "precision", y = "Pathway", hue="Method", data = combinedDF, palette="Blues")
    plt.xlim(0, 1)
    plt.xlabel("Precision")
    ax.figure.savefig("precisionBarplot"+"_"+str("combined")+".png")
    plt.clf()
    
    sns.set_context("paper")#, rc={"font.size":10,"axes.titlesize":12,"axes.labelsize":12})
    plt.figure(figsize=(10,5))
    ax = sns.barplot(x = "precision", y="Method", data = combinedDF, palette="Blues", ci="sd")
    plt.xlim(0, 1)
    plt.xlabel("Precision")
    ax.figure.savefig("precision"+"_"+str("combined")+".png")
    plt.clf()

    #"f1Score"
    sns.set_context("paper")
    ax = sns.barplot(x = "f1Score", y = "Pathway", hue="Method", data = combinedDF, palette="Purples")
    plt.xlim(0, 1)
    plt.xlabel("f1Score")
    ax.figure.savefig("f1ScoreBarplot"+"_"+str("combined")+".png")
    plt.clf()

    sns.set_context("paper")
    plt.figure(figsize=(5,10))
    ax = sns.barplot(y = "f1Score", x ="Method", data = combinedDF, color = "grey", ci="sd")
    plt.ylim(0, 1)
    plt.ylabel("F1 Score")
    plt.xlabel("Methods")
    ax.figure.savefig("f1ScoreStdDev"+"_"+str("combined")+".png")
    plt.clf()

    #recall
    sns.set_context("paper")
    ax = sns.barplot(x = "recall", y = "Pathway", hue="Method", data = combinedDF, palette="OrRd")
    plt.xlim(0, 1)
    plt.xlabel("Recall")
    ax.figure.savefig("RecallBarplot"+"_"+str("combined")+".png")
    plt.clf()

    sns.set_context("paper")
    plt.figure(figsize=(10,5))
    ax = sns.barplot(x = "recall", y="Method", data = combinedDF, palette="OrRd", ci="sd")
    plt.xlim(0, 1)
    plt.xlabel("recall")
    ax.figure.savefig("recallStdDev"+"_"+str("combined")+".png")
    plt.clf()    
    
    #precision-recall curve
    sns.set_context("paper")
    plt.figure(figsize=(10,5))
    fig, ax = plt.subplots()
    sns.scatterplot(x = "recall", y = "precision", data = combinedDF, hue="Method", palette = "Set2", ax = ax)
    #for i in range(combinedDF.shape[0]):
    #    plt.text(x=combinedDF.recall[i],y=combinedDF.precision[i],s=combinedDF.Pathway[i], fontdict=dict(color='red',size=6), bbox=dict(facecolor='white',alpha=0.5))
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    ax.figure.savefig("precision_recall_curve"+"_"+str("Combined")+".png")
    plt.clf()

    #catplot combining precision, recall, F1 score
    sns.set_context("paper")
    plt.figure(figsize=(10,5))
    cDF2 = combinedDF.melt(id_vars = ["Pathway", "Method"], value_vars = ["precision", "recall", "f1Score"])
    cDF2.columns = ["Pathway", "Method", "Metric", "Value"]
    ax = sns.catplot(x = "Value", y="Metric", hue =  "Method", data = cDF2, palette="OrRd", ci="sd", kind="bar")
    plt.xlim(0, 1)
    ax.savefig("all_metrics"+"_"+str("Combined")+".png")
    plt.clf()


if __name__ == '__main__':
    __all__ = ["WikiPathways"]
    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.max_rows', 10000)

    s=WikiPathways()
    testPathways = ['WP4482', 'WP3935', 'WP4868', 'WP49', 'WP127', 'WP195', 'WP205', 'WP231', 'WP286', 'WP23', 'WP4462', 'WP2038', 'WP453', 'WP545', 'WP306'] #['WP254', 'WP23', "WP4495", "WP437",  "WP4495", "WP195", "WP49", "WP286", "WP395", "WP127", "WP364", "WP205", "WP231", "WP366"] # "WP1835", "WP1840",  "WP1919", "WP1836",
    
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
    compareCytoscapeWithProgramOutput()
    createCombinedEdgelists()
