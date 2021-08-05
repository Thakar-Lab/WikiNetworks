import networkx as nx
import re
import urllib
import csv
import itertools as it
import sys
from bs4 import BeautifulSoup
from random import randint, sample, choice
import requests
import binascii
from bioservices import WikiPathways
import numpy as np
from shapely.geometry import *
from scipy.spatial import cKDTree, distance_matrix
import pandas as pd
import fiona
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from string import ascii_uppercase, digits
from unidecode import unidecode
import collections

def id_generator(size=5, chars=ascii_uppercase + digits):
    return ''.join(choice(chars) for _ in range(size))

def givenCodeGetGPML(s, code):
    """download gpml files from wikipathways database"""
    code=code.decode()
    res=s.getPathwayAs(code, filetype="gpml")
    #newres = unidecode(res)
    newres=binascii.a2b_base64(bytes(res,"ascii")) # convert to ascii
    return(newres)

def getCurationTags(s, code):
    url="https://webservice.wikipathways.org/getCurationTags?pwId="+str(code)+"&format=xml"
    res=s.http_get(url + "&format=json")
    return(res)
    
def processCurationTags(curationTagsDict):
    """Processes the output of getCurationTags to get a readable list of warning labels. The user should consider this list of tags as a quality indicator."""
    #get the dictionary of tags
    #print the pathway information
    print("\n***PATHWAY INFORMATION***\n")
    print("Name: ", curationTagsDict['tags'][0]['name'], "\n")
    print("Display Name: ", curationTagsDict['tags'][0]['displayName'], "\n")
    for key in curationTagsDict['tags'][0]['pathway'].keys():
        print(str(key)+": "+curationTagsDict['tags'][0]['pathway'][str(key)]+"\n")
    print("\n***CURATION WARNINGS:***\n")
    for annotation in range(1, len(curationTagsDict['tags'])):
        print("\n***CURATION TAG "+str(annotation)+"***\n")
        print("Name: "+curationTagsDict['tags'][annotation]['name'])
        print("Display Name: "+curationTagsDict['tags'][annotation]['displayName'])
        print("Description: "+curationTagsDict['tags'][annotation]['text'])
    print("*********\n")

def getAnchor(interacts, featureDFs):
    node_to_anchor = {}
    anchor_to_node = {}
    for entry in interacts:
        nodes = entry.find_all('Point')
        anchors = entry.find_all('Anchor')
        node2 = nodes[len(nodes) - 1].get('GraphRef')
        node1 = nodes[0].get('GraphRef')
        if not anchors is None:
            #print("Anchors", anchors)
            i = 0
            for anchor in list(anchors): #
                i = i + 1
                anchor_graphid = anchor.get('GraphId')
                if anchor_graphid is None:
                    anchor_graphid = '_'.join([str(node1), str(node2), str(i)])
                #if not anchor_graphid is None and not node1 is None:
                #    anchor_to_node[anchor_graphid] = node1
                elif not anchor_graphid is None and not node2 is None and node1 is None:
                    #print("Get anchor", anchor_graphid, node2, node1)
                    anchor_to_node[anchor_graphid] = anchor_graphid #node2
                    node1DF = gpd.GeoDataFrame(featureDFs['interactDF'][featureDFs['interactDF'].anchor == str(anchor_graphid)], geometry=featureDFs['interactDF'].node1_coords) # add additional clause to remove matches to nodes that only have defined interactions
                    for feature in ["datanodeDF", "shapeDF"]:
                        distDF_node1 = ckdnearest(node1DF, featureDFs[feature])
                        #print(min(distDF_node1.dist))
                        temp = distDF_node1.matched_Node1.tolist()
                        temp = [x for x in temp if str(x) != 'nan']
                        #print(temp)
                        if len(temp) > 0:
                            anchor_to_node[anchor_graphid] = [x for x in temp]
                            #print(anchor_to_node[anchor_graphid])
                        elif len(temp) == 0:
                            anchor_to_node[anchor_graphid] = anchor_graphid
                    #node1DF = gpd.GeoDataFrame(featureDFs['interactDF'], geometry = featureDFs['interactDF'].node1_coords)
                    #node2DF = gpd.GeoDataFrame(featureDFs['interactDF'], geometry = featureDFs['interactDF'].node2_coords)
                    #distDF_node1 = ckdnearest_interact(node1DF, node2DF)
    #return(node_to_anchor)
    return(anchor_to_node)

def getNodeLoc(features):
    """get coordinates of either labels, shapes or nodes"""
    nodeData = []
    for entry in list(features):
        graphid=entry.get('GraphId')
        if graphid is None:
            graphid = id_generator()
        groupRef=entry.get('GroupRef')
        graphics=entry.find('Graphics')
        centerX=float(graphics.get("CenterX"))
        centerY=float(graphics.get("CenterY"))
        width=float(graphics.get("Width"))
        height=float(graphics.get("Height"))
        textlabel=entry.get('TextLabel')
        color=graphics.get('Color')
        color = ''.join(["#",str(color)])
        if color is None or color is "None":
            color = "#000000"
        #get database identifiers
        databaseInf=entry.find("Xref")
        if not databaseInf is None:
            database=databaseInf.get("Database")
            databaseID=databaseInf.get("ID")
            #if database information is not found, replace with textlabel
            if database is None or database is "":
                database = "Unknown"
            if databaseID is None or databaseID is "":
                databaseID = textlabel
        else:
            database = "Unknown"
            databaseID = textlabel
        nodeData.append([textlabel, graphid, groupRef, centerX, centerY, width, height, color, database, databaseID])
    nodeDF = pd.DataFrame(nodeData, columns = ['textlabel', 'graphid', 'groupref', 'centerX', 'centerY', 'width', 'height', 'color', 'database', 'databaseID'])
    nodeDF = gpd.GeoDataFrame(nodeDF, geometry = gpd.points_from_xy(nodeDF.centerX, nodeDF.centerY))
    return(nodeDF)

def makeFeatureLists(gpml, isFromFile = False):
    """read all lines into a bs4 object using libXML parser"""
    if isFromFile:
        soup = BeautifulSoup(gpml, "xml")    
    else:
        soup = BeautifulSoup(''.join(gpml.decode('utf-8','ignore')), 'xml')

    nodes = soup.find_all('DataNode')
    interacts = soup.find_all('Interaction')
    states = soup.find_all('State')
    labels = soup.find_all('Label')
    groups = soup.find_all('Group')
    shapes = soup.find_all('Shape')
    return([soup, nodes, interacts, states, labels, groups, shapes])

def getInteractionLoc(featureList):
    interacts = featureList[2]
    interactData = []
    node2Data = []
    for entry in interacts:
        graphid=entry.get('GraphId')
        groupRef=entry.get('GroupRef')
        graphics=entry.find('Graphics')
        color=graphics.get('Color')
        if color is None:
            color = "000000" #"Unknown"
        nodes=entry.find_all('Point')
        node1_graphref=nodes[0].get('GraphRef')
        node2_graphref=nodes[len(nodes)-1].get('GraphRef')
        node1_x=float(nodes[0].get('X'))
        node2_x=float(nodes[len(nodes)-1].get('X'))
        node1_y=float(nodes[0].get('Y'))
        node2_y=float(nodes[len(nodes)-1].get('Y'))
        node1_coords = Point(node1_x, node1_y)
        node2_coords = Point(node2_x, node2_y)
        arrow = nodes[len(nodes)-1].get('ArrowHead')
        geometry = MultiPoint([(node1_x, node1_y), (node2_x, node2_y)]).convex_hull
        comments = entry.find_all("Comment")
        anchor = entry.find("Anchor")
        if anchor is None:
            anchor = "None"
        else:
            anchor = anchor.get("GraphId")
        interactData.append([graphid, groupRef, node1_graphref, node1_x, node1_y, node1_coords, node2_graphref, node2_x, node2_y, node2_coords, arrow, geometry, comments, color, anchor])
    interactDF = gpd.GeoDataFrame(pd.DataFrame(interactData, columns=["edgeID", "edgeGroup", "node1_graphref", "node1_x", "node1_y", "node1_coords", "node2_graphref", "node2_x", "node2_y", "node2_coords", "arrow", "geometry", "comments", "color", "anchor"]))
    return(interactDF)

def joinGroups(nodeDF, interactDF):
    grouped_nodes = nodeDF.groupby('groupref')
    #print(pd.DataFrame(grouped_nodes))
    group_node_data = []
    for name, group in grouped_nodes:
        #print("name", name)
        #print("group", group)
        nodelist = group.graphid.tolist()
        #if len(nodelist) == 1:
        #    nodelist = [nodelist]
        groupRef = name
        arrow = "group"
        i = 1
        #print("nodelist", nodelist)
        for node1 in nodelist:
            #print("node1", node1)
            node1_graphref = node1
            node1_x = group.centerX[group.graphid == node1].astype(float)
            node1_y = group.centerY[group.graphid == node1].astype(float)
            #print(["node1 coords", node1_x.tolist(), node1_y.tolist()])
            if not len(node1_x.tolist()) == 0 and not len(node1_y.tolist()) == 0:
                node1_coords = Point(node1_x, node1_y)
                graphid = '_'.join(nodelist)+"_"+str(i)
                i = i + 1
                for node2 in nodelist:
                    if node1 is not node2:
                        #print(node1, node2)
                        node2_graphref = node2
                        node2_x = group.centerX[group.graphid == node2].astype(float)
                        node2_y = group.centerY[group.graphid == node2].astype(float)
                        node2_coords = Point(node2_x, node2_y)
                        geometry = MultiPoint([(node1_x, node1_y), (node2_x, node2_y)])
                        group_node_data.append([graphid, groupRef, node1_graphref, node1_x, node1_y, node1_coords, node2_graphref, node2_x, node2_y, node2_coords, arrow, geometry])
    group_nodeDF = gpd.GeoDataFrame(pd.DataFrame(group_node_data, columns=["edgeID", "edgeGroup", "node1_graphref", "node1_x", "node1_y", "node1_coords", "node2_graphref", "node2_x", "node2_y", "node2_coords", "arrow", "geometry"]))
    interactDF = pd.concat([interactDF, group_nodeDF], sort=True)
    interactDF.reset_index(inplace=True, drop=True)
    return(interactDF)

def ckdnearest2(gdfA, gdfB, gdfB_cols=['Place']):
    A = np.concatenate(
        [np.array(geom.coords) for geom in gdfA.geometry.to_list()])
    B = [np.array(geom.coords) for geom in gdfB.geometry.to_list()]
    B_ix = tuple(itertools.chain.from_iterable(
        [itertools.repeat(i, x) for i, x in enumerate(list(map(len, B)))]))
    B = np.concatenate(B)
    ckd_tree = cKDTree(B)
    dist, idx = ckd_tree.query(A, k=1)
    idx = itemgetter(*idx)(B_ix)
    gdf = pd.concat(
        [gdfA, gdfB.loc[idx, gdfB_cols].reset_index(drop=True),
         pd.Series(dist, name='dist')], axis=1)
    return gdf

def ckdnearest(gdA, gdB, k=1):
    nA = np.array(list(zip(gdA.geometry.x, gdA.geometry.y)))
    nB = np.array(list(zip(gdB.geometry.x, gdB.geometry.y)))
    btree = cKDTree(nB)
    dist, idx = btree.query(nA, k=k)

    gdf = pd.concat(
        [gdA, gdB.loc[idx, gdB.columns != 'geometry'].reset_index(),
         pd.Series(dist, name='dist')], axis=1)

    return gdf

def ckdnearest_interact(gdA, gdB):
    #(intentionally) Does not return matches with distance == 0

    nA = np.array(list(zip(gdA.geometry.x, gdA.geometry.y)))
    nB = np.array(list(zip(gdB.geometry.x, gdB.geometry.y)))

    btree = cKDTree(nB)

    dist, idx = btree.query(nA, k=[1,2])
    idx = np.array([temp[1] for temp in idx])
    dist = np.array([temp[1] for temp in dist])
    #print(idx)
    #print(gdA.columns)
    ids=gdA['edgeID'].tolist()
    ids=[ids[temp] for temp in list(idx)]

    #rename column for gdB edgeID
    gdB = gdB.rename(columns={'edgeID':'matchID', 'node1_graphref': 'matchUpstream', 'node2_graphref': 'matchDownstream'})
    #gdA = gdA.rename(columns={'edgeID':'matchID', 'node1_graphref': 'matchUpstream', 'node2_graphref': 'matchDownstream'})

    gdf = pd.concat(
        [gdA, gdB.loc[idx, gdB.columns != 'geometry'].reset_index(),
         pd.Series(dist, name='dist')], axis=1)

    return gdf

def getFeatureDFs(featureList):

    #featureList: soup, nodes, interacts, states, labels, groups, shapes

    #Make node dataframes
    datanodeDF = getNodeLoc(featureList[1])
    labelDF = getNodeLoc(featureList[4])
    shapeDF = getNodeLoc(featureList[6])

    #make group dataframe
    groupDF = getGroupDF(featureList[5])

    #make interaction dataframe
    interactDF = getInteractionLoc(featureList)
    interactDF = joinGroups(datanodeDF, interactDF)
    #return all feature dataframes
    featureDFs = {}
    featureDFs['datanodeDF']=datanodeDF
    featureDFs['labelDF']=labelDF
    featureDFs['shapeDF']=shapeDF
    featureDFs['groupDF']=groupDF
    featureDFs['interactDF']=interactDF
    featureDFs['interactDF'] = mapUnlabeledGroups(featureDFs)
    return(featureDFs)

def getGroupDF(groupFeatures):
    groupData = []
    #featureList: soup, nodes, interacts, states, labels, groups, shapes
    for entry in groupFeatures:
        groupid = entry.get("GroupId")
        graphid = entry.get("GraphId")
        if graphid is None:
            graphid = groupid
        else:
            graphid = graphid
        groupData.append([groupid, graphid])
    groupDF = pd.DataFrame(groupData, columns = ['groupid', 'graphid'])
    #print(groupDF.head(10))
    return(groupDF)

def mapEndPoints(featureDFs):

    """
        Map the endpoints of all interactions to either nodes, shapes, or labels.
        Caveat: maps to shapes, not their constituent nodes. TODO: add a recursive search as in V2, to map shape IDs to nodes and add edges from upstream to all nodes in shape.
        Does not handle the case where endpoints map to the endpoints of other interactions - this should be another function (TODO).
        Also does not explicitly have a distance cutoff - just returns the minimum distance and the corresponding mapped entity.
        An intuitive reject criterion is as follows:
            reject match if distance between interaction endpoint and center of the node is greater than (distance between any corner point and center + 10% of that distance) (TODO)
            AND
            reject match if an endpoint is explicitly stated in the interaction and the matched node does not match the explicit endpoint (TODO)
    """

    #Add columns with graphids of matched nodes and the closest distance of the matched nodes
    #For node 1:
    featureDFs['interactDF']['matched_Node1_dist'] = np.where(featureDFs['interactDF']['node1_graphref'].isnull(), np.inf, 0)
    featureDFs['interactDF']['matched_Node1'] = np.where(featureDFs['interactDF']['node1_graphref'].isnull(), np.nan, featureDFs['interactDF']['node1_graphref'])

    #For node 2:
    featureDFs['interactDF']['matched_Node2_dist'] = np.where(featureDFs['interactDF']['node2_graphref'].isnull(), np.inf, 0)
    featureDFs['interactDF']['matched_Node2'] = np.where(featureDFs['interactDF']['node2_graphref'].isnull(), np.nan, featureDFs['interactDF']['node2_graphref'])

    #Add column that describes whether the endpoint of an interaction is a shape/label/node/nan
    featureDFs['interactDF']['matched_Node1_cat'] =  "" #np.where(featureDFs['interactDF']['node1_graphref'].isnull(), np.nan, "explicit")
    featureDFs['interactDF']['matched_Node2_cat'] = "" #np.where(featureDFs['interactDF']['node2_graphref'].isnull(), np.nan, "explicit")

    for feature in featureDFs.keys():
        #print(feature)

        if not feature is 'interactDF' and not feature is 'groupDF':

            #Get groupIDs for each node
            groupIDs = featureDFs[feature]['groupref']

            #print(groupIDs.shape)
            if len(groupIDs) >= 1:
                #Map node 1
                node1DF = gpd.GeoDataFrame(featureDFs['interactDF'], geometry=featureDFs['interactDF'].node1_coords)
                distDF_node1 = ckdnearest(node1DF, featureDFs[feature])
                featureDFs['interactDF']['matched_Node1'] = np.where((featureDFs['interactDF']['matched_Node1_dist'] >= distDF_node1['dist']), distDF_node1['graphid'], featureDFs['interactDF']['matched_Node1'])

                featureDFs['interactDF']['matched_Node1_dist'] = np.where((featureDFs['interactDF']['matched_Node1_dist'] >= distDF_node1['dist']), distDF_node1['dist'], featureDFs['interactDF']['matched_Node1_dist'])

                featureDFs['interactDF']['matched_Node1_cat'] = np.where((featureDFs['interactDF']['matched_Node1_dist'] >= distDF_node1['dist']), feature, featureDFs['interactDF']['matched_Node1_cat'])

                #Map node 2
                node2DF = gpd.GeoDataFrame(featureDFs['interactDF'], geometry=featureDFs['interactDF'].node2_coords)
                distDF_node2 = ckdnearest(node2DF, featureDFs[feature])

                featureDFs['interactDF']['matched_Node2'] = np.where((featureDFs['interactDF']['matched_Node2_dist'] >= distDF_node2['dist']), distDF_node2['graphid'], featureDFs['interactDF']['matched_Node2'])

                featureDFs['interactDF']['matched_Node2_dist'] = np.where((featureDFs['interactDF']['matched_Node2_dist'] >= distDF_node2['dist']) , distDF_node2['dist'], featureDFs['interactDF']['matched_Node2_dist'])

                featureDFs['interactDF']['matched_Node2_cat'] = np.where((featureDFs['interactDF']['matched_Node2_dist'] >= distDF_node2['dist']), feature, featureDFs['interactDF']['matched_Node2_cat'])

    return(featureDFs['interactDF'])

def matchRef_alias(x,featureDFs, featureList, anchor_graphids):
    mRef = matchRef(x,featureDFs, featureList, anchor_graphids)
    if not mRef is None:
        return(mRef)
    else:
        return([])

def mapPointsToInteractions(featureDFs,featureList):

    """
        Similar to the function mapEndPoints.
        Instead of searching labels/shapes/datanodes for matches, search the endpoints of OTHER interactions for a match.
        This takes care of branched interactions.
    """
    
    #Create geoDataFrames for endpoint 1 and endpoint 2 of each interaction
    node1DF = gpd.GeoDataFrame(featureDFs['interactDF'], geometry = featureDFs['interactDF'].node1_coords)
    node2DF = gpd.GeoDataFrame(featureDFs['interactDF'], geometry = featureDFs['interactDF'].node2_coords)

    #Check if the first endpoint of any interaction maps to the second endpoint of any ***other*** interaction
    distDF_node1 = ckdnearest_interact(node1DF, node2DF)

    list1 = list(distDF_node1['matchUpstream'])
    list2 = list(node1DF['matched_Node1'])
    anchor_graphids = getAnchor(featureList[2], featureDFs)

    cond1 = [any(item in matchRef_alias(list2[temp], featureDFs, featureList, anchor_graphids) for item in matchRef_alias(list1[temp], featureDFs, featureList, anchor_graphids)) or any(item in matchRef_alias(list2[temp], featureDFs, featureList, anchor_graphids) for item in matchRef_alias(list1[temp], featureDFs, featureList, anchor_graphids))for temp in range(0, len(list2))]

    distDF_node1['dist'] = np.where(cond1, np.inf, distDF_node1['dist'])

    featureDFs['interactDF']['matched_Node1'] = np.where((featureDFs['interactDF']['matched_Node1_dist'] >= distDF_node1['dist']), distDF_node1['matchUpstream'], featureDFs['interactDF']['matched_Node1'])

    featureDFs['interactDF']['matched_Node1_dist'] = np.where(featureDFs['interactDF']['matched_Node1_dist'] >= distDF_node1['dist'], distDF_node1['dist'], featureDFs['interactDF']['matched_Node1_dist'])

    featureDFs['interactDF']['matched_Node1_cat'] = np.where(featureDFs['interactDF']['matched_Node1_dist'] >= distDF_node1['dist'], "edge", featureDFs['interactDF']['matched_Node1_cat'])


    #Check if the first endpoint of any interaction maps to the first endpoint of any ***other*** interaction

    distDF_node2 = ckdnearest_interact(node1DF, node2DF)

    list1 = list(distDF_node2['matchDownstream'])
    list2 = list(node2DF['matched_Node1'])

    cond1 = [any(item in matchRef_alias(list2[temp], featureDFs, featureList, anchor_graphids) for item in matchRef_alias(list1[temp], featureDFs, featureList, anchor_graphids)) or any(item in matchRef_alias(list2[temp], featureDFs, featureList, anchor_graphids) for item in matchRef_alias(list1[temp], featureDFs, featureList, anchor_graphids))for temp in range(0, len(list2))]

    distDF_node2['dist'] = np.where(cond1, np.inf, distDF_node2['dist'])

    featureDFs['interactDF']['matched_Node1'] = np.where((featureDFs['interactDF']['matched_Node1_dist'] >= distDF_node1['dist']), distDF_node1['matchDownstream'], featureDFs['interactDF']['matched_Node1'])

    featureDFs['interactDF']['matched_Node1_dist'] = np.where(featureDFs['interactDF']['matched_Node1_dist'] >= distDF_node1['dist'], distDF_node1['dist'], featureDFs['interactDF']['matched_Node1_dist'])

    featureDFs['interactDF']['matched_Node1_cat'] = np.where(featureDFs['interactDF']['matched_Node1_dist'] >= distDF_node1['dist'], "edge", featureDFs['interactDF']['matched_Node1_cat'])
    
    return(featureDFs['interactDF'])

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
        df = pd.read_csv(df, engine = "python", sep = None)
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

def matchRef(x, featureDFs, featureList, anchor_graphids):
    node_grouprefs = featureDFs['datanodeDF']['groupref'].tolist()
    node_graphids = featureDFs['datanodeDF']['graphid'].tolist()

    shape_grouprefs = featureDFs['shapeDF']['groupref'].tolist()
    shape_graphids = featureDFs['shapeDF']['graphid'].tolist()
    shape_textlabels = featureDFs['shapeDF']['textlabel'].tolist()

    label_grouprefs = featureDFs['labelDF']['groupref'].tolist()
    label_graphids = featureDFs['labelDF']['graphid'].tolist()
    label_textlabels = featureDFs['labelDF']['textlabel'].tolist()

    group_graphids = featureDFs['groupDF']['graphid'].tolist()
    group_groupids = featureDFs['groupDF']['groupid'].tolist()

    textlabels = featureDFs['datanodeDF']['textlabel'].tolist()

    if not x is None:
        if x in anchor_graphids:
            x = anchor_graphids[x]
        else:
            x = x
        if x in node_graphids:
            #preferentially return graphids
            return(textlabels[node_graphids.index(x)])

        elif x in node_grouprefs:
            #match groupref to graphid
            matchGraphID = featureDFs['datanodeDF']['graphid'][featureDFs['datanodeDF']['groupref'] == x]
            return(textlabels[node_grouprefs.index(matchGraphID)])
        elif x in shape_graphids:
            #use shape graphid to find corresponding shape groupref
            temp_groupref = shape_grouprefs[shape_graphids.index(x)]
            if not temp_groupref is None:
                #find nodes that correspond to the shape groupref
                matchGraphID = featureDFs['datanodeDF']['graphid'][featureDFs['datanodeDF']['groupref'] == temp_groupref].tolist()
                #find the textlabels of the matched datanodes
                temp_textlabel = [textlabels[node_graphids.index(element)] for element in matchGraphID]
                #print("Shape ID", x, temp_groupref, matchGraphID, temp_textlabel)
                return(temp_textlabel)
            elif temp_groupref is None and not shape_textlabels[shape_graphids.index(x)] is None:
                return(shape_textlabels[shape_graphids.index(x)])
            else:
                #Ungrouped shape ID
                return(x)
        elif x in label_graphids:
            #Use label graphid to find corresponding label groupref
            matchLabelID = featureDFs['labelDF']['graphid'][featureDFs['labelDF']['graphid'] == x].tolist()
            if not matchLabelID is None:
                #Find textlabels of these matched labels
                temp_textlabel = [label_textlabels[label_graphids.index(element)] for element in matchLabelID]
                return(temp_textlabel)
            else:
                #Ungrouped label ID
                return(x)
        elif x in group_graphids:
            #Use group graphId to find corresponding groupId
            temp_groupref = group_groupids[group_graphids.index(x)]
            if not temp_groupref is None:
                #Find all datanodes that are part of this group
                matchGraphID = featureDFs['datanodeDF']['graphid'][featureDFs['datanodeDF']['groupref'] == temp_groupref].tolist()
                #Find the textlabels of the matched datanodes
                temp_textlabel = [textlabels[node_graphids.index(element)] for element in matchGraphID]
                #print("Group ID: ", x, temp_groupref, matchGraphID, temp_textlabel)
                return(temp_textlabel)
            else:
                #print("Unmatched groupID: ", x)
                return(x)
        else:
            #print("unknown: ", x)
            return(x)
    else:
        #print("Unmapped node", x)
        return(x)

def replaceDist(seriesObj, textlabels):
    if isinstance(seriesObj, list) or isinstance(seriesObj, set):
        return(any([l2 in textlabels for l2 in seriesObj]))
    else:
        seriesObj = [seriesObj]
        return(any([l2 in textlabels for l2 in seriesObj]))

def processInteractDF(featureDFs, featureList):

    """
    Attempt to map node ids to datanodes.
    That is: if a matched node ID is a shape, group or label, check if any datanode has that shape, group or label attribute.
    If so, replace the matched node id with the datanode id.
    The goal of this function is to map every interaction endpoint to a datanode, as far as possible.
    """

    anchor_graphids = getAnchor(featureList[2], featureDFs)
    featureDFs['interactDF']['matched_Node1_textlabel'] = [matchRef(x, featureDFs, featureList, anchor_graphids) for x in featureDFs['interactDF']['matched_Node1'].tolist()]
    featureDFs['interactDF']['matched_Node2_textlabel'] = [matchRef(x,featureDFs, featureList, anchor_graphids) for x in featureDFs['interactDF']['matched_Node2'].tolist()]

    # if endpoints are successfully mapped to nodes, assign a distance of zero to that mapped node

    textlabels = featureDFs['datanodeDF']['textlabel'].tolist()
    featureDFs['interactDF']['matched_Node1_textlabel'] = [matchRef(x, featureDFs, featureList, anchor_graphids) for x in featureDFs['interactDF']['matched_Node1'].tolist()]
    featureDFs['interactDF']['matched_Node2_textlabel'] = [matchRef(x,featureDFs, featureList, anchor_graphids) for x in featureDFs['interactDF']['matched_Node2'].tolist()]

    featureDFs['interactDF']['matched_Node1_dist'] = np.where(featureDFs['interactDF']['matched_Node1_textlabel'].apply(replaceDist, textlabels = textlabels), 0, featureDFs['interactDF']['matched_Node1_dist'])
    featureDFs['interactDF']['matched_Node2_dist'] = np.where(featureDFs['interactDF']['matched_Node2_textlabel'].apply(replaceDist, textlabels = textlabels), 0, featureDFs['interactDF']['matched_Node2_dist'])

    featureDFs['interactDF']['matched_Node1_cat'] = np.where(featureDFs['interactDF']['matched_Node1_textlabel'].isnull(), "unmapped", featureDFs['interactDF']['matched_Node1_cat'])
    featureDFs['interactDF']['matched_Node2_cat'] = np.where(featureDFs['interactDF']['matched_Node2_textlabel'].isnull(), "unmapped", featureDFs['interactDF']['matched_Node2_cat'])
    return(featureDFs['interactDF'])

def flatten(l):
    
    #https://stackoverflow.com/a/2158532
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el

def checkArrow(arrow):
    #arrow_classes = ["Line", "Arrow", "Receptor", "ReceptorRound", "ReceptorSquare", "LigandRound", "LigandSquare", "TBar"]
    #WP vocabulatry: https://www.w3.org/2012/pyRdfa/extract?uri=http://vocabularies.wikipathways.org/wpTypes#
    activationArrow = ["Line", "Arrow", "mim-stimulation", "mim-transcription-translation", "mim-transcription", "mim-translation", "mim-conversion", "mim-catalysis", "mim-binding", "mim-branching-left", "mim-branching-right", "mim-modification", "mim-necessary-stimulation"]#, "group"]
    groupArrow = ["group"]
    inhibitionArrow = ["TBar", "mim-inhibition"]

    if not arrow is None:
        if arrow in activationArrow:
            interactionType= "a"
            signal="a"
        elif arrow in inhibitionArrow:
            interactionType= "i"
            signal="i"
        elif arrow in groupArrow:
            interactionType= "g"
            signal="a"
        else:
            interactionType = "u"
            #print(str(arrow))
            signal="u"
    else:
        arrow = "undefined"
        interactionType="u"
        signal="u"
    return([arrow, interactionType, signal])

def makeGraph(featureDFs, featureList, reWire_Inferred_Groups=False):
    """
        Prepare a networkX graph from the processed interaction dataframe
    """
    anchor_graphids = getAnchor(featureList[2], featureDFs)
    #Add interactions between nodes with textlabels
    #G = nx.from_pandas_edgelist(featureDFs['interactDF'], source = 'matched_Node1_textlabel', target = 'matched_Node2_textlabel', edge_attr='arrow')
    featureDFs['interactDF']['color'] = [''.join(["#",str(temp)]) for temp in featureDFs['interactDF']['color'].tolist()]
    G = nx.DiGraph()
    for node1, node2, arrow, edgeID, color in zip(featureDFs['interactDF']['matched_Node1_textlabel'].tolist(), featureDFs['interactDF']['matched_Node2_textlabel'].tolist(), featureDFs['interactDF']['arrow'].tolist(), featureDFs['interactDF']['edgeID'].tolist(), featureDFs['interactDF']['color'].tolist()):
        node1 = list(flatten([node1]))
        node2 = list(flatten([node2]))
        for n1 in node1:
            for n2 in node2:
                if not arrow is None:
                    arrow, interaction, signal = checkArrow(arrow)
                    G.add_edge(str(n1), str(n2), edgeID = str(edgeID), arrow = str(arrow), interaction = str(interaction), signal=str(signal), color=str(color))
                else:
                    G.add_edge(str(n1), str(n2), edgeID = str(edgeID), arrow = "unknown", interaction = "u", signal="u", color=str(color))

    #remove nodes that are 'None'
    if None in G.nodes():
        G.remove_node(None)
    else:
        if "None" in G.nodes():
            G.remove_node("None")
        else:
            G = G

    #Add nodes that are not in any interaction
    all_datanodes = featureDFs['datanodeDF']['textlabel'].tolist()
    for node in all_datanodes:
        if not node in G.nodes():
            G.add_node(str(node))

    # get list of unlabeled nodes
    unlabeledNodes = set(featureDFs['interactDF'][featureDFs['interactDF'].matched_Node1 == featureDFs['interactDF'].matched_Node1_textlabel]['matched_Node1'].tolist() + featureDFs['interactDF'][featureDFs['interactDF'].matched_Node2 == featureDFs['interactDF'].matched_Node2_textlabel]['matched_Node2'].tolist())
    unlabeledNodes = unlabeledNodes.intersection(set(list(G.nodes())))

    #remove unlabeled nodes
    for node in list(set(unlabeledNodes)):
        node1 = matchRef(node, featureDFs, featureList, anchor_graphids)
        # print(node, node1)
        G = nx.relabel_nodes(G, {node: node1})
        G = passThroughUnlabeled(node, G)
    
    #rewire groups
    node_grouprefs = featureDFs['datanodeDF']['groupref'].tolist()

    textlabels = featureDFs['datanodeDF']['textlabel'].tolist()

    for node1 in textlabels:
        group_node1 = node_grouprefs[textlabels.index(node1)]
        for node2 in textlabels:
            group_node2 = node_grouprefs[textlabels.index(node2)]
            if (not group_node1 is None) & (not group_node2 is None):
                if (group_node1 is group_node2) & (len(group_node1) > 0) & (len(group_node2) > 0) & (not node1 is node2):
                    #print("test",group_node1, group_node2, node1, node2)
                    G.add_edge(str(node1), str(node2), arrow="group")

    if reWire_Inferred_Groups:
        node_grouprefs = set(featureDFs['interactDF']['edgeID'].tolist())
        
        for group in node_grouprefs:
            if "unlabeled" in str(group):
                #find the nodes in the inferred group - concatenate the source and target
                groupNodes = featureDFs['interactDF']['matched_Node1_textlabel'][featureDFs['interactDF']['edgeID'] == group].tolist()
                groupNodes.extend(featureDFs['interactDF']['matched_Node2_textlabel'][featureDFs['interactDF']['edgeID'] == group].tolist())
                groupNodes = set(groupNodes)
                #now add edges from the nodes in that group to all the nodes that are downstream of any member of that group
                #and add edges from all the nodes that are upstream of any member of that group to all the nodes in that group
                for x in groupNodes:
                    upstreamNodes = set(flatten(featureDFs['interactDF']['matched_Node1_textlabel'][featureDFs['interactDF']['matched_Node2_textlabel'] == x]))
                    downstreamNodes = set(flatten(featureDFs['interactDF']['matched_Node1_textlabel'][featureDFs['interactDF']['matched_Node1_textlabel'] == x]))
                    for y in groupNodes:
                        for node1, node2, arrow, edgeID, color in zip(featureDFs['interactDF']['matched_Node1_textlabel'].tolist(), featureDFs['interactDF']['matched_Node2_textlabel'].tolist(), featureDFs['interactDF']['arrow'].tolist(), featureDFs['interactDF']['edgeID'].tolist(), featureDFs['interactDF']['color'].tolist()):
                            node1 = list(flatten([node1]))
                            node2 = list(flatten([node2]))
                            for n1 in node1:
                                for n2 in node2:                                
                                    if (n1 in upstreamNodes) & (x == n2) & (n1 != y):
                                        if not arrow is None:
                                            arrow, interaction, signal = checkArrow(arrow)
                                            G.add_edge(str(n1), str(y), edgeID = str(edgeID), arrow = str(arrow), interaction = str(interaction), signal=str(signal), color=str(color))
                                        else:
                                            G.add_edge(str(n1), str(y), edgeID = str(edgeID), arrow = "unknown", interaction = "u", signal="u", color=str(color))

                                    elif (n2 in downstreamNodes) & (x == n1) & (n2 != y):
                                        if not arrow is None:
                                            arrow, interaction, signal = checkArrow(arrow)
                                            G.add_edge(str(y), str(n2), edgeID = str(edgeID), arrow = str(arrow), interaction = str(interaction), signal=str(signal), color=str(color))
                                        else:
                                            #G.add_edge(n2, y, edgeID = edgeID, arrow = "unknown", interaction = "u", signal="u", color=color)
                                            G.add_edge(str(y), str(n2), edgeID = str(edgeID), arrow = str(arrow), interaction = str("u"), signal=str("u"), color=str(color))
                                                
    #last pass to check node aliases
    for node in G.nodes():
        node1 = matchRef(node, featureDFs, featureList, anchor_graphids)
        G = nx.relabel_nodes(G, {node: node1})
    
    #finally - remove newlines in graph nodes
    for node in G.nodes():
        #node1 = re.sub("\n", " ", node) # 
        node1 = node.strip()
        G = nx.relabel_nodes(G, {node: node1})
    
    #remove unlabeled nodes that are connected to nothing except themselves
    for node in unlabeledNodes:
        if node in G.nodes():
            #print(node)
            upstream = list(G.predecessors(node))
            downstream = list(G.successors(node))
            #print(upstream, downstream)
            if upstream == downstream:
                G.remove_node(node)
            if (node, node) in G.edges():
                G.remove_edge(node, node)
    
    #add node attributes from datanodeDF
    for attr in featureDFs['datanodeDF'].columns:
        attr = str(attr)
        attrDict = {}
        for tl in featureDFs['datanodeDF']["textlabel"].tolist():
            attrDict[tl] = str(featureDFs['datanodeDF'][attr][featureDFs['datanodeDF']["textlabel"] == tl].tolist()[0])
        for gi in featureDFs['datanodeDF']["graphid"].tolist():
            attrDict[gi] = str(featureDFs['datanodeDF'][attr][featureDFs['datanodeDF']["graphid"] == gi].tolist()[0])
        nx.set_node_attributes(G, name= attr, values=attrDict)


    return(G)


def addEdgeAnnotations(annotationFile, graph):
    """
    Go through all edges in the edgelist and add colors from a pre-made annotation file
    """
    #Read annotation file
    annotFile = pd.read_csv(annotationFile)
    #print(annotFile.head(10))

    #Find edges that map to labels
    for edge in G.edges():
        print(G[edge[0]][edge[1]]['color'])

def mapUnlabeledGroups(featureDFs):

    points = [[x, y] for x, y in zip(featureDFs['datanodeDF']['centerX'], featureDFs['datanodeDF']['centerY'])]
    group_node_data=[]
    distDF = distance_matrix(points, points)
    nodelist = featureDFs['datanodeDF']['graphid'].tolist()
    textlabels = featureDFs['datanodeDF']['textlabel'].tolist()

    for node1 in range(0, len(nodelist)):
        graphid1 = nodelist[node1]
        height1 = featureDFs['datanodeDF']['height'].tolist()[node1]
        width1 = featureDFs['datanodeDF']['width'].tolist()[node1]
        for node2 in range(0, len(nodelist)):
            if not node1 is node2:
                graphid2 = nodelist[node2]
                height2 = featureDFs['datanodeDF']['height'].tolist()[node2]
                width2 = featureDFs['datanodeDF']['width'].tolist()[node2]
                dist = distDF[node1, node2]
                cond1 = dist <= ((height1 + height2)/2) + ((height1 + height2)/4)
                #cond2a = dist <= (width1 + width2)/2 + (width1 + width2)/4 #close together horizontally 
                #cond2b = abs(featureDFs['datanodeDF']['centerX'].tolist()[node1] - featureDFs['datanodeDF']['centerX'].tolist()[node2]) <= 0.5*abs(featureDFs['datanodeDF']['centerY'].tolist()[node1] - featureDFs['datanodeDF']['centerY'].tolist()[node2]) #linear/parallel
                #cond2 = cond2a and cond2b
                if cond1: #or cond2:
                    #print("group", height1, height2, width1, width2, dist, textlabels[node1], textlabels[node2])
                    graphid = "_".join([str(graphid1), str(graphid2)]) + "_unlabeled"
                    groupRef = "_".join([str(graphid1), str(graphid2)]) + "_unlabeled"
                    arrow = "group"
                    node1_x = featureDFs['datanodeDF']['centerX'].tolist()[node1]
                    node1_y = featureDFs['datanodeDF']['centerY'].tolist()[node1]
                    node2_x = featureDFs['datanodeDF']['centerX'].tolist()[node2]
                    node2_y = featureDFs['datanodeDF']['centerY'].tolist()[node2]
                    node1_coords = Point(node1_x, node1_y)
                    node2_coords = Point(node2_x, node2_y)
                    geometry = MultiPoint([(node1_x, node1_y), (node2_x, node2_y)]).convex_hull
                    #group_node_data.append([graphid, groupRef, graphid1, node1_x, node1_y, node1_coords, graphid2, node2_x, node2_y, node2_coords, arrow, geometry])
                    group_node_data.append([graphid, groupRef, graphid1, node1_x, node1_y, node1_coords, graphid2, node2_x, node2_y, node2_coords, arrow, geometry, "unlabeled", "unlabeled", None])

    group_nodeDF = gpd.GeoDataFrame(pd.DataFrame(group_node_data, columns=["edgeID", "edgeGroup", "node1_graphref", "node1_x", "node1_y", "node1_coords", "node2_graphref", "node2_x", "node2_y", "node2_coords", "arrow", "geometry", "comments", "color", "anchor"]))

    #group_nodeDF = gpd.GeoDataFrame(pd.DataFrame(group_node_data, columns=["edgeID", "edgeGroup", "node1_graphref", "node1_x", "node1_y", "node1_coords", "node2_graphref", "node2_x", "node2_y", "node2_coords", "arrow", "geometry"]))
    featureDFs['interactDF'] = pd.concat([featureDFs['interactDF'], group_nodeDF], sort=True)
    featureDFs['interactDF'].reset_index(inplace=True, drop=True)
    return(featureDFs['interactDF'])

def passThroughUnlabeled(node, graph):
    """
        Given a node, remove that node and pass signal directly from all upstream nodes to all downstream nodes
    """
    upstream = graph.predecessors(node)
    downstream = graph.successors(node)
    new_edges = [(p,s) for p in upstream for s in downstream]
    graph.remove_node(node)
    graph.add_edges_from(new_edges)
    return(graph)

def edgeListToCSV(G, pathwayID):
    fh=open('_'.join([pathwayID, "edgeList.csv"]),'w')
    fh.write(",".join(["Source", "Target", "edgeID", "signal", "color"]))
    fh.write("\n")
    for edge in G.edges(data=True):
        #print(edge)
        fh.write(",".join(["\""+str(edge[0])+"\"", "\""+str(edge[1])+"\""]))#, edge[2]['edgeID'], edge[2]['interaction'], edge[2]['arrow'], edge[2]['signal']]))
        for k in ['edgeID', 'signal', 'color']: #edge[2].keys():
            if k in edge[2].keys():
                #print(k)
                fh.write(",")
                fh.write(str("\""+str(edge[2][k])+"\""))
            else:
                fh.write(",")
                fh.write("\""+"None"+"\"")
        fh.write("\n")
    fh.close()

def runParsePathway(s, pathwayID):
    #set up processing pipeline, stage processed input
    gpml = givenCodeGetGPML(s, pathwayID.encode('utf-8'))
    featureList = makeFeatureLists(gpml)
    featureDFs = getFeatureDFs(featureList)
    #pipeline
    featureDFs['interactDF'] = mapEndPoints(featureDFs)
    featureDFs['interactDF'] = mapPointsToInteractions(featureDFs, featureList)
    featureDFs['interactDF'] = processInteractDF(featureDFs, featureList)
    #tests
    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #    print("TEST2", featureDFs['interactDF']['matched_Node1'][featureDFs['interactDF']['edgeID'] == "e601a"], featureDFs['interactDF']['matched_Node2'][featureDFs['interactDF']['edgeID'] == "e601a"])
    #    print("***")
    #make diagnostic plots
    makeMatchedDistHist(featureDFs, pathwayID)
    makeEndPointBarPlot(featureDFs, pathwayID)
    #write out graphml
    graph = makeGraph(featureDFs, featureList)
    featureDFs['interactDF'].columns = ['Color' if x=='color' else x for x in featureDFs['interactDF'].columns]
    #print(set(list(graph.nodes)))
    #write out edgelist
    featureDFs['interactDF'].to_csv('_'.join([pathwayID, "interactDF.csv"]))
    featureDFs['datanodeDF'].to_csv('_'.join([pathwayID, "datanodeDF.csv"]))
    edgeListToCSV(graph, pathwayID)
    #Write out graphml
    nx.write_graphml_lxml(graph, ''.join([pathwayID, "_graph.graphml"]))
    return graph

def unitTest(name, pathwayID):
    gpml = open(str(name), 'r') #open("WP23_minimal.gpml", "r")
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

def getNetworkDiffs(manualOut, programOut):
    manual = pd.read_csv(manualOut, sep=None, engine = "python", header=0)
    program = pd.read_csv(programOut, sep=None, engine = "python", header=0)
    print(manualOut)
    print(programOut)
    program=program.loc[:,["Source", "Target", "signal"]]
    #print(manual.columns)
    manual=manual.loc[:,["Source", "Target", "signal"]]

    manualNodes = set(set(manual.Source).union(set(manual.Target)))
    programNodes = set(set(program.Source).union(set(program.Target)))
    #Make edges from manual and program so that we can compare node1, node2, and the annotated signals.
    manList = []
    progList = []
    for i in range(len(manual)):
        manList.append(" ".join([str(list(manual.Source)[i]), str(list(manual.Target)[i]), str(list(manual.signal)[i])]))
    
    for i in range(len(program)):
        if str(list(program.Source)[i]) in manualNodes and str(list(program.Target)[i]) in manualNodes:
            progList.append(" ".join([str(list(program.Source)[i]), str(list(program.Target)[i]), str(list(program.signal)[i])]))
    

    extraInManual = manualNodes.difference(programNodes)
    extraInProgram = programNodes.difference(manualNodes)
    edgesInManual = len(manual)
    edgesInProgram = len(program)
    extraEdgesInProgram = set(progList).difference(set(manList))
    extraEdgesInManual = set(manList).difference(set(progList))
    propTrueNodes = float(len(manualNodes.intersection(programNodes)))/len(manualNodes) #/len(programNodes) #Proportion of true nodes = number(intersection of manual and program nodes)/number(program nodes)
    propTrueEdges = float(len(set(manList).intersection(set(progList))))/len(manList)#/len(progList) #Proportion of true edges = number(intersection of manual and program edges)/number(program edges)
    truePositiveEdges = float(len(set(manList).intersection(set(progList))))
    recall = float(truePositiveEdges)/len(manList)#(truePositiveEdges + len(extraEdgesInManual))
    #print(truePositiveEdges)
    #print(len(extraEdgesInProgram))
    precision = truePositiveEdges/len(progList) #(truePositiveEdges + len(extraEdgesInProgram))
    f1Score = float((2* precision * recall))/(precision + recall)
    #false positive = edge where there is no edge, ie, the edge is in program but not in manual, ie, extra edges in program = "number_of_extra_programEdges"
    #true positive = edge where there is edge in both program and in manual
    #false negative = missing edge in program, ie, extra edges in manual, = "number_of_extra_manualEdges"
    #true negative = correctly predicted no edge when there is no edge, ((n-1)*n)/2 where n = number of nodes
    #True positive rate = true positives/(true positives + false negatives)
    #False positive rate = false positives/(true negatives + false positives)
    trueNegatives = ((len(manualNodes)*(len(manualNodes) - 1))/2.0) - float(len(set(manList).intersection(set(progList))))
    truePositiveRate = truePositiveEdges/(truePositiveEdges + len(extraEdgesInManual))
    if (trueNegatives + float(len(extraEdgesInProgram))) > 0:
        falsePositiveRate = float(len(extraEdgesInProgram))/(trueNegatives + float(len(extraEdgesInProgram)))
    else:
        falsePositiveRate = 0.0
    returnVec = [len(manualNodes),
                len(programNodes),
                #"; ".join([str(temp) for temp in manualNodes]),
                #"; ".join([str(temp) for temp in programNodes]),
                "; ".join([str(temp) for temp in extraInManual]),
                "; ".join([str(temp) for temp in extraInProgram]),
                len(extraInManual),
                len(extraInProgram),
                len(extraEdgesInManual),
                len(extraEdgesInProgram),
                "; ".join(extraEdgesInManual),
                "; ".join(extraEdgesInProgram),
                propTrueNodes,
                propTrueEdges,
                truePositiveEdges,
                recall,
                precision,
                truePositiveRate,
                falsePositiveRate,
                f1Score]
    return(returnVec)

def evaluateProgramOutput():
    programOuts = ["WP49_edgeList.csv", "WP127_edgeList.csv", "WP195_edgeList.csv", "WP205_edgeList.csv", "WP231_edgeList.csv", "WP286_edgeList.csv", "WP23_edgeList.csv", "WP3935_edgeList.csv", "WP4482_edgeList.csv", "WP4868_edgeList.csv"] #, "test_1_edgeList.csv"
    manualOuts = ["WP49_edgeList_edited.csv", "WP127_edgeList_edited.csv", "WP195_edgeList_edited.csv", "WP205_edgeList_edited.csv",  "WP231_edgeList_edited.csv",  "WP286_edgeList_edited.csv", "WP23_edgeList_edited.csv", "WP3935_edgelist_edited.csv", "WP4482_edgelist_edited.csv", "WP4868_edgelist_edited.csv"] #, "test_1_edgeList_edited.csv"
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
    programOuts = ["WP49_edgeList.csv", "WP127_edgeList.csv", "WP195_edgeList.csv", "WP205_edgeList.csv", "WP231_edgeList.csv", "WP286_edgeList.csv", "WP23_edgeList.csv", "WP3935_edgeList.csv", "WP4868_edgeList.csv", "WP4482_edgeList.csv"]
    manualOuts = ["WP49_cytoscape.csv", "WP127_cytoscape.csv", "WP195_cytoscape.csv", "WP205_cytoscape.csv", "WP231_cytoscape.csv", "WP286_cytoscape.csv", "WP23_cytoscape.csv", "WP3935_cytoscape.csv", "WP4868_cytoscape.csv", "WP4482_cytoscape.csv"] 
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

def evaluateCytoscapeOutput():
    programOuts = ["WP49_cytoscape.graphml", "WP127_cytoscape.graphml", "WP195_cytoscape.graphml", "WP205_cytoscape.graphml", "WP231_cytoscape.graphml", "WP286_cytoscape.graphml", "WP23_cytoscape.graphml", "WP3935_cytoscape.graphml", "WP4868_cytoscape.graphml", "WP4482_cytoscape.graphml"]
    #Generate edgelists from cytoscape graphmls
    for i in range(len(programOuts)):
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
        
        nx.write_edgelist(graph, programOuts[i], delimiter="\t", data = ["targetArrows"]) # 
        print(programOuts[i])
        tempDF = pd.read_csv(programOuts[i], header=None, sep="\t", escapechar = "\\")
        tempDF.columns = ["Source", "Target", "signal"]
        tempDF.signal = [checkArrow(arrow)[1] for arrow in tempDF["signal"]]
        tempDF.to_csv(programOuts[i], quoting=csv.QUOTE_ALL)
        #print(tempDF)
    
    manualOuts = ["WP49_edgeList_edited.csv", "WP127_edgeList_edited.csv", "WP195_edgeList_edited.csv", "WP205_edgeList_edited.csv",  "WP231_edgeList_edited.csv",  "WP286_edgeList_edited.csv", "WP23_edgeList_edited.csv", "WP3935_edgelist_edited.csv", "WP4868_edgelist_edited.csv", "WP4482_edgelist_edited.csv"]
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
    programOuts = ["WP49_cytoscape.graphml", "WP127_cytoscape.graphml", "WP195_cytoscape.graphml", "WP205_cytoscape.graphml", "WP231_cytoscape.graphml", "WP286_cytoscape.graphml", "WP23_cytoscape.graphml"]
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
        tempDF = pd.read_csv(programOuts[i], header=None, sep="\t")
        tempDF.columns = ["Source", "Target", "signal"]
        tempDF.signal = [checkArrow(arrow)[1] for arrow in tempDF["signal"]]
        tempDF.to_csv(programOuts[i], quoting=csv.QUOTE_ALL)
        #print(tempDF)
    
    manualOuts = ["WP49_edgeList_edited.csv", "WP127_edgeList_edited.csv", "WP195_edgeList_edited.csv", "WP205_edgeList_edited.csv",  "WP231_edgeList_edited.csv",  "WP286_edgeList_edited.csv", "WP23_edgeList_edited.csv"]
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
    #cytoscapeFairplay = pd.read_csv("netDiffResult_cytoscape_fairplay.tsv", sep= "\t", index_col=0)
    #cytoscapeFairplay["Method"] = "cytoscape_fairplay"
    wikinetworks = pd.read_csv("netDiffResult_wikinetworks.tsv", sep= "\t", index_col=0)
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

