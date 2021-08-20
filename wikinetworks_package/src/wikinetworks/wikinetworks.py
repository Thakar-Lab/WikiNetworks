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
    
