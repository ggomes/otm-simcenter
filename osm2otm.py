from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon

import time
import numpy as np
import geopandas as gpd
import networkx as nx
from shapely.geometry import Point
import requests
from lxml import etree
import math
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import pickle

MilesToKilometers = 1.609344

road_param_types = {
    'residential':      {'id': 0, 'capacity': 2000, 'speed': 100, 'jam_density': 100},
    'tertiary':         {'id': 1, 'capacity': 2001, 'speed': 101, 'jam_density': 101},
    'tertiary_link':    {'id': 2, 'capacity': 2002, 'speed': 102, 'jam_density': 102},
    'unclassified':     {'id': 3, 'capacity': 2003, 'speed': 103, 'jam_density': 103},
    'secondary':        {'id': 4, 'capacity': 2004, 'speed': 104, 'jam_density': 104},
}

default_lanes_each_direction = {
    'residential': 1,
    'tertiary': 1,
    'secondary': 1,
    'primary': 1,
    'unclassified': 1
}

def new_link_id():
    global max_link_id
    max_link_id = max_link_id + 1
    return max_link_id

def adhoc_correct_discrepancy(link):

    if link['id'] == 415803770:
        # turns: left | | | lanes 3
        # OSM: https://www.openstreetmap.org/way/415803770#map=19/37.87693/-122.28277
        # Google: https://www.google.com/maps/@37.8768752,-122.2828014,76m/data=!3m1!1e3
        # Correction: Turn includes the parking lane. Ignore it
        link['turn_lanes'] = 'left||'

    if link['id'] == 415876791:
        # turns: left | through | through lanes 2
        # OSM: https://www.openstreetmap.org/way/415876791#map=19/37.87453/-122.26411
        # Google: https://www.google.com/maps/@37.8745003,-122.2648584,78a,35y,94.72h,32.94t/data=!3m1!1e3
        # Correction: Missing lane
        link['lanes'] = 3
        link['lanes_backward'] = 2
        link['turn_lanes_backward'] = '|'

    if link['id'] == 574381942:
        # Hearst @ LeConte
        # turns: left | through | through lanes 2
        # OSM: https://www.openstreetmap.org/way/574381942#map=19/37.87458/-122.26374
        # Google: https://www.google.com/maps/@37.8745017,-122.2643212,88m/data=!3m1!1e3
        # Correction: lanes -> 3
        link['lanes'] = 3
        link['lanes_backward'] = 2
        link['turn_lanes_backward'] = '|'

def project_geometry(geometry, crs=None, to_crs=None, to_latlong=False):
    """
    Project a shapely Polygon or MultiPolygon from lat-long to UTM, or
    vice-versa

    Parameters
    ----------
    geometry : shapely Polygon or MultiPolygon
        the geometry to project
    crs : dict
        the starting coordinate reference system of the passed-in geometry,
        default value (None) will set settings.default_crs as the CRS
    to_crs : dict
        if not None, just project to this CRS instead of to UTM
    to_latlong : bool
        if True, project from crs to lat-long, if False, project from crs to
        local UTM zone

    Returns
    -------
    tuple
        (geometry_proj, crs), the projected shapely geometry and the crs of the
        projected geometry
    """
    default_crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

    if crs is None:
        crs = default_crs

    gdf = gpd.GeoDataFrame()
    gdf.crs = crs
    gdf.gdf_name = 'geometry to project'
    gdf['geometry'] = None
    gdf.loc[0, 'geometry'] = geometry
    gdf_proj = project_gdf(gdf, to_crs=to_crs, to_latlong=to_latlong)
    geometry_proj = gdf_proj['geometry'].iloc[0]
    return geometry_proj, gdf_proj.crs

def project_gdf(gdf, to_crs=None, to_latlong=False):
    """
    Project a GeoDataFrame to the UTM zone appropriate for its geometries'
    centroid.

    The simple calculation in this function works well for most latitudes, but
    won't work for some far northern locations like Svalbard and parts of far
    northern Norway.

    Parameters
    ----------
    gdf : GeoDataFrame
        the gdf to be projected
    to_crs : dict
        if not None, just project to this CRS instead of to UTM
    to_latlong : bool
        if True, projects to latlong instead of to UTM

    Returns
    -------
    GeoDataFrame
    """

    default_crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

    assert len(gdf) > 0, 'You cannot project an empty GeoDataFrame.'
    start_time = time.time()

    # if gdf has no gdf_name attribute, create one now
    if not hasattr(gdf, 'gdf_name'):
        gdf.gdf_name = 'unnamed'

    # if to_crs was passed-in, use this value to project the gdf
    if to_crs is not None:
        projected_gdf = gdf.to_crs(to_crs)

    # if to_crs was not passed-in, calculate the centroid of the geometry to
    # determine UTM zone
    else:
        if to_latlong:
            # if to_latlong is True, project the gdf to latlong
            latlong_crs = default_crs
            projected_gdf = gdf.to_crs(latlong_crs)
            # log('Projected the GeoDataFrame "{}" to default_crs in {:,.2f} seconds'.format(gdf.gdf_name, time.time()-start_time))
        else:
            # else, project the gdf to UTM
            # if GeoDataFrame is already in UTM, just return it
            if (gdf.crs is not None) and ('+proj=utm ' in gdf.crs):
                return gdf

            # calculate the centroid of the union of all the geometries in the
            # GeoDataFrame
            avg_longitude = gdf['geometry'].unary_union.centroid.x

            # calculate the UTM zone from this avg longitude and define the UTM
            # CRS to project
            utm_zone = int(math.floor((avg_longitude + 180) / 6.) + 1)
            utm_crs = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(utm_zone)

            # project the GeoDataFrame to the UTM CRS
            projected_gdf = gdf.to_crs(utm_crs)
            # log('Projected the GeoDataFrame "{}" to UTM-{} in {:,.2f} seconds'.format(gdf.gdf_name, utm_zone, time.time()-start_time))

    projected_gdf.gdf_name = gdf.gdf_name
    return projected_gdf

def project_graph(G, to_crs=None):
    """
    Project a graph from lat-long to the UTM zone appropriate for its geographic
    location.

    Parameters
    ----------
    G : networkx multidigraph
        the networkx graph to be projected
    to_crs : dict
        if not None, just project to this CRS instead of to UTM

    Returns
    -------
    networkx multidigraph
    """

    G_proj = G.copy()
    start_time = time.time()

    # create a GeoDataFrame of the nodes, name it, convert osmid to str
    nodes, data = zip(*G_proj.nodes(data=True))
    gdf_nodes = gpd.GeoDataFrame(list(data), index=nodes)
    gdf_nodes.crs = G_proj.graph['crs']
    gdf_nodes.gdf_name = '{}_nodes'.format(G_proj.name)

    # create new lat/lon columns just to save that data for later, and create a
    # geometry column from x/y
    gdf_nodes['lon'] = gdf_nodes['x']
    gdf_nodes['lat'] = gdf_nodes['y']
    gdf_nodes['geometry'] = gdf_nodes.apply(lambda row: Point(row['x'], row['y']), axis=1)
    # log('Created a GeoDataFrame from graph in {:,.2f} seconds'.format(time.time()-start_time))

    # project the nodes GeoDataFrame to UTM
    gdf_nodes_utm = project_gdf(gdf_nodes, to_crs=to_crs)

    # extract data for all edges that have geometry attribute
    edges_with_geom = []
    for u, v, key, data in G_proj.edges(keys=True, data=True):
        if 'geometry' in data:
            edges_with_geom.append({'u':u, 'v':v, 'key':key, 'geometry':data['geometry']})

    # create an edges GeoDataFrame and project to UTM, if there were any edges
    # with a geometry attribute. geom attr only exists if graph has been
    # simplified, otherwise you don't have to project anything for the edges
    # because the nodes still contain all spatial data
    if len(edges_with_geom) > 0:
        gdf_edges = gpd.GeoDataFrame(edges_with_geom)
        gdf_edges.crs = G_proj.graph['crs']
        gdf_edges.gdf_name = '{}_edges'.format(G_proj.name)
        gdf_edges_utm = project_gdf(gdf_edges, to_crs=to_crs)

    # extract projected x and y values from the nodes' geometry column
    start_time = time.time()
    gdf_nodes_utm['x'] = gdf_nodes_utm['geometry'].map(lambda point: point.x)
    gdf_nodes_utm['y'] = gdf_nodes_utm['geometry'].map(lambda point: point.y)
    gdf_nodes_utm = gdf_nodes_utm.drop('geometry', axis=1)
    # log('Extracted projected node geometries from GeoDataFrame in {:,.2f} seconds'.format(time.time()-start_time))

    # clear the graph to make it a blank slate for the projected data
    start_time = time.time()
    edges = list(G_proj.edges(keys=True, data=True))
    graph_name = G_proj.graph['name']
    G_proj.clear()

    # add the projected nodes and all their attributes to the graph
    G_proj.add_nodes_from(gdf_nodes_utm.index)
    attributes = gdf_nodes_utm.to_dict()
    for label in gdf_nodes_utm.columns:
        nx.set_node_attributes(G_proj, name=label, values=attributes[label])

    # add the edges and all their attributes (including reconstructed geometry,
    # when it exists) to the graph
    for u, v, key, attributes in edges:
        if 'geometry' in attributes:
            row = gdf_edges_utm[(gdf_edges_utm['u']==u) & (gdf_edges_utm['v']==v) & (gdf_edges_utm['key']==key)]
            attributes['geometry'] = row['geometry'].iloc[0]

        # attributes dict contains key, so we don't need to explicitly pass it here
        G_proj.add_edge(u, v, **attributes)

    # set the graph's CRS attribute to the new, projected CRS and return the
    # projected graph
    G_proj.graph['crs'] = gdf_nodes_utm.crs
    G_proj.graph['name'] = '{}_UTM'.format(graph_name)
    if 'streets_per_node' in G.graph:
        G_proj.graph['streets_per_node'] = G.graph['streets_per_node']
    # log('Rebuilt projected graph in {:,.2f} seconds'.format(time.time()-start_time))
    return G_proj

def overpass_request(data, timeout=180):
    url = 'http://overpass-api.de/api/interpreter'
    headers = {'User-Agent': 'Python OSMnx package (https://github.com/gboeing/osmnx)', 'Accept-Encoding': 'gzip, deflate', 'Accept': '*/*', 'Connection': 'keep-alive', 'referer': 'Python OSMnx package (https://github.com/gboeing/osmnx)', 'Accept-Language': 'en'}
    response = requests.post(url, data=data, timeout=timeout, headers=headers)
    response_json = response.json()
    return response_json
#
# def query_json():
#
#     # Dwight and 6th
#     # west = -122.29541
#     # east = -122.29246
#     # south = 37.8596
#     # north = 37.86111
#
#     west = -122.2981
#     north = 37.8790
#     east = -122.2547
#     south = 37.8594
#
#     infrastructure = 'way["highway"]'
#     timeout = 180
#     osm_filter = '["area"!~"yes"]["highway"!~"cycleway|footway|path|pedestrian|steps|track|corridor|proposed|construction|bridleway|abandoned|platform|raceway|service"]["motor_vehicle"!~"no"]["motorcar"!~"no"]["access"!~"private"]["service"!~"parking|parking_aisle|driveway|private|emergency_access"]'
#     maxsize = ''
#
#     # turn bbox into a polygon and project to local UTM
#     polygon = Polygon([(west, south), (east, south), (east, north), (west, north)])
#     geometry, crs_proj = project_geometry(polygon)
#
#     if isinstance(geometry, Polygon):
#         geometry_proj_consolidated_subdivided = MultiPolygon([geometry])
#
#     geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
#
#     response_jsons = []
#
#     for poly in geometry:
#
#         west, south, east, north = poly.bounds
#         query_template = '[out:json][timeout:{timeout}]{maxsize};({infrastructure}{filters}({south:.6f},{west:.6f},{north:.6f},{east:.6f});>;);out;'
#         query_str = query_template.format(north=north, south=south,
#                                           east=east, west=west,
#                                           infrastructure=infrastructure,
#                                           filters=osm_filter,
#                                           timeout=timeout, maxsize=maxsize)
#         response_json = overpass_request(data={'data': query_str}, timeout=timeout)
#         response_jsons.append(response_json)
#
#     return response_jsons

def read_way(element):
    # OSM tags considered:
    #     'name'
    #     'highway',
    #     'lanes',
    #     'lanes:forward',
    #     'lanes:backward',
    #     'lanes:both_ways',
    #     'turn',
    #     'turn:forward',
    #     'turn:backward',
    #     'turn:both_ways',
    #     'turn:lanes',
    #     'turn:lanes:forward',
    #     'turn:lanes:backward',
    #     'turn:lanes:both_ways',
    #     'maxspeed',
    #     'name',
    #     'oneway'

    if 'tags' not in element:
        print("ERROR: LINK WITH NO TAGS")
        return None

    tags = element['tags']

    nodes = element['nodes']
    link = {
        'id': element['id'],
        'length': 0,
        'start_node_id': nodes[0],
        'end_node_id': nodes[-1],
        'roadparam': 0,
        'nodes': nodes
    }

    # name ..................
    if 'name' in tags:
        link['name'] = tags['name']
    else:
        link['name'] = ''

    # highway ..................
    if 'highway' in tags:
        link['highway'] = tags['highway']
    else:
        link['highway'] = 'default_highway'

    # oneway ..................
    if 'oneway' in tags:
        if tags['oneway'] == 'yes':
            link['flip'] = False
            link['bidirectional'] = False

        if tags['oneway'] == 'no':
            link['flip'] = False
            link['bidirectional'] = True

        if tags['oneway'] == '-1':
            link['flip'] = True
            link['bidirectional'] = False
    else:
        link['flip'] = False
        if 'lanes:backward' in tags:
            link['bidirectional'] = True
        else:
            link['bidirectional'] = False

    # max speed ..................
    if 'maxspeed' in tags:
        x = tags['maxspeed'].split(" ")

        if len(x)<2:
            link['maxspeed_kph'] = float(x[0])*MilesToKilometers
        else:
            if x[1].lower()=='kph':
                link['maxspeed_kph'] = float(x[0])
            elif x[1].lower()=='mph':
                link['maxspeed_kph'] = float(x[0])*MilesToKilometers
            else:
                print("ERROR UNKNOWN UNITS")
    else:
        link['maxspeed_kph'] = 50

    # lanes, lanes_backward ..........................
    if link['bidirectional']:

        # (T,T,T),(F,T,T)
        if 'lanes:forward' in tags and 'lanes:backward' in tags:
            link['lanes'] = int(tags['lanes:forward'])
            link['lanes_backward'] = int(tags['lanes:backward'])

        # (T,T,F)
        elif 'lanes' in tags and 'lanes:forward' in tags:
            link['lanes'] = int(tags['lanes:forward'])
            link['lanes_backward'] = int(tags['lanes']) - int(tags['lanes:forward'])

        # (T,F,T)
        elif 'lanes' in tags and 'lanes:backward' in tags:
            link['lanes_backward'] = int(tags['lanes:backward'])
            link['lanes'] = int(tags['lanes']) - int(tags['lanes:backward'])

        # (T,F,F)
        elif 'lanes' in tags:

            # stupid, but see https://wiki.openstreetmap.org/wiki/Key:lanes#Assumptions
            if (int(tags['lanes']) % 2) == 0:
                link['lanes'] = int(int(tags['lanes']) / 2)
            else:
                link['lanes'] = int(tags['lanes'])
            link['lanes_backward'] = link['lanes']

        # (F,T,F)
        elif 'lanes:forward' in tags:
            link['lanes'] = int(tags['lanes:forward'])
            link['lanes_backward'] = link['lanes']

        # (F,F,T)
        elif 'lanes:backward' in tags:
            link['lanes'] = int(tags['lanes:backward'])
            link['lanes_backward'] = link['lanes']

        # (F,F,F)
        else:
            link['lanes'] = np.nan
            link['lanes_backward'] = np.nan

    else:
        link['lanes_backward'] = 0
        if 'lanes' in tags:
            link['lanes'] = int(tags['lanes'])
        elif 'lanes:forward' in tags:
            link['lanes'] = int(tags['lanes:forward'])
        else:
            link['lanes'] = np.nan

    # turns_lanes, turns_lanes_backward ................

    # check no clash in turn, turn:forward, turn:backward, turn:both_ways
    if 'turn' in tags and 'turn:forward' in tags:
        print("ERROR 'turn' in tags and 'turn:forward' in tags")

    if 'turn' in tags and 'turn:backward' in tags:
        print("ERROR 'turn' in tags and 'turn:backward' in tags")

    if 'turn' in tags and 'turn:both_ways' in tags:
        print("ERROR 'turn' in tags and 'turn:both_ways' in tags")

    if 'turn:both_ways' in tags and 'turn:forward' in tags:
        print("ERROR 'turn:both_ways' in tags and 'turn:forward' in tags")

    if 'turn:both_ways' in tags and 'turn:backward' in tags:
        print("ERROR 'turn:both_ways' in tags and 'turn:backward' in tags")

    # if no clash, copy turn and turn:both_ways to turn:forward and turn:backward
    if 'turn' in tags:
        tags['turn:forward'] = tags['turn']
        tags['turn:backward'] = tags['turn']

    if 'turn:both_ways' in tags:
        tags['turn:forward'] = tags['turn:both_ways']
        tags['turn:backward'] = tags['turn:both_ways']

    # check no clash turn:lanes, turn:lanes:forward, turn:lanes:backward, turn:lanes:both_ways
    if 'turn:lanes' in tags and 'turn:lanes:forward' in tags:
        print("ERROR 'turn:lanes' in tags and 'turn:lanes:forward' in tags")

    if 'turn:lanes' in tags and 'turn:lanes:backward' in tags:
        print("ERROR 'turn:lanes' in tags and 'turn:lanes:backward' in tags")

    if 'turn:lanes' in tags and 'turn:lanes:both_ways' in tags:
        print("ERROR 'turn:lanes' in tags and 'turn:lanes:both_ways' in tags")

    if 'turn:lanes:both_ways' in tags and 'turn:lanes:forward' in tags:
        print("ERROR 'turn:lanes:both_ways' in tags and 'turn:lanes:forward' in tags")

    if 'turn:lanes:both_ways' in tags and 'turn:lanes:backward' in tags:
        print("ERROR 'turn:lanes:both_ways' in tags and 'turn:lanes:backward' in tags")

    # if no clash, copy turn and turn:lanes:both_ways to turn:lanes:forward and turn:lanes:backward
    if 'turn:lanes' in tags:
        tags['turn:lanes:forward'] = tags['turn:lanes']
        tags['turn:lanes:backward'] = tags['turn:lanes']

    if 'turn:lanes:both_ways' in tags:
        tags['turn:lanes:forward'] = tags['turn:lanes:both_ways']
        tags['turn:lanes:backward'] = tags['turn:lanes:both_ways']

    # check clash turn:lanes:* and turn:*
    if 'turn:lanes:forward' in tags and 'turn:forward' in tags:
        if tags['turn:lanes:forward'] != tags['turn:forward']:
            print("ERROR 'turn:lanes:forward' in tags and 'turn:forward' in tags")

    if 'turn:lanes:backward' in tags and 'turn:backward' in tags:
        if tags['turn:lanes:backward'] != tags['turn:backward']:
            print("ERROR 'turn:lanes:backward' in tags and 'turn:backward' in tags")

    # copy turn to turn lanes
    if 'turn:forward' in tags:
        if math.isnan(link['lanes']):
            tags['turn:lanes:forward'] = '*'  # undefined
        else:
            tags['turn:lanes:forward'] = "|".join([tags['turn:forward'] for i in range(0,link['lanes'])])

    if 'turn:backward' in tags:
        if math.isnan(link['lanes_backward']):
            tags['turn:lanes:backward'] = '*'  # undefined
        else:
            tags['turn:lanes:backward'] = "|".join([tags['turn:backward'] for i in range(0, link['lanes_backward'])])

    # now all of the information is in tags['turn:lanes:forward'] and tags['turn:lanes:backward']

    # copy to link
    if 'turn:lanes:forward' in tags:
        link['turn_lanes'] = tags['turn:lanes:forward']
    else:
        if math.isnan(link['lanes']):
            link['turn_lanes'] = '*'  # undefined
        else:
            link['turn_lanes'] = '|'*(link['lanes']-1)

    if link['bidirectional'] and 'turn:lanes:backward' in tags:
        link['turn_lanes_backward'] = tags['turn:lanes:backward']
    else:
        if math.isnan(link['lanes_backward']):
            link['turn_lanes_backward'] = '*'  # undefined
        else:
            link['turn_lanes_backward'] = '|'*(link['lanes_backward']-1)

    # resolve undefined number of lanes
    if math.isnan(link['lanes']):
        if link['turn_lanes'] == '*':
            link['lanes'] = default_lanes_each_direction[link['highway']]
        else:
            link['lanes'] = len(link['turn_lanes'].split('|'))

    if math.isnan(link['lanes_backward']):
        if link['turn_lanes_backward'] == '*':
            link['lanes_backward'] = default_lanes_each_direction[link['highway']]
        else:
            link['lanes_backward'] = len(link['turn_lanes_backward'].split('|'))

    # correct discrepancies between turns and lanes
    if len(link['turn_lanes'].split('|'))!=link['lanes']:
        adhoc_correct_discrepancy(link)

    if link['lanes_backward']!=0 and len(link['turn_lanes_backward'].split('|'))!=link['lanes_backward']:
        print('ERROR: backward id=',link['id'],' turn=',link['turn_lanes_backward'],' lanes=',link['lanes_backward'])
        adhoc_correct_discrepancy(link)

    return link

def read_node(element):
    node = {}
    node['id'] = element['id']
    node['y'] = element['lat']
    node['x'] = element['lon']
    return node

def parse_jsons(jsons):
    # make sure we got data back from the server requests
    elements = []
    for json in jsons:
        elements.extend(json['elements'])

    # extract nodes and paths from the downloaded osm data
    links = {}
    nodes = {}
    for element in elements:

        if element['type'] == 'way':
            link = read_way(element)
            links[link['id']] = link

        elif element['type'] == 'node':
            node = read_node(element)
            nodes[node['id']] = node

        else:
            print('unknown')

    return links, nodes

def split_streets(links):

    # collect internal and external nodes
    internal_nodes=set()
    external_nodes=set()
    for link_id,link in links.items():
        my_nodes=link['nodes']
        external_nodes.update([my_nodes[-1],my_nodes[0]])
        internal_nodes.update(my_nodes[1:-1])

    # iterate through nodes that are both internal and external
    for node_id in internal_nodes.intersection(external_nodes):

        # links that have this as an internal node
        split_links = [link for link in links.values() if node_id in link['nodes'][1:-1]]

        for link in split_links:
            loc = link['nodes'].index(node_id)

            # create new link representing the first part of the link
            new_link = link.copy()
            new_link['id'] = new_link_id()
            new_link['nodes']= link['nodes'][:loc+1].copy()
            new_link['end_node_id']= node_id
            new_link['turn_lanes'] = '|'*(link['lanes']-1)
            links[new_link['id']] = new_link

            # fix link
            link['nodes'] = link['nodes'][loc:]
            link['start_node_id']= node_id
            link['turn_lanes_backward'] = '|'*(link['lanes_backward']-1)

def flip_wrong_way_links(links):

    flip_links = [link for link in links.values() if link['flip']]

    if len(flip_links)>0:
        print("ERROR: FLIP LINKS ARE NOT IMPLEMENTED!!!!")

    # remove flip_links attribute
    for link in links.values():
        del link['flip']

def expand_bidirectional_links(links):

    bi_links = [link for link in links.values() if link['bidirectional']]

    for link in bi_links:

        start_node_id = link['start_node_id']
        end_node_id = link['end_node_id']

        existing_backward_links = [link for link in links.values() if link['start_node_id']==end_node_id and link['end_node_id']==start_node_id]
        if len(existing_backward_links)>0:
            print("ERROR: FOUND A BACKWARD LINK")

        backward_link = link.copy()

        backward_link['id'] = new_link_id()
        backward_link['start_node_id'] = end_node_id
        backward_link['end_node_id'] = start_node_id
        backward_link['lanes'] = link['lanes_backward']
        backward_link['nodes'] = list(reversed(link['nodes']))
        backward_link['turn_lanes'] = link['turn_lanes_backward']

        # add to links
        links[backward_link['id']] = backward_link

    # remove backward attributes
    for link in links.values():
        del link['lanes_backward']
        del link['turn_lanes_backward']

    print(len(links))

def find_direction(start_node,end_node):
    return math.atan2(end_node['y'] - start_node['y'], end_node['x'] - start_node['x']) * 180 / math.pi

def create_road_connections(links):

    all_road_conns = []

    for link in links.values():

        road_conns = []

        # get next links
        end_node = nodes[link['end_node_id']]
        next_links = [link for link in links.values() if link['id'] in end_node['out_links']]

        # ignore U turns
        next_links = [x for x in next_links if x['end_node_id']!=link['start_node_id']]

        # not needed when there are no turning options
        if len(next_links)<2:
            continue

        # trivial case
        if link['lanes']==1:
            for next_link in next_links:
                road_conns.append({ 'in_link': link['id'] ,
                                    'in_lanes' : [] ,
                                    'out_link' : next_link['id'] ,
                                    'out_lanes' : [] })
            continue

        # compute angles
        d_in = find_direction(nodes[link['nodes'][-2]],nodes[link['nodes'][-1]])

        relative_angle = []
        for next_link in next_links:
            from_node = nodes[next_link['nodes'][0]]
            to_node = nodes[next_link['nodes'][1]]
            d_out = find_direction(from_node,to_node)
            relative_angle.append( ( next_link['id'] , d_out-d_in ) )
        relative_angle = sorted(relative_angle,key=lambda x:x[1])

        # extract all of the turns we have
        turn_lanes = link['turn_lanes'].split('|')
        turn_to_lanes = {}
        for i in range(len(turn_lanes)):
            for t in turn_lanes[i].split(';'):
                if t=='':
                    continue
                if t in turn_to_lanes.keys():
                    turn_to_lanes[t].append(i+1)
                else:
                    turn_to_lanes[t] = [i+1]

        if 'left' in turn_to_lanes:
            road_conns.append({'in_link':link['id'],
                              'in_lanes': turn_to_lanes['left'],
                              'out_link': relative_angle[0][0],
                              'out_lanes':[]})
            next_links.remove(links[relative_angle[0][0]])
            del relative_angle[0]

        if 'right' in turn_to_lanes:
            road_conns.append({'in_link':link['id'],
                               'in_lanes':turn_to_lanes['right'],
                               'out_link':relative_angle[len(relative_angle)-1][0],
                               'out_lanes':[]})
            next_links.remove(links[relative_angle[len(relative_angle)-1][0]])
            del relative_angle[len(relative_angle)-1]


        if ('through' in turn_to_lanes) and len(relative_angle)==1:
            road_conns.append({'in_link':link['id'],
                               'in_lanes':turn_to_lanes['through'],
                               'out_link':relative_angle[0][0],
                               'out_lanes':[]})
            next_links.remove(links[relative_angle[0][0]])
            del relative_angle[0]

        if len(next_links)==0:
            continue

        unmarked_lanes = [i+1 for i in range(len(turn_lanes)) if turn_lanes[i]=='']

        for next_link in next_links:
            road_conns.append({'in_link':link['id'],
                               'in_lanes':unmarked_lanes,
                               'out_link':next_link['id'],
                               'out_lanes':[]})

        # check that all lanes are covered
        lanes_covered = set()
        for road_conn in road_conns:
            if len(road_conn['in_lanes'])==0:
                lanes_covered.update( range(1,link['lanes']+1) )
            else:
                lanes_covered.update(set(road_conn['in_lanes']))

        if lanes_covered != set(range(1,link['lanes']+1)):
            print('ERROR Not all lanes are covered with road connections')

        # check that all next_links are represented
        to_links = set([road_conn['out_link'] for road_conn in road_conns])
        next_links=set([x['id'] for x in links.values() if x['id'] in end_node['out_links'] and x['end_node_id']!=link['start_node_id']])
        if to_links!=next_links:
            print('ERROR Not all outgoing links are covered with road connections')

        all_road_conns.extend(road_conns)

    return all_road_conns

def convert_latlon_to_meters(nodes):
    centroid = np.mean([[v['x'], v['y']] for v in nodes.values()], axis=0)
    for node_id, node in nodes.items():
        dx, dy = latlong2meters(node['y'], node['x'], centroid[1], centroid[0])
        node['x'] = dx
        node['y'] = dy

def latlong2meters(lat, lon, clat, clon):
    R = 6378137    # Radius of earth in meters

    lat = lat*math.pi / 180
    lon = lon*math.pi / 180

    clat = clat*math.pi / 180
    clon = clon*math.pi / 180

    dx = math.acos( 1 - math.pow(math.cos(clat), 2) * (1-math.cos(lon-clon)) ) * R
    if lon<clon:
        dx = -dx
    dy = (lat-clat) * R

    return dx, dy

def compute_link_lengths(links, nodes):
    for link in links.values():
        start_node = nodes[link['start_node_id']]
        end_node = nodes[link['end_node_id']]
        link['length'] = math.sqrt( math.pow(end_node['x']-start_node['x'],2) + math.pow(end_node['y']-start_node['y'],2) )

# def remove_simple_nodes(links, nodes):
#
#     simple_nodes = [node['id'] for node in nodes.values() if len(node['in_links'])==1 and len(node['out_links'])==1]
#
#
#     # discard superfluous nodes
#     for node_id in simple_nodes:
#         if node_id in nodes:
#             del nodes[node_id]
#         for link in [x for x in links.values() if node_id in x['nodes']]:
#             link['nodes'].remove(node_id)
#
#     ### REMOVE SIMPLE NODES FROM LINK.NODES

def create_otm_scenario():
    scenario = etree.Element("scenario")
    etree.SubElement(scenario, "commodities")
    models = etree.SubElement(scenario, "models")
    return scenario

def add_veh_types(scenario, veh_types):
    commodities = next(scenario.iter("commodities"))
    for data in veh_types:
        etree.SubElement(commodities, "commodity", data)

def add_network(scenario, links, nodes, external_nodes, road_conns):

    # Base XML data
    network = etree.SubElement(scenario, "network")

    # road params
    road_params = etree.SubElement(network, "road_params")
    for road_type_name, road_type_params in road_param_types.items():
        etree.SubElement(road_params, "road_param", {
            'id': str(road_type_params['id']),
            'capacity': str(road_type_params['capacity']),
            'speed': str(road_type_params['speed']),
            'jam_density': str(road_type_params['jam_density'])
        })

    # get all node positions
    node_id_map = {}
    node_set = etree.SubElement(network, 'nodes')
    node_id = 0
    for node_osmid in external_nodes:
        node = nodes[node_osmid]
        node_id_map[node['id']] = node_id
        etree.SubElement(node_set, 'node', {
            'id': str(node_id),
            'x': '{:.2f}'.format(node['x']),
            'y': '{:.2f}'.format(node['y'])
        })
        node_id += 1

    link_id_map = {}
    link_set = etree.SubElement(network, "links")
    link_id = 0
    for link_osmid, link in links.items():
        link_id_map[str(link['id'])] = link_id

        elink = etree.SubElement(link_set, 'link', {
            'id': str(link_id),
            'length': '{:.2f}'.format(link['length']),                       # WHAT UNITS IS THIS IN?
            'full_lanes':  str(link['lanes']),                                         # HOW TO GET THE LANES?
            'start_node_id':  str(node_id_map[link['start_node_id']]),
            'end_node_id':  str(node_id_map[link['end_node_id']]),
            'roadparam': str(link['roadparam'])
        })
        link_id += 1

        epoints = etree.SubElement(elink, "points")
        for node_id in link['nodes']:
            node = nodes[node_id]
            etree.SubElement(epoints, "point", {
                'x': '{:.2f}'.format(node['x']),
                'y': '{:.2f}'.format(node['y'])
            })


    # road connections
    rc_id = 0
    roadconnections = etree.SubElement(network, "roadconnections")

    for road_conn in road_conns:
        rc = etree.SubElement(roadconnections, "roadconnection", {
            "id": str(rc_id),
            "in_link": str(road_conn['in_link']),
            "out_link": str(road_conn['out_link']),
        })
        if 'in_link_lanes' in road_conn:
            rc['in_link_lanes'] = '{}#{}'.format(min(road_conn['in_link_lanes']),max(road_conn['in_link_lanes']))
        if 'out_link_lanes' in road_conn:
            rc['out_link_lanes'] = '{}#{}'.format(min(road_conn['out_link_lanes']),max(road_conn['out_link_lanes']))
        rc_id += 1

    print('asdf')
    # # SUBNETWORK DATA
    # subnetworks = SubElement(scenario, "subnetworks")
    # # zip all edge ids to their respective subnetwork ids
    # subnet_pairs = [(graph.edges[edge]['link_id'], graph.edges[edge]['subnetwork_id']) for edge in graph.edges]
    # distinct_subnets = set(list(zip(*subnet_pairs))[1])
    #
    # for subnet_id in distinct_subnets:
    #     subnet_links = sorted([i[0] for i in subnet_pairs if i[1] == subnet_id])
    #     subnetwork = SubElement(subnetworks, "subnetwork", {"id": str(subnet_id)})
    #     subnetwork.text = ",".join([str(i) for i in subnet_links])

def to_xml(scenario):
    # print(etree.tostring(scenario, pretty_print=True))
    with open("blabla.xml", "wb") as xml_file:
        xml_file.write(etree.tostring(scenario, pretty_print=True))
    # with open("blabla.xml", "wt") as xml_file:
    #     xml_file.write(_pretty_tostring(scenario))

def plot_graph(scenario, bbox=None, fig_height=6, fig_width=None, margin=0.02,
               axis_off=True, equal_aspect=False, bgcolor='w', show=True,
               save=False, close=True, file_format='png', filename='temp',
               dpi=300, annotate=False, node_color='#66ccff', node_size=15,
               node_alpha=1, node_edgecolor='none', node_zorder=1,
               edge_color='#999999', edge_linewidth=1, edge_alpha=1,
               use_geom=True):
    """
    Plot a networkx spatial graph.

    Parameters
    ----------
    G : networkx multidigraph
    bbox : tuple
        bounding box as north,south,east,west - if None will calculate from
        spatial extents of data. if passing a bbox, you probably also want to
        pass margin=0 to constrain it.
    fig_height : int
        matplotlib figure height in inches
    fig_width : int
        matplotlib figure width in inches
    margin : float
        relative margin around the figure
    axis_off : bool
        if True turn off the matplotlib axis
    equal_aspect : bool
        if True set the axis aspect ratio equal
    bgcolor : string
        the background color of the figure and axis
    show : bool
        if True, show the figure
    save : bool
        if True, save the figure as an image file to disk
    close : bool
        close the figure (only if show equals False) to prevent display
    file_format : string
        the format of the file to save (e.g., 'jpg', 'png', 'svg')
    filename : string
        the name of the file if saving
    dpi : int
        the resolution of the image file if saving
    annotate : bool
        if True, annotate the nodes in the figure
    node_color : string
        the color of the nodes
    node_size : int
        the size of the nodes
    node_alpha : float
        the opacity of the nodes
    node_edgecolor : string
        the color of the node's marker's border
    node_zorder : int
        zorder to plot nodes, edges are always 2, so make node_zorder 1 to plot
        nodes beneath them or 3 to plot nodes atop them
    edge_color : string
        the color of the edges' lines
    edge_linewidth : float
        the width of the edges' lines
    edge_alpha : float
        the opacity of the edges' lines
    use_geom : bool
        if True, use the spatial geometry attribute of the edges to draw
        geographically accurate edges, rather than just lines straight from node
        to node

    Returns
    -------
    fig, ax : tuple
    """

    network = next(scenario.iter("network"))
    nodes = next(network.iter("nodes"))
    links = next(network.iter("links"))

    # get bounding box
    xs = [float(node.attrib['x']) for node in nodes.iter("node")]
    ys = [float(node.attrib['y']) for node in nodes.iter("node")]

    north = max(ys)
    south = min(ys)
    east = max(xs)
    west = min(xs)

    # if caller did not pass in a fig_width, calculate it proportionately from
    # the fig_height and bounding box aspect ratio
    bbox_aspect_ratio = 1 #(north-south)/(east-west)
    if fig_width is None:
        fig_width = fig_height / bbox_aspect_ratio

    # create the figure and axis
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor=bgcolor)
    ax.set_facecolor(bgcolor)

    # draw the edges as lines from node to node
    lines = []

    for link in links.iter("link"):
        points = [[float(point.attrib['x']),float(point.attrib['y'])] for point in next(link.iter("points")).iter("point")]

        i = 0
        num_points = len(points)
        while i < num_points-1:
            lines.append([points[i], points[i+1]])
            i += 1

    # add the lines to the axis as a linecollection
    lc = LineCollection(lines, colors=edge_color, linewidths=edge_linewidth, alpha=edge_alpha, zorder=2)
    ax.add_collection(lc)

    # scatter plot the nodes
    # node_Xs = [float(x) for _, x in G.nodes(data='x')]
    # node_Ys = [float(y) for _, y in G.nodes(data='y')]
    # ax.scatter(node_Xs, node_Ys, s=node_size, c=node_color, alpha=node_alpha, edgecolor=node_edgecolor, zorder=node_zorder)

    # set the extent of the figure
    margin_ns = (north - south) * margin
    margin_ew = (east - west) * margin
    ax.set_ylim((south - margin_ns, north + margin_ns))
    ax.set_xlim((west - margin_ew, east + margin_ew))

    # configure axis appearance
    xaxis = ax.get_xaxis()
    yaxis = ax.get_yaxis()

    xaxis.get_major_formatter().set_useOffset(False)
    yaxis.get_major_formatter().set_useOffset(False)

    # if axis_off, turn off the axis display set the margins to zero and point
    # the ticks in so there's no space around the plot
    # if axis_off:
    #     ax.axis('off')
    #     ax.margins(0)
    #     ax.tick_params(which='both', direction='in')
    #     xaxis.set_visible(False)
    #     yaxis.set_visible(False)
    #     fig.canvas.draw()

    # if equal_aspect:
    #     # make everything square
    #     ax.set_aspect('equal')
    #     fig.canvas.draw()
    # else:
    #     # if the graph is not projected, conform the aspect ratio to not stretch the plot
    #     if G.graph['crs'] == settings.default_crs:
    #         coslat = np.cos((min(node_Ys) + max(node_Ys)) / 2. / 180. * np.pi)
    #         ax.set_aspect(1. / coslat)
    #         fig.canvas.draw()

    # annotate the axis with node IDs if annotate=True
    # if annotate:
    #     for node, data in G.nodes(data=True):
    #         ax.annotate(node, xy=(data['x'], data['y']))

    # save and show the figure as specified
    plt.show()
    # return fig, ax

##############################################################

# 1. query osm
# jsons = query_json()
# with open('rgiom.p','wb') as file:
#     pickle.dump( jsons, file)
# with open('rgiom.p','rb') as file:
#     jsons = pickle.load(file)

# 2. parse osm
links, nodes = parse_jsons(jsons)

max_link_id = max([x for x in links.keys()]) + 1

# split links when they cross another street
split_streets(links)


# node inlinks and outlinks
for node_id,node in nodes.items():
    node['out_links']=[link_id for link_id,link in links.items() if link['start_node_id']==node_id]
    node['in_links']=[link_id for link_id,link in links.items() if link['end_node_id']==node_id]

# check that there are no more shared nodes
internal_nodes=set()
external_nodes=set()
for link_id,link in links.items():
    my_nodes=link['nodes']
    external_nodes.update([my_nodes[-1],my_nodes[0]])
    internal_nodes.update(my_nodes[1:-1])

# CHECK: all internal nodes should be disconnected
for node_id in internal_nodes:
    node = nodes[node_id]
    if len(node['out_links'])!=0 or len(node['in_links'])!=0:
        print("ERROR 2094bj 240")

# join links between simple external nodes
external_to_internal_nodes = set()
for node_id in external_nodes:
    node = nodes[node_id]
    if len(node['out_links'])==1 and len(node['in_links'])==1:
        in_link = links[node['in_links'][0]]
        out_link = links[node['out_links'][0]]
        if in_link['lanes']==out_link['lanes'] and \
            in_link['lanes_backward']==out_link['lanes_backward'] and \
            in_link['flip']==out_link['flip'] and \
            in_link['bidirectional']==out_link['bidirectional'] and \
            in_link['turn_lanes']==out_link['turn_lanes'] and \
            in_link['turn_lanes_backward']==out_link['turn_lanes_backward'] and \
            in_link['highway']==out_link['highway']:

            in_link['nodes'].extend(out_link['nodes'][1:])
            in_link['end_node_id'] = in_link['nodes'][-1]

            del links[node['out_links'][0]]
            external_to_internal_nodes.add(node_id)
            node['out_links'] = []
            node['in_links'] = []

            end_node = nodes[out_link['end_node_id']]
            end_node['in_links'].remove(out_link['id'])
            end_node['in_links'].append(in_link['id'])

external_nodes = external_nodes.difference(external_to_internal_nodes)
internal_nodes.update(external_to_internal_nodes)

# 3. deal with oneway tag
flip_wrong_way_links(links)
expand_bidirectional_links(links)

# 4. convert latlong to meters
convert_latlon_to_meters(nodes)

# 5. compute link lengths
compute_link_lengths(links, nodes)

# create road connections
road_conns = create_road_connections(links)

# 6. remove internal nodes
# remove_simple_nodes(links, nodes)

# 7. convert to otm
scenario = create_otm_scenario()

# vehicle types
add_veh_types(scenario, [
    {'id': '0', 'name': 'type0', 'pathfull': 'false'}
])

# network
add_network(scenario, links, nodes, external_nodes, road_conns)

# # demands
# source_links = get_source_links(graph)
#
# # pathless demands
# dt = 300
# vph = [0, 4, 3, 2, 5]
# convert.add_pathless_demand(scenario, veh_types[0], source_links[0], dt, vph)
#
# # splits
# split_nodes = convert.get_split_nodes(graph)
#
# node_io = convert.get_node_io(split_nodes[0], graph)
# link_out_to_split_profile = dict
# convert.add_split(split_nodes[0], veh_types[0], node_io.link_in[0], dt, link_out_to_split_profile)

# # pathfull demands
# convert.add_pathfull_demand(scenario,)

to_xml(scenario)

plot_graph(scenario)

print('done')