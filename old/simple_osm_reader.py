
# coding: utf-8

# In[206]:

 
import json
from io import BytesIO
import xml.dom.minidom
from xml.etree.ElementTree import Element, SubElement, tostring


# In[207]:


def pretty_tostring(xml_data):
    xml_data = xml.dom.minidom.parseString(tostring(xml_data))
    xml_pretty_str = xml_data.toprettyxml(indent="  ")
    return xml_pretty_str


# In[208]:


def read_way(element):
    if 'tags' in element:
        tags = element['tags']

        highway = tags['highway']

        if 'oneway' in tags:
            lanes = tags['lanes']
        else:
            lanes = 1

        if 'maxspeed' in tags:
            maxspeed = tags['maxspeed']
        else:
            maxspeed = 20

        if 'name' in tags:
            name = tags['name']
        else:
            name = 'unnamed'

        if 'oneway' in tags:
            oneway = tags['oneway']
        else:
            oneway = 'no'

    else:
        highway = ''
        lanes = 1
        maxspeed = 20
        name = 'unnamed'
        oneway = 'no'
        
        
    if 'nodes' in element:
        nodes = element['nodes']
    else:
        print("NO NODES!!!!")
    
    link = {}
    link['id'] = element['id']
    link['length'] = 0                       # assign later
    link['start_node_id'] = nodes[0]
    link['end_node_id'] = nodes[-1]
    link['roadparam'] = 0
    link['full_lanes'] = lanes
    link['nodes'] = nodes
    return link


# In[209]:


def read_node(element):
    node = {}
    node['id'] = element['id']
    node['x'] = element['lat']
    node['y'] = element['lon']
    return node


# In[210]:


with open("openstreetbrowser.osm.json", "r") as read_file:
    data = json.load(read_file)


# In[211]:


links = {}
nodes = {}
for element in data['elements']:
    
    id = element['id']
    
    # NOTE: Keep only primary, secondary, tertiary (maybe)
    if element['type']=='way':
        link = read_way(element)
        links[link['id']] = link
         
    elif element['type']=='node':
        node = read_node(element)  
        nodes[node['id']] = node    
        
    else:
        print('unknown')


# In[212]:


# link points
discard_nodes = []
for link_id, link in links.items():
    
    my_nodes = link['nodes']
    
    # discard internal nodes
    if len(my_nodes)>2:
        discard_nodes.extend(my_nodes[1:-2])
    
    # add node positions to link
    points = []
    link['points'] = points
    for my_node_id in my_nodes:
        node = nodes[my_node_id]
        point = {}
        point['x'] = node['x']
        point['y'] = node['y']
        points.append(point)
            
# discard superfluous nodes
for key in discard_nodes:
    if key in nodes:
        del nodes[key]
        


# In[213]:


scenario = Element('scenario')
network = SubElement(scenario,'network')

xlinks = SubElement(network,'links')
for link_id,link in links.items():
    xlink = SubElement(xlinks,'link',{
        'id':str(link['id']),
        'full_lanes':str(link['full_lanes']),
        'length':str(link['length']),
        'start_node_id':str(link['start_node_id']),
        'end_node_id':str(link['end_node_id']),
        'roadparam':'0'
        })
    
    xpoints = SubElement(xlink,'points')
    for point in link['points']:
        xpoint = SubElement(xpoints,'point',{
            'x':str(point['x']),
            'y':str(point['y'])
        })
        
    
xnodes = SubElement(network,'links')
for node_id,node in nodes.items():
    xnode = SubElement(xnodes,'node',{
        'id':str(node['id']),
        'x':str(node['x']),
        'y':str(node['y'])
    })
        
# ET.dump(scenario)

with open("bla.xml", "wt") as xml_file:
    xml_file.write(pretty_tostring(scenario))

