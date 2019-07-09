import osm_loader
from lxml import etree

road_param_types = {
    'residential':      {'id': 0, 'capacity': 2000, 'speed': 100, 'jam_density': 100},
    'tertiary':         {'id': 1, 'capacity': 2001, 'speed': 101, 'jam_density': 101},
    'tertiary_link':    {'id': 2, 'capacity': 2002, 'speed': 102, 'jam_density': 102},
    'unclassified':     {'id': 3, 'capacity': 2003, 'speed': 103, 'jam_density': 103},
    'secondary':        {'id': 4, 'capacity': 2004, 'speed': 104, 'jam_density': 104},
}

class OTM:

    def __init__(self, cfg=''):
        self.config_file = cfg
        self.scenario = {}

    def load_from_osm(self,west=-122.2981,north=37.8790,east=-122.2547,south=37.8594,fixes={}):
        self.scenario = osm_loader.load_from_osm(west=west,north=north,east=east,south=south,fixes=fixes)

    def save_to_xml(self, outfile):

        scenario=etree.Element('scenario')
        etree.SubElement(scenario,'commodities')
        models = etree.SubElement(scenario, 'models')
        model = etree.SubElement(models,'model',{'type':'ctm','name':'myctm','is_default':'true'})
        etree.SubElement(model,'model_params',{'sim_dt':'2','max_cell_length':'100'})

        # vehicle types
        commodities=next(scenario.iter('commodities'))
        etree.SubElement(commodities,'commodity',{'id':'0','name':'type0','pathfull':'false'})

        # network
        network=etree.SubElement(scenario,'network')

        # road params
        road_params=etree.SubElement(network,'roadparams')
        for road_type_name,road_type_params in road_param_types.items():
            etree.SubElement(road_params,'roadparam',{
                'id':str(road_type_params['id']),
                'capacity':str(road_type_params['capacity']),
                'speed':str(road_type_params['speed']),
                'jam_density':str(road_type_params['jam_density'])
            })

        # get all node positions
        # node_id_map={}
        node_set=etree.SubElement(network,'nodes')
        node_id=0
        for node_osmid in self.scenario['external_nodes']:
            node = self.scenario['nodes'][node_osmid]
            # node_id_map[node['id']]=node_osmid
            etree.SubElement(node_set,'node',{
                'id':str(node_osmid),
                'x':'{:.2f}'.format(node['x']),
                'y':'{:.2f}'.format(node['y'])
            })
            node_id+=1

        link_id_map={}
        link_set=etree.SubElement(network,'links')
        link_id=0
        for link_osmid,link in self.scenario['links'].items():
            link_id_map[str(link['id'])]=link_id

            if link['start_node_id'] not in self.scenario['nodes'].keys():
                print('ERROR: link[start_node_id] not in nodes.keys()')

            if link['end_node_id'] not in self.scenario['nodes'].keys():
                print('ERROR: link[end_node_id] not in nodes.keys()')

            elink=etree.SubElement(link_set,'link',{
                'id':str(link_osmid),
                'length':'{:.2f}'.format(link['length']),  # WHAT UNITS IS THIS IN?
                'full_lanes':str(link['lanes']),  # HOW TO GET THE LANES?
                'start_node_id':str(link['start_node_id']),
                'end_node_id':str(link['end_node_id']),
                'roadparam':str(link['roadparam'])
            })
            link_id+=1

            epoints=etree.SubElement(elink,'points')
            for node_id in link['nodes']:
                node=self.scenario['nodes'][node_id]
                etree.SubElement(epoints,'point',{
                    'x':'{:.2f}'.format(node['x']),
                    'y':'{:.2f}'.format(node['y'])
                })

        # road connections
        rc_id=0
        roadconnections=etree.SubElement(network,'roadconnections')

        for road_conn in self.scenario['road_conns']:
            rc=etree.SubElement(roadconnections,'roadconnection',{
                'id':str(rc_id),
                'in_link':str(road_conn['in_link']),
                'out_link':str(road_conn['out_link']),
            })
            if 'in_link_lanes' in road_conn:
                rc['in_link_lanes']='{}#{}'.format(min(road_conn['in_link_lanes']),max(road_conn['in_link_lanes']))
            if 'out_link_lanes' in road_conn:
                rc['out_link_lanes']='{}#{}'.format(min(road_conn['out_link_lanes']),max(road_conn['out_link_lanes']))
            rc_id+=1

        # actuators
        actuators=etree.SubElement(scenario,'actuators')

        nodes_with_signals = [node for node in self.scenario['nodes'].values() if node['type']=='traffic_signals']
        actuator_id = 0
        for node in nodes_with_signals:

            in_links = [self.scenario['links'][link_id] for link_id in node['in_links']]

            road_conns = [road_conn for road_conn in self.scenario['road_conns'] if road_conn['in_link'] in node['in_links']]


            ### FIGURE OUT PHASES / ROAD CONNECTIONS

            actuator=etree.SubElement(actuators,'actuator',{'id':str(actuator_id),'type':'signal'})
            etree.SubElement(actuator,'actuator_target',{'type':'node','id':str(node['id'])})
            signal = etree.SubElement(actuator,'signal')
            actuator_id += 1

        # # SUBNETWORK DATA
        # subnetworks = SubElement(scenario, 'subnetworks')
        # # zip all edge ids to their respective subnetwork ids
        # subnet_pairs = [(graph.edges[edge]['link_id'], graph.edges[edge]['subnetwork_id']) for edge in graph.edges]
        # distinct_subnets = set(list(zip(*subnet_pairs))[1])
        #
        # for subnet_id in distinct_subnets:
        #     subnet_links = sorted([i[0] for i in subnet_pairs if i[1] == subnet_id])
        #     subnetwork = SubElement(subnetworks, 'subnetwork', {'id': str(subnet_id)})
        #     subnetwork.text = ','.join([str(i) for i in subnet_links])

        with open(outfile,'wb') as xml_file:
            xml_file.write(etree.tostring(scenario,pretty_print=True))

