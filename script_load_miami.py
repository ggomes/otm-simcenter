from otm import OTM

otm = OTM()
otm.load_from_osm(
    north=25.9841,south=25.1923,east=-80.1107,west=-80.8593,
    simplify_roundabouts=False,
    fixes={
    }
)

X = otm.get_link_table()
X.sort_values(by='travel_time', ascending=True, inplace=True)
print(X)

otm.merge_nodes([])

print(X)
otm.save_to_xml('miami.xml')

print('DONE')




# ERROR: id= 651292620  turn= left  lanes= 2
# ERROR: id= 651633290  turn= left  lanes= 2
# ERROR: id= 651633296  turn= left  lanes= 2
# ERROR: backward id= 651633296  turn= left  lanes= 2
# ERROR: id= 651633298  turn= left  lanes= 2
# ERROR: backward id= 651633298  turn= left  lanes= 2
# ERROR: id= 651634719  turn= left  lanes= 2
# ERROR: backward id= 651634719  turn= left  lanes= 2
# ERROR: id= 651634720  turn= left  lanes= 2
# ERROR: backward id= 651634720  turn= left  lanes= 2


# ERROR 'turn:lanes:both_ways' in tags and 'turn:lanes:forward' in tags
# OSM:
# Google:
# Correction:

