from otm import OTM
import pandas as pd

# north=61.2597,south=61.0672,east=-149.6302,west=-150.0446,

otm = OTM()
otm.load_from_osm(
    north=61.2597,south=61.0672,east=-149.6302,west=-150.0446,
    simplify_roundabouts=False,
    fixes={
        # ERROR 'turn:lanes:both_ways' in tags and 'turn:lanes:forward' in tags
        # OSM: https://www.openstreetmap.org/way/651292620#map=18/61.13014/-149.88261
        # Google: https://www.google.com/maps/@61.1302853,-149.8829239,341m/data=!3m1!1e3
        651292620 : [('lanes:both_ways',None),('turn:lanes:both_ways',None)]
    }
)

X = otm.get_link_table()
X.sort_values(by='travel_time', ascending=True, inplace=True)

print(X)
otm.save_to_xml('anchorage.xml')

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

