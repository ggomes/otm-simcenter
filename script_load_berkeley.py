from otm import OTM

otm = OTM()
otm.load_from_osm(
    west=-122.2981,north=37.8790,east=-122.2547,south=37.8594,
    simplify_roundabouts=True,
    fixes={
        # turns: left | | | lanes 3
        # OSM: https://www.openstreetmap.org/way/415803770#map=19/37.87693/-122.28277
        # Google: https://www.google.com/maps/@37.8768752,-122.2828014,76m/data=!3m1!1e3
        # Correction: Turn includes the parking lane. Ignore it
        415803770:[('turn_lanes','left||')],

        # turns: left | through | through lanes 2
        # OSM: https://www.openstreetmap.org/way/415876791#map=19/37.87453/-122.26411
        # Google: https://www.google.com/maps/@37.8745003,-122.2648584,78a,35y,94.72h,32.94t/data=!3m1!1e3
        # Correction: Missing lane
        415876791:[('lanes',3),('lanes_backward',2),('turn_lanes_backward','|')],

        # Hearst @ LeConte
        # turns: left | through | through lanes 2
        # OSM: https://www.openstreetmap.org/way/574381942#map=19/37.87458/-122.26374
        # Google: https://www.google.com/maps/@37.8745017,-122.2643212,88m/data=!3m1!1e3
        # Correction: lanes -> 3
        574381942:[('lanes',3),('lanes_backward',2),('turn_lanes_backward','|')]
    }
)


otm.save_to_xml('filename.xml')

print('DONE')