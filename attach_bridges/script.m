
clear
close all

% load json
fname = 'EQ mag7.5_data.json';
fid = fopen(fname);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
X = jsondecode(str);
clear fname fid raw str

% extract information
num_bridges = numel(X.Bridge);
bridges = repmat(struct('id',nan,'lat',nan,'lng',nan),1,num_bridges);
for i=1:num_bridges
    bridges(i).id = str2double(X.Bridge(i).Site.StructureNumber);
    bridges(i).lat = X.Bridge(i).Site.Location.Latitude;
    bridges(i).lng = X.Bridge(i).Site.Location.Longitude;
end
clear X sites

% load xml
otm = OTMWrapper('anchorage.xml',false);

otm.show_network(5);


figure(otm.fig)
plot([bridge.lat],[bridges.lng],'o')


% % snapping
% for i=1:num_bridges
% 
%     lat = bridges(i).lat;
%     lng = bridges(i).lng;
%     
%     
%     for j=1:otm.api.scenario().get_num_links()
%     
% end


