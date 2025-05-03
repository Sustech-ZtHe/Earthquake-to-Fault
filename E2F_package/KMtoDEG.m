function points=KMtoDEG(points,la0,lo0)
if ~isempty(points)
    equator_distance_degrees = (2 * pi * 6371.0) / 360.0; 
    evla_delta_latitude =  points(:,2) ./ equator_distance_degrees; 
    points(:,2) = la0 + evla_delta_latitude; 
    evlo_delta_longitude = points(:,1) ./ (equator_distance_degrees * cosd(la0)); 
    points(:,1) = lo0 + evlo_delta_longitude; 
end