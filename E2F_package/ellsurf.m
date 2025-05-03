function ellsurf(ell,la0,lo0)
    equator_distance_degrees = (2 * pi * 6371.0) / 360.0; 
    evla_delta_latitude =  ell(:,2) ./ equator_distance_degrees; 
    ell(:,2) = la0 + evla_delta_latitude; 
    evlo_delta_longitude = ell(:,1) ./ (equator_distance_degrees * cosd(la0)); 
    ell(:,1) = lo0 + evlo_delta_longitude; 
    set(gcf, 'Renderer', 'opengl');
    for i=1:length(unique(ell(:,end)))
         num=find(ell(:,end)==i);
         f=convhull(ell(num,1:3));
         patch('vertices',abs([ell(num,1:3)]), ...
               'faces',f,'facecolor','red','facealpha',0.1,'linestyle','none');hold on
         
    end
end

