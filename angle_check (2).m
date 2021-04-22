function [angle] = angle_check(polygon, ncorner)
% Checks angle of an interior corner of a polygon
% if output<0, interior angle is greater than pi
corners = polygon.Vertices(:,1) + 1i*polygon.Vertices(:,2); % corners (or vertices) of the polygon
a1 = circshift(corners,-1) - corners; a1 = a1./abs(a1);
a2 = circshift(corners,1) - corners; a2 = a2./abs(a2);

% get unit vectors
v = a1(ncorner);
u = a2(ncorner);

%angle 
angle_1 = atan2(imag(v),real(v))-atan2(imag(u),real(u));
angle = mod(angle_1,2*pi);
end