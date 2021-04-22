function [b1, b2, n1, n2] = new_boundary_normal(polygon, ncorner, dist)

corners = polygon.Vertices(:,1) + 1i*polygon.Vertices(:,2); % corners (or vertices) of the polygon
a1 = circshift(corners,-1) - corners; a1 = a1./abs(a1);
a2 = circshift(corners,1) - corners; a2 = a2./abs(a2);
b1 = corners(ncorner) + dist*a1(ncorner);
b2 = corners(ncorner) + dist*a2(ncorner);
n1 = a1(ncorner)*1i;
n2 = -a2(ncorner)*1i;
