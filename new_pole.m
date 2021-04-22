function pole = new_pole(polygon, ncorner, dist, theta, interior)
% polygon = a polygon object, as created with the Matlab function polyshape
% ncorner = the number (or index) of the corner where the new pole should be created
% dist = distance of the new pole from the corner
% theta = angle between the new pole and the bissectrice of the corner
% interior = boolean variable indicating "true" for a new pole inside the
%            polygon and false outside the polygon

corners = polygon.Vertices(:,1) + 1i*polygon.Vertices(:,2); % corners (or vertices) of the polygon
side1 = corners - circshift(corners,1); % arrow from previous corner to corner
side2 = circshift(corners,-1) - corners; % arrow from corner to next corner
phi = (wrapTo2Pi(angle(-side1./side2))-2*pi)/2;% angle bisecting the corner (inward)
dir = side2./abs(side2).*exp(1i*phi); % unitary arrow bisecting the corner (inward)

if interior
    pole = corners(ncorner) + dist*dir(ncorner)*exp(1i*theta); 
else
    pole = corners(ncorner) - dist*dir(ncorner)*exp(1i*theta);
end
