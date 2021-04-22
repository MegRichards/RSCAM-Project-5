function [P, normal] = discretization_boundary(polygon, N)

corners = polygon.Vertices(:,1) + 1i*polygon.Vertices(:,2); % corners (or vertices) of the polygon
alpha = linspace(1/N,1-1/N,N);
P = corners*(1-alpha) + circshift(corners,-1)*alpha;
normal = reshape(repmat(((P(:,2) - P(:,1))*1i),1,N),size(P,1)*size(P,2),1);
normal = (normal./abs(normal));
P = reshape(P,size(P,1)*size(P,2),1);
%% Plot
%scaling_normal = 0.1;
%figure,
%for i = 1:length(P)
%plot(real(P(i)), imag(P(i)), '.', 'color', 'k'), hold on
%quiver(real(P(i)),imag(P(i)),real(scaling_normal*normal(i)),scaling_normal*imag(normal(i)))
end
