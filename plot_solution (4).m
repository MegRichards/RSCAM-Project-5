function plot_solution(polygon, type1, type2, varargin)
xv = polygon.Vertices(:,1);
yv = polygon.Vertices(:,2);
[X,Y] = meshgrid(linspace(-1,1,100),linspace(-1,1,100));
if or(type2 == "exterior", type2 == "transmission")
    [X,Y] = meshgrid(linspace(-1.5,1.5,100),linspace(-1.5,1.5,100));
end
[in,~] = inpolygon(X,Y,xv,yv);
In = nan(size(in)); In(in)= 1;
Out = nan(size(in)); Out(boolean(-(in-1)))= 1;
Z = X+1i*Y;
%plot polygon
corners = polygon.Vertices;
corners(end+1,:) = corners(1,:);

if type1 == "Laplace"
    if type2 == "interior"
        A = varargin{1};
        weights = varargin{2};
        u = zeros(size(Z));
        for i=1:size(Z,2)
            u(:,i) = A(Z(:,i))*weights;
        end
        figure, plot(corners(:,1),corners(:,2),'k', 'LineWidth',1); hold on,
        contour(X,Y,u.*In,30);
        %figure, fig = imagesc(u.*In);
        %set(fig,'AlphaData',~isnan(In))
        
    elseif type2 == "exterior"
        u = zeros(size(Z));
        A = varargin{1};
        weights = varargin{2};
        H = varargin{3};
        for i=1:size(Z,2)
            u(:,i) = A(Z(:,i))*weights + H(Z(:,i));
        end
        figure, plot(corners(:,1),corners(:,2),'k', 'LineWidth',1); hold on,
        contour(X,Y,u.*Out,30);
        %figure, fig = imagesc(u.*Out);
        %set(fig,'AlphaData',~isnan(Out))
        
    elseif type2 == "transmission"
        u_in = zeros(size(Z));
        u_ex = zeros(size(Z));
        A = varargin{1};
        weights = varargin{2};
        H = varargin{3};
        A_in = A{1}; A_ex = A{2};
        weights_in = weights{1}; weights_ex = weights{2};
        for i=1:size(Z,2)
            u_in(:,i) = A_in(Z(:,i))*weights_in;
            u_ex(:,i) = A_ex(Z(:,i))*weights_ex + H(Z(:,i));
        end
        figure, plot(corners(:,1),corners(:,2),'k', 'LineWidth',1); hold on,
        contour(X,Y,u_in.*In,30);
        contour(X,Y,u_ex.*Out,30);

        %figure, fig = imagesc(X(1,:),Y(:,1),flipud(real(u_ex).*Out));
        %imagesc(X(1,:),Y(:,1),flipud(real(u_in).*In));
        %set(fig,'AlphaData',~isnan(flipud(Out))), hold on;
        %plot(corners(:,1),-corners(:,2),'k', 'LineWidth',1); hold on,
    end
    
elseif type1 == "Helmholtz"
    if type2 == "exterior"
        u = zeros(size(Z));
        A = varargin{1};
        weights = varargin{2};
        H = varargin{3};
        for i=1:size(Z,2)
            u(:,i) = A(Z(:,i))*weights + H(Z(:,i));
        end
        %contour(X,Y,real(u).*Out,30);
        %u(u>2) = 1;
        %u(u<-2) = -1;
        figure, fig = imagesc(X(1,:),Y(:,1),flipud(real(u).*Out));
        set(fig,'AlphaData',~isnan(flipud(Out))), hold on;
        plot(corners(:,1),-corners(:,2),'k', 'LineWidth',1); hold on,
    end
end