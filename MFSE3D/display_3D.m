function display_3D(rho)
mymap = zeros(500,3);
for i = 1:100
    mymap(i,1) = 0;
    mymap(i,2) = i/100;
    mymap(i,3) = 1;
end
for i = 101:200
    mymap(i,1) = 0;
    mymap(i,2) = 1;
    mymap(i,3) = 1-(i-100)/100;
end
for i = 201:300
    mymap(i,1) = (i-200)/100;
    mymap(i,2) = 1;
    mymap(i,3) = 0;
end
for i = 301:400
    mymap(i,1)=1;
    mymap(i,2)=1-(i-300)/200;
    mymap(i,3)=0;
end
for i = 401:500
    mymap(i,1)=1;
    mymap(i,2)=0.5-(i-400)/200;
    mymap(i,3)=0;
end
[nely,nelx,nelz] = size(rho);
rho_temp = repmat(rho(:)',6,1);
rho_temp = rho_temp(:);
freedofs = [];
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
F = zeros(nelx*nely*nelz*6,4); V = zeros(nelx*nely*nelz*8,3);
temp = 0;
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.5)
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                V(8*temp+1:8*(temp+1),:) = vert;
                F(6*temp+1:6*(temp+1),:) = face;               
            else
                freedofs = [freedofs,6*temp+1:6*(temp+1)];
            end
            face = face + 8;
            temp = temp + 1;
        end
    end
end
F(freedofs,:) = []; rho_temp(freedofs,:) = [];
patch('Faces',F,'Vertices',V,'FaceVertexCData',rho_temp,'FaceColor','flat','EdgeColor','none');
hold on;
colormap(mymap); hco = colorbar; caxis([0 1]);
set(hco,'Ticks',0:0.1:1,'ylim',[0 1]);
set(hco,'YTickLabel',num2str(get(hco,'YTick')','%.1f'))
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end