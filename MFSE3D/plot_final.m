function plot_final(plot_materials)
[nely,nelx,nelz] = size(plot_materials);
tolne = nelx*nely*nelz;
tolnn = (nelx+1)*(nely+1)*(nelz+1);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:) + 1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],tolne,1);
index = zeros(tolne,8);
Z3 = plot_materials; num_nodes = 0; center_crop = 0.5;
for i = 1:tolne
    for j = 1:8
        index(i,j)=((edofMat(i,3*(j-1)+1))-1)/3+1;
    end
    if Z3(i) > center_crop
        num_nodes = num_nodes + 1;
    end
end
indexm = zeros(num_nodes,8); num_elements = 0;
for i = 1:tolne  
    if Z3(i) > center_crop
        num_elements = num_elements + 1;
        indexm(num_elements,:) = index(i,:);
    end
end
mindex = indexm(:); tolnnn = (1:tolnn); M1 = zeros(tolnn,1);
Mindex = intersect(tolnnn,mindex);
MMindex = setdiff(tolnnn,Mindex)'; M1(Mindex,:) = 1; M1(MMindex,:) = 0;
M1 = reshape(M1,nely+1,nelx+1,nelz+1);
plot_materials = M1;
plot_materials = smooth3(plot_materials,'box',3);
p1 = patch(isosurface(plot_materials,0.5),'FaceColor','b','EdgeColor','none');
p2 = patch(isocaps(plot_materials,0.5),'FaceColor','b','EdgeColor','none');
isonormals(plot_materials,p1);
isonormals(plot_materials,p2);
p = [p1,p2]; direction = [1 0 0];
rotate(p,direction,-90);
lighting flat;
view([23,34]);
axis vis3d tight;axis equal tight off;
camlight('right');pause(1e-6);
end