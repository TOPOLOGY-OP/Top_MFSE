function [eIntopMat,ptIntopMat] = MFSE2D(nptx,npty,refine,corlencoe,type)
ptdist = 1;
if type == 1 % ISOTROPIC
    corlen = corlencoe*min(nptx,npty)*ptdist;
else % ANISOTROPIC
    corlen_x = corlencoe*nptx*ptdist; corlen_y = corlencoe*npty*ptdist; 
end
elelen = ptdist/refine; nelx = refine*nptx; nely = refine*npty;
tolne = nelx*nely;  tolnpt = nptx*npty;
%% BUILD CORRELATION MATRIX
[Xpt, Ypt] = meshgrid((0.5:1:nptx)*ptdist, (npty-0.5:-1:0.5)*ptdist);
Xpt = Xpt(:); Ypt = Ypt(:);
corMat = zeros(tolnpt,tolnpt);
for i = 1:size(corMat,1)
    for j = i+1:size(corMat,2)
        if type == 1
            corMat(i,j) = exp(-(((Xpt(j)-Xpt(i))^2+(Ypt(j)-Ypt(i))^2)/corlen^2));
        else
            corMat(i,j) = exp(-((Xpt(j)-Xpt(i))^2/corlen_x^2+(Ypt(j)-Ypt(i))^2/corlen_y^2));
        end
    end
end
corMat = corMat+corMat';
for i = 1:size(corMat, 1)
    corMat(i,i) = 1;
end
%% DO SERIES EXPANSION OF THE MATERIAL FIELD
if size(corMat,1)<1e4
    [eigfunMat, eigvalMat] = eig(corMat);
else
    [eigfunMat, eigvalMat] = eigs(corMat,1500);
end
eigvalVec = diag(eigvalMat);
[eigvalVec, eigsortind] = sort(eigvalVec, 'descend');
neig = 0; tmpsum = 0.;
while tmpsum<(1-1e-4)*sum(abs(eigvalVec))
    neig = neig+1;
    tmpsum = tmpsum+eigvalVec(neig);
end
EXPANMat = sparse(1:neig, 1:neig, eigvalVec(1:neig).^(-1/2), neig, neig)...
    *eigfunMat(:,eigsortind(1:neig))'; clear eigfunMat;
% COMPUTE PHI ON ELEMENTS AND MATERIAL-FIELD POINTS
[Xe, Ye] = meshgrid((0.5:1:nelx)*elelen, (nely-0.5:-1:0.5)*elelen);
Xe = Xe(:); Ye = Ye(:);
eIntopMat = zeros(neig, tolne);
grsize = min(round(tolnpt/1), tolne); ngr = ceil(tolne/grsize);
for igr = 1:ngr
    eind = (igr-1)*grsize+1:min(igr*grsize, tolne);
    Xe_sub = Xe(eind); Ye_sub = Ye(eind);
    if type == 1
        eptvals = exp(-(((repmat(Xpt',length(eind),1)-repmat(Xe_sub, 1, tolnpt)).^2 ...
            +(repmat(Ypt',length(eind),1)-repmat(Ye_sub, 1, tolnpt)).^2)/corlen^2))';
    else
        eptvals = exp(-(((repmat(Xpt',length(eind),1)-repmat(Xe_sub, 1, tolnpt)).^2/corlen_x^2 ...
            +(repmat(Ypt',length(eind),1)-repmat(Ye_sub, 1, tolnpt)).^2/corlen_y^2)))';
    end
    eptvals(abs(eptvals)<1e-9)=0;
    eIntopMat(:,eind) = EXPANMat*eptvals;
end
ptIntopMat = EXPANMat*corMat'; clear corMat;
end