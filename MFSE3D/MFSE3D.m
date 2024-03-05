function [eIntopMat,ptIntopMat] = MFSE3D(nptx,npty,nptz,refine,corlence_x,corlence_y,corlence_z)
ptdist = 1;
%% 各向异性相关长度 %%
corlen_x = max(1,corlence_x * nptx);
corlen_y = max(1,corlence_y * npty);
corlen_z = max(1,corlence_z * nptz);
%% 各向同性相关长度 %%
% corlen_x = min(min(corlen_x,corlen_y),corlen_z);
% corlen_y = corlen_x; corlen_z = corlen_x;
elelen = ptdist/refine;
nelx = refine*nptx; nely = refine*npty; nelz = refine*nptz;
tolne = nelx*nely*nelz;
tolnpt = nptx*npty*nptz;
% [Ypt, Xpt ,Zpt] = meshgrid((npty-0.5:-1:0.5)*ptdist, (0.5:1:nptx)*ptdist, (0.5:1:nptz)*ptdist);
[Xpt, Ypt ,Zpt] = meshgrid((0.5:1:nptx)*ptdist, (npty-0.5:-1:0.5)*ptdist,  (0.5:1:nptz)*ptdist);

Xpt = Xpt(:); Ypt = Ypt(:); Zpt = Zpt(:);
corMat = zeros(tolnpt,tolnpt);
for i = 1:size(corMat,1)
    for j = i + 1:size(corMat,2)
            corMat(i,j)=exp(-((Xpt(j)-Xpt(i))^2/corlen_x^2+(Ypt(j)-Ypt(i))^2/corlen_y^2+(Zpt(j)-Zpt(i))^2/corlen_z^2));
    end
end
corMat = corMat + corMat';
for i = 1:size(corMat, 1)
    corMat(i,i) = 1;
end

%% 材料场级数展开
if size(corMat,1) < 1e4
    [eigfunMat, eigvalMat] = eig(corMat);
else
    [eigfunMat, eigvalMat] = eigs(corMat,1500);
end
eigvalVec = diag(eigvalMat);
[eigvalVec, eigsortind] = sort(eigvalVec, 'descend');
neig = 0; tmpsum = 0.;
while tmpsum<(1-1e-6)*sum(abs(eigvalVec))
    neig = neig + 1;
    tmpsum = tmpsum+eigvalVec(neig);
end
EXPANMat = sparse(1:neig, 1:neig, eigvalVec(1:neig).^(-1/2), neig, neig)...
    *eigfunMat(:,eigsortind(1:neig))';
clear eigfunMat;

%% 计算插值阵
[Xe, Ye ,Ze] = meshgrid((0.5:1:nelx)*elelen, (nely-0.5:-1:0.5)*elelen, (0.5:1:nelz)*elelen);
Xe = Xe(:); Ye = Ye(:); Ze = Ze(:);

% [Ye, Xe, Ze] = meshgrid((nely-0.5:-1:0.5)*elelen, (0.5:1:nelx)*elelen, (0.5:1:nelz)*elelen);
% Xe = Xe(:); Ye = Ye(:); Ze = Ze(:);

%% 组装eIntopMat
eIntopMat = zeros(neig, tolne);
grsize = min(round(tolnpt/20), tolne);
ngr = ceil(tolne/grsize);
for igr = 1:ngr %% 构建插值矩阵（节省内存）
    eind = (igr-1)*grsize+1:min(igr*grsize, tolne);
    Xe_sub = Xe(eind); Ye_sub = Ye(eind); Ze_sub = Ze(eind);
    eptvals = exp(-((repmat(Xpt',length(eind),1)-repmat(Xe_sub, 1, tolnpt)).^2/corlen_x^2 ...
        +(repmat(Ypt',length(eind),1)-repmat(Ye_sub, 1, tolnpt)).^2/corlen_y^2 ...
        +(repmat(Zpt',length(eind),1)-repmat(Ze_sub, 1, tolnpt)).^2/corlen_z^2))';
    eptvals(abs(eptvals) < 1e-9) = 0;
    eIntopMat(:,eind) = EXPANMat*eptvals;
end
ptIntopMat = EXPANMat*corMat';
end