function Top_mfse3D(nptx,npty,nptz,refine,volfra,corlencoe)
%% DEFINE PARAMETERS
maxloop = 1000;
tolx = 0.01;
displayflag = 2;
solver = 1;
%% SETUP PARAMETERS
E0 = 200e5; Emin = E0*1e-9; nu = 0.3; ptdist = 1;
elelen = ptdist/refine;
nelx = refine*nptx; nely = refine*npty; nelz = refine*nptz;
tolne = nelx*nely*nelz; tolnpt = nptx*npty*nptz; tolvol = tolne*elelen^3;
fprintf([' Number of material-field  points:%10i \n'...
    ' Number of finite elements:%10i\n'],tolnpt,tolne);
%% SETUP BOUNDARY CONDITIONS
il = nelx; jl = 0; kl = 0:nelz;                         % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
[jf,kf] = meshgrid(1:nely+1,1:nelz+1);                  % Coordinates
fixednid = (kf-1)*(nely+1)*(nelx+1)+jf;                 % Node IDs
temp = 1:(nelx+1)*(nely+1);
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2; 3*temp(:)];
%% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-100,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
KE = lk_H8(nu);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
%% MATERIAL FIELD SERIES EXPANSION
[eIntopMat,ptIntopMat] = MFSE3D(nptx,npty,nptz,refine,corlencoe,corlencoe,corlencoe);
save('MFSE3D.mat','eIntopMat','ptIntopMat');
load('MFSE3D.mat');
%% INITIALIZE DESIGN VARIABLES
beta = 0.5; penal = 3;
x = (-log(1/volfra-1)/beta)*ones(1,tolnpt)/ptIntopMat; x = x';
%% INITIALIZE ITERATION
loop = 0; obj = 0.; neig = length(x);
change = 1.; ichange = 1; n = neig;
xmin = -1000*ones(n,1); xmax = 1000*ones(n,1);
low = xmin; upp = xmax;
xold1 = x;  xold2 = x; clf;
Obj = []; Volf = [];
%% START ITERATION
while (change > tolx || beta<20) && loop < maxloop
    loop = loop + 1; objold = obj;
    %% MATERIAL FIELD CONVERSION
    ePhi = eIntopMat'*x;
    [ePhiProj,edproj] = THRESHOLD(ePhi,beta);
    %% DO FINITE ELEMENT ANALYSIS
    sK = reshape(KE(:)*(Emin+ePhiProj'.^penal*(E0-Emin)), 24*24*tolne, 1);
    K = sparse(iK,jK,sK); K = (K + K')/2;
    % ·´³ý·¨
    if solver == 1
        U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    else
        [U(freedofs,:),~,~] = pcg(K(freedofs,freedofs),F(freedofs,:),1e-7,5000,diag(diag(K(freedofs,freedofs))));
    end
    %% EVALUATE OPTIMIZATION RESPONSES AND DO SENSITIVITY ANALYSIS
    obj = F'*U;
    ce = sum((U(edofMat)*KE).*U(edofMat),2);
    dcdx = eIntopMat*(-penal*E0*ePhiProj.^(penal-1).*ce.*edproj);
    vol = sum(ePhiProj*elelen^3);
    voldgdx = eIntopMat*(edproj*elelen^3);
    %% SCALE THE OBJECTIVE FUNCTION
    if loop == 1; c_scale = obj/10; end
    dcdx = dcdx/c_scale;
    %% UPDATE DESIGN VARIABLES
    m = 1; 
    cc = 10000*ones(m,1); d = zeros(m,1); a0 = 1; a = zeros(m,1);
    fval = zeros(m, 1); fval(1) = 100*(vol/tolvol-volfra); 
    dfdx = zeros(m, n); dfdx(1,:) = 100*voldgdx/tolvol;
    [xmma,~,~,~,~,~,~,~,~,low,upp]=...
        mmasub(m,n,loop,x,xmin,xmax,xold1,xold2, ...
        obj,dcdx,fval,dfdx,low,upp,a0,a,cc,d);
    xold2 = xold1; xold1 = x; x = xmma;   
    % TUNE PROJECTION PARAMETER
    change = abs(obj-objold)/obj;
    if change < 0.005 && loop > 30
        ichange = ichange+1;
    else
        ichange = 1;
    end
    if mod(ichange,3) == 0
        beta = min(beta * 1.1,20);
    end
    %% PRINT RESULTS AND PLOT DENSITIES
    fprintf([' It.:%5i Obj.:%9.4f Vol:%7.4f numdesvars :%5i' ...
        ' beta:%5.1f ch.:%6.3f\n'],...
        loop,obj,vol/tolvol,neig,beta,change);
    plot_xPhys = reshape(ePhiProj,[nely,nelx,nelz]);
    plot_xPhys = cat(3,flip(plot_xPhys,3),plot_xPhys);
    figure(1),clf
    if displayflag
        display_3D(plot_xPhys);
    else
        plot_final(plot_xPhys);
    end
    figure(2); clf;
    Obj = cat(2,Obj,obj);
    Volf = cat(2,Volf,vol/tolvol);
    plotConvergence(Obj,Volf);    
end
figure(3);
if displayflag
    plot_final(plot_xPhys);
else
    display_3D(plot_xPhys);
end
%% SAVE RESULT
save('result.mat','ePhiProj');
end
