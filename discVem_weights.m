function [AD,bD,A] = discVem_weights(mesh,f,alpha,varargin)
% read fields of mesh
verts = mesh.verts;
elems = mesh.elems;
edges2elems = mesh.edges2elems;
edges       = mesh.edges;
bdNode      = mesh.bdNode;
NE    = mesh.NE;
Ne    = mesh.Ne;
NeBnd = mesh.NeBnd;
NeInt = mesh.NeInt;
nVerts= mesh.Nverts;

% coefficient $\rho$
if(size(varargin,2)==1)  % read $\rho$ (if given)
    rho = varargin{1};
    isCont = rho.isCont;
    rho    = rho.Fn;
    if(isa(rho,'function_handle')) % evaluate \rho if function is given
        % evaluation of \rho at centroids of elements
        centroid = nan(NE,2);
        for el_id = 1:NE
            vert_ids = elems{el_id};      % Global IDs of the vertices of this element
            locVert = verts(vert_ids, :); % Coordinates of the vertices of this element
            area_components = locVert(:,1) .* locVert([2:end,1],2) - locVert([2:end,1],1) .* locVert(:,2);
            area = 0.5 * abs(sum(area_components));
            centroid(el_id,:) = sum((locVert + locVert([2:end,1],:)) .* repmat(area_components,1,2)) / (6*area);
        end
        rhoEl = rho(centroid(:,1),centroid(:,2)); % eval at incenter
        if(isCont) 
           % evaluation of \rho at midpoints of edges
            midptsEdges = (verts(edges(:,1),:)+verts(edges(:,2),:))/2;
            rhoEdge = rho(midptsEdges(:,1),midptsEdges(:,2));
        end
    else                                % vector with one entry per element
        rhoEl = rho;
        isCont = 0;
    end
else                                    % default: \rho = 1
    rhoEl = ones(NE,1);
    isCont = 1;
end

% contributions of \int_E \rho \nabla \pi\phi_i \cdot \nabla \pi\phi_j
vemOut = vem(mesh,f,rhoEl);
projec = vemOut.projector;
diamEl = vemOut.H;
A1     = vemOut.K;
bD     = vemOut.f;
clear vemOut

% contributions of edge integrals
% allocate vectors
sizeVecB = 0;
for edgeID = 1:Ne                    % loop over all edges
    if(numel(edges2elems{edgeID})==2)
        sizeVecB = sizeVecB + 4*sum(nVerts(edges2elems{edgeID}));
    else
        sizeVecB = sizeVecB + 2*sum(nVerts(edges2elems{edgeID}));
    end
end
sizeVecC = 16*NeInt+4*NeBnd;
iiB = nan(sizeVecB,1); jjB = nan(sizeVecB,1); ssB = nan(sizeVecB,1); ptB = 0;% allocate space
iiC = nan(sizeVecC,1); jjC = nan(sizeVecC,1); ssC = nan(sizeVecC,1); ptC = 0;% for matrices
% global dofs
globDof = [0; cumsum(mesh.Nverts)]; %
% contributions of edge integrals
for edgeID = 1:Ne                     % loop over all edges
    currEdge = edges(edgeID,:);
    elem1  = edges2elems{edgeID};     % read element K+, K-
    normal = verts(currEdge(2),:)-verts(currEdge(1),:); % normal (no unitaria)
    normal = [normal(2), -normal(1)];
    if(numel(elem1)==2)              % interior edge
        weight = rhoEl(elem1(1))*rhoEl(elem1(2))/(rhoEl(elem1(1))+rhoEl(elem1(2)));
        e1 = elems{elem1(1)};
        e2 = elems{elem1(2)};
        inde1 = find(e1==currEdge(1));
        inde2 = find(e2==currEdge(1));
        if(e1(mod(inde1(1),numel(e1))+1) ~= currEdge(2))
            tem = e1; % intercambio K+ y K-
            e1 = e2;
            e2 = tem;
            tem = inde1;
            inde1 = inde2;
            inde2 = tem;
            elem1 = flipud(elem1);
        end
        dofe1 = (globDof(elem1(1))+1):globDof(elem1(1)+1);
        dofe2 = (globDof(elem1(2))+1):globDof(elem1(2)+1);
        reorder1 = [inde1:numel(e1), 1:inde1-1];
        dofe1 = dofe1(reorder1);
        reorder2 = [inde2:-1:1, numel(e2):-1:(inde2+1)];
        dofe2 = dofe2(reorder2);
        % bilinear form b1
        pi1 = projec{elem1(1)}; pi1 = pi1(:,reorder1);
        pi2 = projec{elem1(2)}; pi2 = pi2(:,reorder2);
        v1 = (normal*pi1(2:3,:)/(2*diamEl(elem1(1))))';
        v2 = (normal*pi2(2:3,:)/(2*diamEl(elem1(2))))';
        if(isCont) % cont. rho
            M1 = (rhoEdge(edgeID)/2)*[v1 v1 -v1 -v1; v2 v2 -v2 -v2]; %[B11 -B11; -B22 B22];
        else       % disc. rho
            M1 = weight*[v1 v1 -v1 -v1; ...
                         v2 v2 -v2 -v2];
        end        
        adv = numel(M1);
        [yy,xx] = meshgrid([dofe1(1:2) dofe2(1:2)],[dofe1 dofe2]);
        iiB(ptB+1:ptB+adv) = xx(:);
        jjB(ptB+1:ptB+adv) = yy(:);
        ssB(ptB+1:ptB+adv) = M1(:);
        ptB = ptB + adv;
        % bilinear form b2
        B11 = [1/3 1/6; 1/6 1/3];
        M2  = 2*weight*[B11 -B11; -B11 B11];
        adv = 16;
        [yy,xx] = meshgrid([dofe1(1:2) dofe2(1:2)]);
        iiC(ptC+1:ptC+adv) = xx(:);
        jjC(ptC+1:ptC+adv) = yy(:);
        ssC(ptC+1:ptC+adv) = M2(:);
        ptC = ptC + adv;
    else                                           % boundary edge
        e1 = elems{elem1(1)};
        cen = mean(verts(e1,:),1);               
        if(dot(cen-verts(currEdge(1),:),normal)>0)
            normal = -normal;
        end
        inde1 = find(e1==currEdge(1));
        dofe1 = (globDof(elem1(1))+1):globDof(elem1(1)+1);
        reorder1 = [inde1:numel(e1), 1:inde1-1];
        e1 = e1(reorder1);
        dofe1 = dofe1(reorder1);
        if(e1(2)~=currEdge(2))            
            dofe1 = dofe1([end 1:end-1]); 
            reorder1 = reorder1([end 1:end-1]);
        end                               
        % bilinear form b1
        pi1 = projec{elem1(1)}; pi1 = pi1(:,reorder1);
        v1 = (normal*pi1(2:3,:)/(2*diamEl(elem1(1))))';
        if(isCont) % cont. rho
            M1 = [v1 v1]*rhoEdge(edgeID);
        else          % disc. rho
            M1 = [v1 v1]*rhoEl(elem1(1));
        end
        adv = numel(M1);
        [yy,xx] = meshgrid(dofe1(1:2),dofe1);
        iiB(ptB+1:ptB+adv) = xx(:);
        jjB(ptB+1:ptB+adv) = yy(:);
        ssB(ptB+1:ptB+adv) = M1(:);
        ptB = ptB + adv;
        % bilinear form b2
        M2 = rhoEl(elem1(1))*[1/3 1/6; 1/6 1/3];
        adv = 4;
        [yy,xx] = meshgrid(dofe1(1:2));
        iiC(ptC+1:ptC+adv) = xx(:);
        jjC(ptC+1:ptC+adv) = yy(:);
        ssC(ptC+1:ptC+adv) = M2(:);
        ptC = ptC + adv;
    end
end
M1 = sparse(iiB,jjB,ssB);
M2 = sparse(iiC,jjC,ssC);

% assemble matrix
A = A1-M1-M1'+alpha*M2;
clear A1 M1 M2

% boundary conditions
Ndof  = globDof(end);
% find boundary dof
bdryDof = [];
for j = 1:NE
    nodes = elems{j};
    tem = find(ismember(nodes,bdNode));
    if(numel(tem)>0)
        bdryDof = [bdryDof; tem+globDof(j)]; %#ok<AGROW>
    end
end
bdidx = zeros(Ndof,1); bdidx(bdryDof) = 1;
Tbd = spdiags(bdidx,0,Ndof,Ndof);
T   = spdiags(1-bdidx,0,Ndof,Ndof);
AD  = T*A*T + Tbd;
bD(bdryDof) = 0;