function vemOut = vem(mesh,rhs,rhoEl)
% standard vem
verts = mesh.verts;
elems = mesh.elems;
n_polys = 3;
NT = length(elems);
globDof = [0; cumsum(mesh.Nverts)];
n_dofs = globDof(end);
% allocate outputs
proj = cell(NT,1);
ar   = zeros(NT,1);
diam = zeros(NT,1);
cent = zeros(NT,2);
% indices for sparse matrices
nVerts = cellfun(@numel,elems);    % vertices per element
numEntries = sum(nVerts.^2);
ii = nan(numEntries,1);
jj = nan(numEntries,1);
ss = nan(numEntries,1);
pt = 0;
F = zeros(n_dofs, 1); % Forcing vector
linear_polynomials = {[0,0], [1,0], [0,1]}; % Impose an ordering on the linear polynomials
mod_wrap = @(x, a) mod(x-1, a) + 1; % A utility function for wrapping around a vector
for el_id = 1:NT
    dofel = (globDof(el_id)+1):globDof(el_id+1);
    vert_ids = elems{el_id}; % Global IDs of the vertices of this element
    locVert = verts(vert_ids, :); % Coordinates of the vertices of this element
    n_sides = length(vert_ids); % Start computing the geometric information
    area_components = locVert(:,1) .* locVert([2:end,1],2) - locVert([2:end,1],1) .* locVert(:,2);
    area = 0.5 * abs(sum(area_components));
    centroid = sum((locVert + locVert([2:end,1],:)) .* repmat(area_components,1,2)) / (6*area);
    diameter = 0; % Compute the diameter by looking at every pair of vertices
    for i = 1:(n_sides-1)
        for j = (i+1):n_sides
            diameter = max(diameter, norm(locVert(i, :) - locVert(j, :)));
        end
    end
    D = zeros(n_sides, n_polys); D(:, 1) = 1;
    B = zeros(n_polys, n_sides); B(1, :) = 1/n_sides;
    for vertex_id = 1:n_sides
        vert = locVert(vertex_id, :); % This vertex and its neighbours
        prev = locVert(mod_wrap(vertex_id - 1, n_sides), :);
        next = locVert(mod_wrap(vertex_id + 1, n_sides), :);
        vertex_normal = [next(2) - prev(2), prev(1) - next(1)]; % Average of normals on edges
        for poly_id = 2:n_polys % Only need to loop over non-constant polynomials
            poly_degree = linear_polynomials{poly_id};
            monomial_grad = poly_degree / diameter; % Gradient of a linear polynomial is constant
            D(vertex_id, poly_id) = dot((vert - centroid), poly_degree) / diameter;
            B(poly_id, vertex_id) = 0.5 * dot(monomial_grad, vertex_normal);
        end
    end
    projector = (B*D) \ B; % Compute the local Ritz projector to polynomials
    stabilising_term = (eye(n_sides) - D * projector)' * (eye(n_sides) - D * projector);
    G = B*D; G(1, :) = 0;
    local_stiffness = (projector' * G * projector + stabilising_term)*rhoEl(el_id); 
    adv = numel(local_stiffness);
    [xx,yy] = meshgrid(dofel);
    ii(pt+1:pt+adv) = xx(:);
    jj(pt+1:pt+adv) = yy(:);
    ss(pt+1:pt+adv) = local_stiffness(:);
    pt = pt + adv;
    F(dofel) = F(dofel) + rhs(centroid) * area / n_sides;
    cent(el_id,:) = centroid;
    diam(el_id)   = diameter;
    ar(el_id)     = area;
    proj{el_id}   = projector;
end
K = sparse(ii,jj,ss,n_dofs,n_dofs);
vemOut = struct('projector',{proj},'centroid',cent,'area',ar,'K',K,'H',diam,'f',F);
end