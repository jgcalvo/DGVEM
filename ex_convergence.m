%% ex_convergence.m
% Convergence test for the weighted SIPG virtual element discretization of
%       -div(rho grad u) = f   in (0,1)^2,      u = 0   on the boundary,
% with exact solution u(x,y) = sin(pi x) sin(pi y).
%
% Three families of polygonal meshes are used (triangles, hexagons and
% Voronoi cells) and four coefficients rho:
%    0: rho = 1                                       (constant)
%    1: rho = 1 + eta*x*y                             (smooth, large gradient)
%    2: rho = 2 + cos(10 pi x) cos(10 pi y)           (oscillatory)
%    3: rho = 1 if x < 1/2, rho = eta otherwise       (discontinuous)
% Case 3 is an interface problem for which the solution above is still exact,
% since d u/d x vanishes on x = 1/2 and hence the normal flux rho du/dx is
% continuous across the interface. For that same reason the jump across the
% interface is almost inactive, and case 4 is provided as a sharper test: its
% exact solution has a genuine kink at x = 1/2, with a nonzero normal flux.
%    4: rho as in case 3, u = phi(x) sin(pi y) with phi linear on the left and
%       quadratic on the right, chosen so that phi and rho phi' are continuous
% The three families of meshes are fitted to the interface, so that rho is
% constant on each element as assumed by the theory: the structured triangular
% meshes already are, make_voronoi mirrors the column of generating points next
% to the interface, and make_hexagonal cuts the cells that the line crosses.
%
% Three relative errors are reported. The first two are computed from the
% elliptic projection Pi u_h of the discrete solution, which is the only
% pointwise information available for a virtual element function:
%    e_L2  = |u - Pi u_h|_{L2} / |u|_{L2}                 (expected O(h^2))
%    e_H1  = |grad u - grad Pi u_h|_{L2} / |grad u|_{L2}  (expected O(h))
% The third one compares degrees of freedom with the nodal interpolant u_I,
%    e_inf = |u_I - u_h|_inf / |u_I|_inf                  (expected O(h^2))
%
% Output: ex_convergence.pdf (log-log plots)
% written to this folder; ex_convergence.mat keeps the raw data
%
% AUTHOR: Juan G. Calvo, 2026

clear; close all

%% parameters
families = {'tria','hexa','voro'};                  % mesh families
famName  = {'Triangles','Hexagons','Voronoi'};
levels   = { [41 65 97 161 257], ...                % tria_NN_*
             [12 17 24  34  48], ...                % hexf_NN_*
             [24 34 48  68  96] };                  % voro_NN_*
rhoCases = 0:4;                                     % coefficients (drop any of
                                                    % them when reporting results)
eta      = 1e1;                                     % contrast
rhoName  = {'(i) $\rho=1$', ...
            sprintf('(ii) $\\rho=1+%gxy$',eta), ...  % the value of eta
            '(iii) $\rho=2+\cos(10\pi x)\cos(10\pi y)$', ...
            '(iv) interface, smooth $u$', ...
            '(v) interface, kink in $u$'};
rhoName  = rhoName(rhoCases+1);
outDir   = './';                                    % results are written here
gammah   = 1e1;                                     % penalty parameter gamma_h

nF = numel(families); nR = numel(rhoCases); nL = numel(levels{1});
hh   = nan(nF,nR,nL);   eL2 = nan(nF,nR,nL);
eH1  = nan(nF,nR,nL);   einf= nan(nF,nR,nL);

%% loop over mesh families, coefficients and refinement levels
for iF = 1:nF
    for iR = 1:nR
        [rho,rhs,uex,gux,guy] = coefficient(rhoCases(iR),eta);
        for iL = 1:nL
            mesh = loadMesh(families{iF},levels{iF}(iL));

            % discretization
            [A,b] = discVem_weights(mesh,rhs,gammah,rho);
            u     = A\b;

            % projector, centroids and diameters (needed for the L2/H1 errors)
            vemOut = vem(mesh,rhs,rho.Fn(mesh.bary(:,1),mesh.bary(:,2)));

            % nodal interpolant of the exact solution (one dof per element/vertex)
            g  = [0; cumsum(mesh.Nverts)];
            uI = nan(g(end),1);
            for j = 1:mesh.NE
                d = (g(j)+1):g(j+1);
                uI(d) = uex(mesh.verts(mesh.elems{j},1),mesh.verts(mesh.elems{j},2));
            end

            [l2,h1] = projErrors(mesh,vemOut,u,uex,gux,guy);
            hh(iF,iR,iL)   = max(mesh.h);
            eL2(iF,iR,iL)  = l2;
            eH1(iF,iR,iL)  = h1;
            einf(iF,iR,iL) = norm(uI-u,Inf)/norm(uI,Inf);
            fprintf('%s  rho=%d  h=%.4f  NE=%6d  e_L2=%.3e  e_H1=%.3e  e_inf=%.3e\n', ...
                    families{iF},rhoCases(iR),hh(iF,iR,iL),mesh.NE, ...
                    l2,h1,einf(iF,iR,iL));
        end
    end
end

%% convergence rates (least squares fit of log(err) against log(h))
rL2 = nan(nF,nR); rH1 = nan(nF,nR); rInf = nan(nF,nR);
for iF = 1:nF
    for iR = 1:nR
        L = log(squeeze(hh(iF,iR,:)));
        p = polyfit(L,log(squeeze(eL2 (iF,iR,:))),1); rL2 (iF,iR) = p(1);
        p = polyfit(L,log(squeeze(eH1 (iF,iR,:))),1); rH1 (iF,iR) = p(1);
        p = polyfit(L,log(squeeze(einf(iF,iR,:))),1); rInf(iF,iR) = p(1);
    end
end

writeConvergence(hh,eL2,eH1,rL2,rH1,rInf, ...
                 families,famName,rhoCases,rhoName,outDir);

save('ex_convergence.mat','hh','eL2','eH1','einf','rL2','rH1','rInf', ...
     'families','famName','levels','rhoCases','rhoName','eta','gammah')
fprintf('\nWritten to %s: ex_convergence.pdf\n',outDir);

%% ---------------------------------------------------------------- helpers
function [eL2,eH1] = projErrors(mesh,vemOut,u,uex,gux,guy)
% projErrors computes the relative L2 and broken H1 errors between the exact
% solution and the elliptic projection Pi u_h, which is affine on each element.
% Each polygon is split into triangles from its centroid and integrated with a
% six point rule of degree four.
bar = [0.108103018168070 0.445948490915965 0.445948490915965
       0.445948490915965 0.108103018168070 0.445948490915965
       0.445948490915965 0.445948490915965 0.108103018168070
       0.816847572980459 0.091576213509771 0.091576213509771
       0.091576213509771 0.816847572980459 0.091576213509771
       0.091576213509771 0.091576213509771 0.816847572980459];
wq  = [0.223381589678011*ones(3,1); 0.109951743655322*ones(3,1)];

g = [0; cumsum(mesh.Nverts)];
nL2 = 0; dL2 = 0; nH1 = 0; dH1 = 0;
for el = 1:mesh.NE
    V  = mesh.verts(mesh.elems{el},:);  nv = size(V,1);
    c  = vemOut.centroid(el,:);  hE = vemOut.H(el);
    a  = vemOut.projector{el}*u((g(el)+1):g(el+1));   % coefficients in the
    gP = [a(2) a(3)]/hE;                              % scaled monomial basis
    cb = mean(V,1);                                   % centroid used to split
    for s = 1:nv
        T = [cb; V(s,:); V(mod(s,nv)+1,:)];
        ar = abs(det([T(2,:)-T(1,:); T(3,:)-T(1,:)]))/2;
        P  = bar*T;                                   % quadrature points
        uq = uex(P(:,1),P(:,2));
        pq = a(1) + a(2)*(P(:,1)-c(1))/hE + a(3)*(P(:,2)-c(2))/hE;
        gx = gux(P(:,1),P(:,2));  gy = guy(P(:,1),P(:,2));
        nL2 = nL2 + ar*sum(wq.*(uq-pq).^2);
        dL2 = dL2 + ar*sum(wq.*uq.^2);
        nH1 = nH1 + ar*sum(wq.*((gx-gP(1)).^2+(gy-gP(2)).^2));
        dH1 = dH1 + ar*sum(wq.*(gx.^2+gy.^2));
    end
end
eL2 = sqrt(nL2/dL2);  eH1 = sqrt(nH1/dH1);
end

function mesh = loadMesh(family,N)
% loadMesh reads a mesh of the given family, adding the fields required by
% the virtual element routines when they are not already stored
switch family
    case 'tria', fname = ['./meshes/tria_NN_' num2str(N) '.mat'];
    case 'hexa', fname = ['./meshes/hexf_NN_' num2str(N) '.mat'];
    case 'voro', fname = ['./meshes/voro_NN_' num2str(N) '.mat'];
end
mesh = load(fname); mesh = mesh.mesh;
end
