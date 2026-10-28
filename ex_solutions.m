%% ex_solutions.m
% Companion figures of ex_convergence.m:
%
%   ex_solutions.pdf     the discrete solution u_h on a coarse mesh of each of
%                        the three families, so that the elements are visible
%                        and the three discretizations can be compared, for the
%                        coefficients (i), (iii) and (v); the last column is the
%                        cut of case (v) along y = 1/2, where the kink of the
%                        exact solution at the interface is apparent and the
%                        jumps of the discontinuous solution can be seen
%   ex_coefficients.pdf  the coefficients themselves, as the method sees them,
%                        that is, the value of rho at the centroid of each
%                        element; the constant case (i) is omitted
%
% Cases (i), (ii) and (iv) share the exact solution sin(pi x) sin(pi y), so only
% (i) is shown; (iv) and (v) share the same coefficient, so ex_coefficients.pdf
% has three panels.
%
% AUTHOR: Juan G. Calvo, 2026

clear; close all

%% parameters
coarse   = {'tria_NN_13','hexf_NN_6','voro_NN_12'}; % one coarse mesh per family
famName  = {'Triangles','Hexagons','Voronoi'};
N        = 17;                                      % hexf_NN_N for the coefficients
eta      = 1e1;                                     % contrast
gammah   = 1e1;                                     % penalty parameter gamma_h
outDir   = './';                                    % results are written here
rhoName  = {'(i) $\rho=1$', ...
            sprintf('(ii) $\\rho=1+%gxy$',eta), ...   % the value of eta, not the symbol
            '(iii) $\rho=2+\cos(10\pi x)\cos(10\pi y)$', ...
            '(iv) interface, smooth $u$', ...
            '(v) interface, kink in $u$'};

%% discrete solutions: one row per family of meshes
cases = [0 2 4];                                    % coefficients to be shown
figure('Position',[100 100 1150 900])
tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
for iF = 1:3
    mesh = load(['./meshes/' coarse{iF} '.mat']); mesh = mesh.mesh;
    g = [0; cumsum(mesh.Nverts)];
    % faces and vertices of the discontinuous space: one vertex per degree of
    % freedom, so that the elements are drawn independently of each other
    nv    = mesh.Nverts;
    faces = nan(mesh.NE,max(nv));
    xy    = nan(g(end),2);
    for el = 1:mesh.NE
        faces(el,1:nv(el)) = (g(el)+1):g(el+1);
        xy((g(el)+1):g(el+1),:) = mesh.verts(mesh.elems{el},:);
    end
    for iC = 1:3
        iR = cases(iC);
        [rho,rhs,uex] = coefficient(iR,eta);
        [A,b] = discVem_weights(mesh,rhs,gammah,rho);
        u = A\b;
        ax = nexttile;
        patch('Faces',faces,'Vertices',[xy u],'FaceVertexCData',u, ...
              'FaceColor','interp','EdgeColor',[0.2 0.2 0.2],'LineWidth',0.5)
        view(-40,25); grid on; box on
        xlim([0 1]); ylim([0 1]); zlim([0 max(u)*1.05]); pbaspect([1 1 0.85])
        % the domain is the unit square, so the horizontal ticks carry no
        % information and are dropped: the panels are small once printed
        set(ax,'FontSize',8,'XTickLabel',[],'YTickLabel',[])
        if(iF == 1), title(rhoName{iR+1},'Interpreter','latex','FontSize',9); end
        if(iC == 1)
            text(ax,-0.22,0.5,famName{iF},'Units','normalized','Rotation',90, ...
                 'HorizontalAlignment','center','FontSize',10)
        end
        if(iF == 3 && iC == 1)          % orientation, on one panel only
            xlabel('$x$','Interpreter','latex','FontSize',9)
            ylabel('$y$','Interpreter','latex','FontSize',9)
        end
        if(iR == 4), uK = u; rhoK = rho; uexK = uex; end
    end

    % cut of case (v) along y = 1/2, drawn one element at a time so that the
    % discontinuities of u_h across the elements are visible
    vemOut = vem(mesh,rhs,rhoK.Fn(mesh.bary(:,1),mesh.bary(:,2)));
    ax = nexttile; hold on; box on
    hCut = [];
    for el = 1:mesh.NE
        P  = mesh.verts(mesh.elems{el},:);  np = size(P,1);
        % where the boundary of the element meets the line y = 1/2; vertices
        % lying on it are included, so that an element having an edge on the
        % line also contributes and the jump across that edge becomes visible
        xc = [];
        for i = 1:np
            A1 = P(i,:);  B1 = P(mod(i,np)+1,:);
            dA = A1(2)-0.5;  dB = B1(2)-0.5;
            if(abs(dA) < 1e-12)
                xc(end+1) = A1(1);                              %#ok<AGROW>
            elseif((dA > 0 && dB < -1e-12) || (dA < 0 && dB > 1e-12))
                xc(end+1) = A1(1)+dA/(dA-dB)*(B1(1)-A1(1));     %#ok<AGROW>
            end
        end
        if(numel(xc) < 2 || max(xc)-min(xc) < 1e-12), continue, end
        xc = [min(xc) max(xc)];
        a  = vemOut.projector{el}*uK((g(el)+1):g(el+1));
        c  = vemOut.centroid(el,:);  hE = vemOut.H(el);
        vv = a(1) + a(2)*(xc-c(1))/hE + a(3)*(0.5-c(2))/hE;
        hCut = plot(xc,vv,'-','Color',[0 0.447 0.741],'LineWidth',1.1);
    end
    xs = linspace(0,1,801);
    hEx = plot(xs,uexK(xs,0.5*ones(size(xs))),'k--','LineWidth',0.9);
    xline(0.5,':','Color',[0.5 0.5 0.5])
    set(ax,'FontSize',8); xlim([0 1]); ylim([0 0.56])
    xlabel('$x$','Interpreter','latex','FontSize',9)
    if(iF == 1)
        title('(v) cut along $y=1/2$','Interpreter','latex','FontSize',9)
        legend([hCut hEx],{'$\Pi^{\nabla}u_h$','exact'},'Interpreter','latex', ...
               'Location','southeast','FontSize',7)
    end
    fprintf('%-10s NE = %4d, h = %.3f\n',famName{iF},mesh.NE,max(mesh.h));
end
% rasterized: the elements are drawn as independent patches and a vector
% export leaves visible seams between them
exportgraphics(gcf,[outDir 'ex_solutions.pdf'],'ContentType','image','Resolution',400)

%% plot coefficients, elementwise constant as used by the method
fname = ['./meshes/hexf_NN_' num2str(N) '.mat'];
mesh = load(fname); mesh = mesh.mesh;
g  = [0; cumsum(mesh.Nverts)];
nv = mesh.Nverts;
faces = nan(mesh.NE,max(nv));  xy = nan(g(end),2);
for el = 1:mesh.NE
    faces(el,1:nv(el)) = (g(el)+1):g(el+1);
    xy((g(el)+1):g(el+1),:) = mesh.verts(mesh.elems{el},:);
end
figure('Position',[100 100 1080 320])
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
label = {rhoName{2},rhoName{3},'(iv) and (v)'};
for k = 1:3
    rho = coefficient(k,eta);                      % cases (ii), (iii) and (iv)
    rv  = rho.Fn(mesh.bary(:,1),mesh.bary(:,2));
    nexttile
    patch('Faces',faces,'Vertices',xy,'FaceVertexCData',rv, ...
          'FaceColor','flat','EdgeColor','none')
    axis equal tight; box on
    xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex')
    title(label{k},'Interpreter','latex','FontSize',10)
    colorbar
end
exportgraphics(gcf,[outDir 'ex_coefficients.pdf'],'ContentType','vector')
fprintf('Written to %s: ex_solutions and ex_coefficients, as .pdf\n',outDir);
