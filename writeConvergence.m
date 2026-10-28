function writeConvergence(hh,eL2,eH1,rL2,rH1,rInf, ...
                          families,famName,rhoCases,rhoName,outDir)
%writeConvergence  Figure and tables of the convergence study
%   writeConvergence(...) writes ex_convergence.pdf
%   into outDir, and prints the rates to the console.
%
%   It is called by ex_convergence.m once the errors have been computed, and
%   can be called again on the contents of ex_convergence.mat to redraw the
%   figure and rewrite the table without repeating the computation.
%
%   AUTHOR: Juan G. Calvo, 2026
nF = numel(families);  nR = numel(rhoCases);

%% figure: log-log plots, one column per mesh family
mk = {'o-','s-','^-','d-','v-'};
figure('Position',[100 100 1100 660])
for iF = 1:nF
    for row = 1:2
        if row == 1, E = eL2; slope = 2; lab = '$\|u-\Pi u_h\|_{0}/\|u\|_{0}$';
        else,        E = eH1; slope = 1; lab = '$|u-\Pi u_h|_{1,h}/|u|_{1}$'; end
        subplot(2,nF,(row-1)*nF+iF); hold on; box on
        for iR = 1:nR
            loglog(squeeze(hh(iF,iR,:)),squeeze(E(iF,iR,:)),mk{iR},'LineWidth',1.1, ...
                   'MarkerSize',5,'MarkerFaceColor','w');
        end
        set(gca,'XScale','log','YScale','log')
        % reference slope
        hR = squeeze(hh(iF,1,:)); eR = squeeze(E(iF,1,:));
        c  = eR(end)/hR(end)^slope;
        loglog(hR,0.35*c*hR.^slope,'k--','LineWidth',0.8)
        text(hR(2),0.35*c*hR(2)^slope,sprintf('  $h^%d$',slope), ...
             'Interpreter','latex','VerticalAlignment','top')
        xlabel('$h$','Interpreter','latex')
        if iF == 1, ylabel(lab,'Interpreter','latex'); end
        if row == 1, title(famName{iF},'Interpreter','latex'); end
        if row == 1 && iF == 1, hLeg = flipud(findobj(gca,'Type','line')); end
        axis tight
    end
end
% one legend for the whole figure, below the plots
lg = legend(hLeg(1:nR),rhoName,'Interpreter','latex','Orientation','horizontal', ...
            'NumColumns',min(nR,3),'FontSize',9);
lg.Units = 'normalized';  lg.Position = [0.5-lg.Position(3)/2 0.005 lg.Position(3:4)];
exportgraphics(gcf,[outDir 'ex_convergence.pdf'],'ContentType','vector')

%% console table
fprintf('\n%s\n','=== estimated convergence rates ===');
fprintf('%-10s %-6s %10s %10s %10s\n','family','rho','rate L2','rate H1','rate inf');
for iF = 1:nF
    for iR = 1:nR
        fprintf('%-10s %-6d %10.2f %10.2f %10.2f\n', ...
                families{iF},rhoCases(iR),rL2(iF,iR),rH1(iF,iR),rInf(iF,iR));
    end
end
end
