% This script creates a combined plot of the errors w.r.t. epsilon for all
% of the systems considered. Run only when the systems have all been solved
% and the necessary matrices containing the errors have been saved.

f = figure
f.Position = [1, 1, 1120, 840]
% fontsize(f, 10, "points")
tiledlayout(2, 2, 'TileSpacing','compact');

nexttile

three = load("3state.mat");

epsl_range = three.epsl_range;
normErr0 = three.normErr0;
normErr1 = three.normErr1;
normErr2 = three.normErr2;

loglog(epsl_range, normErr0, 'ob', 'MarkerFaceColor', 'b');
hold on
loglog(epsl_range, normErr1, 'sr', 'MarkerFaceColor', 'r');
loglog(epsl_range, normErr2, 'dk', 'MarkerFaceColor', 'k');
plot(epsl_range,0.0005*epsl_range.^1,'-b')
plot(epsl_range,0.0007*epsl_range.^2,'-r')
plot(epsl_range,0.0024*epsl_range.^3,'-k')
hold off
if length(epsl_range) > 1
    set(gca,'xlim',[min(epsl_range) max(epsl_range)]);
end
xlabel('\epsilon');
ylabel('Error');
title('3-state system')

nexttile

eight = load("8state.mat");

epsl_range = eight.epsl_range;
normErr0 = eight.normErr0;
normErr1 = eight.normErr1;
normErr2 = eight.normErr2;

a = loglog(epsl_range, normErr0, 'ob', 'MarkerFaceColor', 'b');
hold on
b = loglog(epsl_range, normErr1, 'sr', 'MarkerFaceColor', 'r');
c = loglog(epsl_range, normErr2, 'dk', 'MarkerFaceColor', 'k');
plot(epsl_range, 0.21*epsl_range.^1, '-b');
plot(epsl_range, 0.17*epsl_range.^2, '-r');
plot(epsl_range, 0.28*epsl_range.^3, '-k');
hold off
if length(epsl_range) > 1
    set(gca,'xlim',[min(epsl_range) max(epsl_range)]);
end
xlabel('\epsilon');
ylabel('Error');
title('8-state system');
hleglines = [a, b, c];
hleg = legend(hleglines, 'O(1) Appr. Err.',...
    'O(\epsilon) Appr. Err.',...
    'O(\epsilon^2) Appr. Err.');
set(hleg,'Location','SouthEast');

nexttile

one = load("1Dlat.mat");

epsl_range = one.epsl_range;
normErr0 = one.normErr0;
normErr1 = one.normErr1;
normErr2 = one.normErr2;

loglog(epsl_range, normErr0, 'ob', 'MarkerFaceColor', 'b');
hold on
loglog(epsl_range, normErr1, 'sr', 'MarkerFaceColor', 'r');
loglog(epsl_range, normErr2, 'dk', 'MarkerFaceColor', 'k');
plot(epsl_range,0.56*epsl_range.^1,'-b')
plot(epsl_range,18.36*epsl_range.^2,'-r')
plot(epsl_range,707.1*epsl_range.^3,'-k')
hold off
if length(epsl_range) > 1
    set(gca,'xlim',[min(epsl_range) max(epsl_range)]);
end
xlabel('\epsilon');
ylabel('Error');
title('1D lattice model')

nexttile

two = load("2Dlat.mat");

epsl_range = two.epsl_range;
normErr0 = two.normErr0;
normErr1 = two.normErr1;
normErr2 = two.normErr2;

loglog(epsl_range, normErr0, 'ob', 'MarkerFaceColor', 'b');
hold on
loglog(epsl_range, normErr1, 'sr', 'MarkerFaceColor', 'r');
loglog(epsl_range, normErr2, 'dk', 'MarkerFaceColor', 'k');
d = plot(epsl_range,0.20*epsl_range.^1,'-b')
e = plot(epsl_range,1.14*epsl_range.^2,'-r')
g = plot(epsl_range,6.27*epsl_range.^3,'-k')
hold off
if length(epsl_range) > 1
    set(gca,'xlim',[min(epsl_range) max(epsl_range)]);
end
xlabel('\epsilon');
ylabel('Error');
title('2D lattice model');
hleglines = [d, e, g];
hleg = legend(hleglines, 'C\times\epsilon',...
    'C\times\epsilon^2',...
    'C\times\epsilon^3');
set(hleg,'Location','SouthEast');

fontsize(f, scale=1.5);

exportgraphics(f, 'Combined.pdf')