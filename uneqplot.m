% This script creates a plot of the errors w.r.t. epsilon for the 8-state
% theoretical system, but using different initial vectors. Only run this
% when the relevant systems have been run and the matrices containing the
% errors have been saved.

f = figure
f.Position = [1, 1, 1080, 1960]
% fontsize(f, 10, "points")
tiledlayout(2, 2, 'TileSpacing','compact');

nexttile

three = load("8state.mat");

epsl_range = three.epsl_range;
normErr0 = three.normErr0;
normErr1 = three.normErr1;
normErr2 = three.normErr2;

loglog(epsl_range, normErr0, 'ob', 'MarkerFaceColor', 'b');
hold on
loglog(epsl_range, normErr1, 'sr', 'MarkerFaceColor', 'r');
loglog(epsl_range, normErr2, 'dk', 'MarkerFaceColor', 'k');
plot(epsl_range, 0.21*epsl_range.^1, '-b');
plot(epsl_range, 0.17*epsl_range.^2, '-r');
plot(epsl_range, 0.28*epsl_range.^3, '-k');
hold off
if length(epsl_range) > 1
    set(gca,'xlim',[min(epsl_range) max(epsl_range)]);
end
xlabel('\epsilon');
ylabel('Error');
title('With O(\epsilon^2) p_{ini}')

nexttile

eight = load("8stateeqo1.mat");

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
title('With O(\epsilon) p_{ini}');
hleglines = [a, b, c];
hleg = legend('O(1) Appr. Err.',...
    'O(\epsilon) Appr. Err.',...
    'O(\epsilon^2) Appr. Err.');
set(hleg,'Location','SouthEast');

nexttile

one = load("8stateuneq.mat");

epsl_range = one.epsl_range;
normErr0 = one.normErr0;
normErr1 = one.normErr1;
normErr2 = one.normErr2;

loglog(epsl_range, normErr0, 'ob', 'MarkerFaceColor', 'b');
hold on
loglog(epsl_range, normErr1, 'sr', 'MarkerFaceColor', 'r');
loglog(epsl_range, normErr2, 'dk', 'MarkerFaceColor', 'k');
plot(epsl_range, 0.21*epsl_range.^1, '-b');
plot(epsl_range, 0.17*epsl_range.^2, '-r');
plot(epsl_range, 0.28*epsl_range.^3, '-k');
hold off
if length(epsl_range) > 1
    set(gca,'xlim',[min(epsl_range) max(epsl_range)]);
end
xlabel('\epsilon');
ylabel('Error');
title('With O(1) p_{ini}')

nexttile

two = load("8statenoneq.mat");

epsl_range = two.epsl_range;
normErr0 = two.normErr0;
normErr1 = two.normErr1;
normErr2 = two.normErr2;

loglog(epsl_range, normErr0, 'ob', 'MarkerFaceColor', 'b');
hold on
loglog(epsl_range, normErr1, 'sr', 'MarkerFaceColor', 'r');
loglog(epsl_range, normErr2, 'dk', 'MarkerFaceColor', 'k');
d = plot(epsl_range, 0.21*epsl_range.^1, '-b');
e = plot(epsl_range, 0.17*epsl_range.^2, '-r');
g = plot(epsl_range, 0.28*epsl_range.^3, '-k');
hold off
if length(epsl_range) > 1
    set(gca,'xlim',[min(epsl_range) max(epsl_range)]);
end
xlabel('\epsilon');
ylabel('Error');
title('With no quasi-equilibration');
hleglines = [d, e, g];
hleg = legend(hleglines, 'C\times\epsilon',...
    'C\times\epsilon^2',...
    'C\times\epsilon^3');
set(hleg, 'Location', 'SouthEast')
fontsize(f, scale=1.5);

exportgraphics(f, 'noneq.pdf')