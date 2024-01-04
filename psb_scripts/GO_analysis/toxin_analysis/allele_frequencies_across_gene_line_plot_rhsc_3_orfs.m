%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a line plot that will show the 
% depths of each sample across the gene.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create matrices for data on all clades.
allele_depths_all_T = alleledepths_all.';
positions_all = alleledepths_all(1,:);
frequencies_all = allele_depths_all_T(:, 2:end);

% Create matrices for data on a specific clade.
alleledepths_clade_T = alleledepths_clade.';
positions_clade = alleledepths_clade(1,:);
frequencies_clade = alleledepths_clade_T(:, 2:end);

% Plot two bar graphs on the same figure.
figure
tiledlayout(3,1);
nexttile;
plot(positions_clade, frequencies_clade);
xregion(3350557,3354828);
xline(3354828);
xregion(3354828,3355520);
xline(3354828);xline(3355520);
xregion(3355609,3356232);
xline(3355609);xline(3356232);
xregion(3356424,3356621);
xline(3356424);xline(3356621);
xregion(3356914,3357006);
xline(3356914);xline(3357006);
xregion(3357003,3357290);
xline(3357003); xline(3357290);
xregion(3357428,3357604);
xline(3357428); xline(3357604);
xregion(3358062,3358802);
xline(3358062); xline(3358802);
xregion(3358855,3359484);
xline(3358855); xline(3359484);

title('rhsC 3 sample depths across the genome (E clade only)');

nexttile;
plot(positions_all, frequencies_all);
xregion(3350557,3354828);
xline(3354828);
xregion(3354828,3355520);
xline(3354828);xline(3355520);
xregion(3355609,3356232);
xline(3355609);xline(3356232);
xregion(3356424,3356621);
xline(3356424);xline(3356621);
xregion(3356914,3357006);
xline(3356914);xline(3357006);
xregion(3357003,3357290);
xline(3357003); xline(3357290);
xregion(3357428,3357604);
xline(3357428); xline(3357604);
xregion(3358062,3358802);
xline(3358062); xline(3358802);
xregion(3358855,3359484);
xline(3358855); xline(3359484);
title('rhsC 3 sample depths across the genome (all clades except E clade)');

% Compute a mean for all samples
M = mean(frequencies_all,2);
M_c = mean(frequencies_clade,2);
nexttile;
plot(positions_all, M);
hold on;
plot(positions_all, M_c);
xregion(3350557,3354828);
xline(3354828);
xregion(3354828,3355520);
xline(3354828);xline(3355520);
xregion(3355609,3356232);
xline(3355609);xline(3356232);
xregion(3356424,3356621);
xline(3356424);xline(3356621);
xregion(3356914,3357006);
xline(3356914);xline(3357006);
xregion(3357003,3357290);
xline(3357003); xline(3357290);
xregion(3357428,3357604);
xline(3357428); xline(3357604);
xregion(3358062,3358802);
xline(3358062); xline(3358802);
xregion(3358855,3359484);
xline(3358855); xline(3359484);

title('rhsC 3 sample depths across the genome averaged');