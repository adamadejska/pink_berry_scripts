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
title('RhsC 4 sample depths across the genome (8 clade only)');

nexttile;
plot(positions_all, frequencies_all);
title('RhsC 4 sample depths across the genome (all clades except 8 clade)');

% Compute a mean for all samples
M = mean(frequencies_all,2);
M_c = mean(frequencies_clade,2);
nexttile;
plot(positions_all, M);
hold on;
plot(positions_all, M_c);
title('RhsC 4 sample depths across the genome averaged');