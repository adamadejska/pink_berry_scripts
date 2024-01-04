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
xregion(3362633,3366313);
xline(3366313);
xregion(3366315,3367145);
xline(3366315);xline(3367145);
xregion(3368130,3368339);
xline(3368130);xline(3368339);
xregion(3368740,3370062);
xline(3368740);xline(3370062);
xregion(3370053,3370376);
xline(3370053);xline(3370376);
xregion(3370712,3372031);
xline(3370712); xline(3372031);
xregion(3372010,3373149);
xline(3372010);xline(3373149);
xregion(3373426,3375078);
title('WapA 2 sample depths across the genome (E clade only)');

nexttile;
plot(positions_all, frequencies_all);
xregion(3362633,3366313);
xline(3366313);
xregion(3366315,3367145);
xline(3366315);xline(3367145);
xregion(3368130,3368339);
xline(3368130);xline(3368339);
xregion(3368740,3370062);
xline(3368740);xline(3370062);
xregion(3370053,3370376);
xline(3370053);xline(3370376);
xregion(3370712,3372031);
xline(3370712); xline(3372031);
xregion(3372010,3373149);
xline(3372010);xline(3373149);
xregion(3373426,3375078);
title('WapA 2 sample depths across the genome (all clades except E clade)');

% Compute a mean for all samples
M = mean(frequencies_all,2);
M_c = mean(frequencies_clade,2);
nexttile;
plot(positions_all, M);
hold on;
plot(positions_all, M_c);
xregion(3362633,3366313);
xline(3366313);
xregion(3366315,3367145);
xline(3366315);xline(3367145);
xregion(3368130,3368339);
xline(3368130);xline(3368339);
xregion(3368740,3370062);
xline(3368740);xline(3370062);
xregion(3370053,3370376);
xline(3370053);xline(3370376);
xregion(3370712,3372031);
xline(3370712); xline(3372031);
xregion(3372010,3373149);
xline(3372010);xline(3373149);
xregion(3373426,3375078);
title('WapA 2 sample depths across the genome averaged');