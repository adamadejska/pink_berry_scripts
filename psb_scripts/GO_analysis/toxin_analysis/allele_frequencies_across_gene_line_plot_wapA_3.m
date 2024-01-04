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
xregion(5049505,5058900);
xline(5058900);
xregion(5046919,5049462);
xline(5046919);xline(5049462);
xregion(5042095,5046771);
xline(5042095);xline(5046771);
xregion(5041002,5041919);
xline(5041002);xline(5041919);
xregion(5039083,5041005);
xline(5039083);xline(5041005);
xregion(5037116,5039074);
xline(5037116); xline(5039074);

title('WapA 3 sample depths across the genome (8 clade only)');

nexttile;
plot(positions_all, frequencies_all);
xregion(5049505,5058900);
xline(5058900);
xregion(5046919,5049462);
xline(5046919);xline(5049462);
xregion(5042095,5046771);
xline(5042095);xline(5046771);
xregion(5041002,5041919);
xline(5041002);xline(5041919);
xregion(5039083,5041005);
xline(5039083);xline(5041005);
xregion(5037116,5039074);
xline(5037116); xline(5039074);

title('WapA 3 sample depths across the genome (all clades except 8 clade)');

% Compute a mean for all samples
M = mean(frequencies_all,2);
M_c = mean(frequencies_clade,2);
nexttile;
plot(positions_all, M);
hold on;
plot(positions_all, M_c);
xregion(5049505,5058900);
xline(5058900);
xregion(5046919,5049462);
xline(5046919);xline(5049462);
xregion(5042095,5046771);
xline(5042095);xline(5046771);
xregion(5041002,5041919);
xline(5041002);xline(5041919);
xregion(5039083,5041005);
xline(5039083);xline(5041005);
xregion(5037116,5039074);
xline(5037116); xline(5039074);

title('WapA 3 sample depths across the genome averaged');