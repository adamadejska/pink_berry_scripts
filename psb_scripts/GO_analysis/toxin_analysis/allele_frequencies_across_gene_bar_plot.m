%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a stacked bar graph that will show the 
% allele frequencies across the gene.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create matrices for data on all clades.
allelefrequencies_all_T = allelefrequencies_all.';
positions_all = allelefrequencies_all(1,:);
frequencies_all = allelefrequencies_all_T(:, 2:end);

% Create matrices for data on a specific clade.
allelefrequencies_clade_T = allelefrequencies_clade.';
positions_clade = allelefrequencies_clade(1,:);
frequencies_clade = allelefrequencies_clade_T(:, 2:end);

% Plot two bar graphs on the same figure.
figure
tiledlayout(2,1);
nexttile;
ba = bar(positions_clade, frequencies_clade, "stacked", 'FaceColor','flat');
title('WapA 3 allele frequencies across the genome (8 clade only)');
ba(1).CData = [0.5 0.5 0.5];  % gray (N)
ba(2).CData = [1 0 0];  % red (A)
ba(3).CData = [0.9290 0.6940 0.1250];  % orange (T)
ba(4).CData = [0 1 0];  % green (C)
ba(5).CData = [1 0 1];  % magenta (G)
ba(6).CData = [0 0 0];  % black (other)

nexttile;
ba2 = bar(positions_all, frequencies_all, "stacked", 'FaceColor','flat');
ylim([0 1]);
title('WapA 3 allele frequencies across the genome (all clades except the 8 clade)');
ba2(1).CData = [0.5 0.5 0.5];  % gray (N)
ba2(2).CData = [1 0 0];  % red (A)
ba2(3).CData = [0.9290 0.6940 0.1250];  % orange (T)
ba2(4).CData = [0 1 0];  % green (C)
ba2(5).CData = [1 0 1];  % magenta (G)
ba2(6).CData = [0 0 0];  % black (other)