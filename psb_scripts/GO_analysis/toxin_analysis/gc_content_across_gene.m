%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a line plot that will show the 
% GC content across the gene.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create matrices for data on a specific clade.
allelefrequencies_clade_T = allelefrequencies_clade.';
positions_clade = allelefrequencies_clade(1,:);
frequencies_clade = allelefrequencies_clade_T(:, 2:end);

% Frequencies are in order: N, A, T, C, G, other
gc = frequencies_clade(:,4) + frequencies_clade(:,5);
at = frequencies_clade(:,2) + frequencies_clade(:,3);

n = 10;             % Number of elements to create the mean over
s1 = size(gc, 1);      % Find the next smaller multiple of n
m  = s1 - mod(s1, n);
y  = reshape(gc(1:m), n, []);     % Reshape x to a [n, m/n] matrix
avg_gc = transpose(sum(y, 1) / n*100);  % Calculate the mean over the 1st dim

s1 = size(at, 1);      % Find the next smaller multiple of n
m  = s1 - mod(s1, n);
y  = reshape(at(1:m), n, []);     % Reshape x to a [n, m/n] matrix
avg_at = transpose(sum(y, 1) / n*100);  % Calculate the mean over the 1st dim

%%%
% Create matrices for data on a specific clade.
allelefrequencies_all_T = allelefrequencies_all.';
positions_all = allelefrequencies_all(1,:);
frequencies_all = allelefrequencies_all_T(:, 2:end);

% Frequencies are in order: N, A, T, C, G, other
gc_all = frequencies_all(:,4) + frequencies_all(:,5);
at_all = frequencies_all(:,2) + frequencies_all(:,3);

s1 = size(gc_all, 1);      % Find the next smaller multiple of n
m  = s1 - mod(s1, n);
y  = reshape(gc_all(1:m), n, []);     % Reshape x to a [n, m/n] matrix
avg_gc_all = transpose(sum(y, 1) / n*100);  % Calculate the mean over the 1st dim

s1 = size(at_all, 1);      % Find the next smaller multiple of n
m  = s1 - mod(s1, n);
y  = reshape(at_all(1:m), n, []);     % Reshape x to a [n, m/n] matrix
avg_at_all = transpose(sum(y, 1) / n*100);  % Calculate the mean over the 1st dim

legend('GC content', 'AT content');
positions = [0:10:(length(avg_at_all)*10)-1];

% Plot two graphs on the same figure.
figure
tiledlayout(2,1);
ax1 = nexttile;
plot(positions, avg_gc);
hold on
plot(positions, avg_at);
title('WapA 2 GC content across the genome (E clade only)');
legend(ax1,{'GC%','AT%'})

ax2 = nexttile;
plot(positions, avg_gc_all);
hold on
plot(positions, avg_at_all);
title('WapA 2 GC content across the genome (all clades except E clade)');
legend(ax2,{'GC%','AT%'})
