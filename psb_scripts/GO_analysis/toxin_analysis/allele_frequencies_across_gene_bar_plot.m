allelefrequencies_T = allelefrequencies.';
positions = allelefrequencies(1,:);
frequencies = allelefrequencies_T(:, 2:end);

figure
ba = bar(positions, frequencies, "stacked", 'FaceColor','flat');
ba(1).CData = [0.5 0.5 0.5];  % gray (N)
ba(2).CData = [1 0 0];  % red (A)
ba(3).CData = [0.9290 0.6940 0.1250];  % orange (T)
ba(4).CData = [0 1 0];  % green (C)
ba(5).CData = [1 0 1];  % magenta (G)
ba(6).CData = [0 0 0];  % black (other)