[shmar, simar] = QM([NX; NY]', newtris);

figure(1);
h1 = histogram(simar, 'Normalization','probability');
%xlim([0 10]);
xlabel('Quality measure');
ylabel('Percentage of triangles');

figure(2);
histogram(shmar, 'Normalization', 'probability');
xlim([0 0.5]);
xlabel('Quality measure');
ylabel('Percentage of triangles');
