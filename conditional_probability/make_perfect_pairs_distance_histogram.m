% Create a distance histograms based on perfect pairs.
% Choose pairs of loci that comply to the definition of perfectly vertically descendent loci on a tree.
% Calculate the distance between the loci in a pair and plot as a histogram.
% Make sure to normalize the results to the overall data and then plot in on a log-log scale

% Import the raw counts dataset.
A = importdata('conditional_probability_dataset_raw_counts.txt');

A = importdata('pairs_dataset_raw_counts_F_clade_only_4plus_extended.txt');
A = A.data;

% Find perfect pairs
length(find(A(:,3)==min(A(:,4),A(:,5))))
s=find(A(:,3)==min(A(:,4),A(:,5)));
% Calculate distance between perfect pairs
r=A(s,2)-A(s,1);

% Plot the distance as a histogram
figure(10);hist(r,0:200:50000)
nh=hist(r,0:200:50000);
figure(11);plot(log(100+(0:200:50000)),log(nh)) % Log-log hist graph
hold on;x=5:1:11;plot(x,7.5-.46*(x-5),'r-') % Best fit line

% Check the distance distribution of all the pairs in the raw data
r0=A(:,2)-A(:,1);
nh0=hist(r0,0:200:50000); % Plot the histogram

% Plot all vs perfect distance distribution
plot((0:200:50000),nh0,'b')
hold on; plot((0:200:50000),nh, 'g')

% Normalize the perfect distribution on the all pairs density
plot((0:200:50000), nh./nh0, 'b')
hold on; plot(log(100+(0:200:50000)), log(nh./nh0), 'r') % log-log normalized density
hold on; x=5:1:11;plot(x,-1.2-0.25*(x-5),'b-') % Best fit line (slope=-0.25)
