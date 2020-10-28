
% Load in tree distance values as a 1xM double object using options when
% loading the data in. (numeric matrix)
% Load in  tree names as a 1xN string array object
% Create a tree object using Neighbour Joining
njtree = seqneighjoin(distances16,'equivar',cluster16names)
% Set the root to the outgroup (sequence of only WT positions)
sel = getbyname(njtree, 'outgroup');
njtree2 = reroot(njtree,sel)
