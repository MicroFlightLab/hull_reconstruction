function [Result,maph] = countRepRow(M)
% count number of repeated rows
[Mu,ia,ic] = unique(M, 'rows', 'stable');           % Unique Values By Row, Retaining Original Order
h = accumarray(ic, 1);                              % Count Occurrences
maph = h(ic);                                       % Map Occurrences To ‘ic’ Values
Result = [Mu, h];

end

