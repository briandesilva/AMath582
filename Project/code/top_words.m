function result = top_words(U, n)
%     returns a matrix, result, with the top n words of each topic in
%     descending order
%     Each column corresponds to a topic
%     U is the basis matrix returned by NMF

[~, IX] = sort(U,'descend');

result = IX(1:n,:);

end