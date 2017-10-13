% tf_idf: takes in a matrix histmat, which has word histograms of documents
% as its columns, ie the i,jth entry of histmat corresponds to the number of
% times the ith word in your dictionary appears in the jth document of your 
% corpus.
% Returns a modified histogram matrix with the term frequency-inverse document 
% frequency of each word in each document replacing the raw frequency.
function modified_histmat = tf_idf(histmat)

%     Initialize modified_histmat with adjusted term frequencies
    max_freq = max(histmat);
%     modified_histmat = bsxfun(@rdivide,histmat .* (.5),max_freq) + (.5);
    modified_histmat = bsxfun(@rdivide,histmat,max_freq);
    
%     Compute inverse document frequency (idf)
    num_occurrences = sum(logical(histmat),2);
    idf = log(size(histmat,2) ./ num_occurrences);
    
    
%     Compute tf-idf
    modified_histmat = bsxfun(@times,modified_histmat,idf);









end