% % % Iterative docSearch % % %
% Uses identical code to docSearch, but chooses the number of topics to use
% dynamically. Iterates over different numbers of topics (2 through k)
% multiple times, averaging the distance between the search document and
% other documents in the specified directory. It then returns the documents
% with the lowest average distances from the search document.


% Queries the user for a directory to search then for a series of files to search.
% The file entered must lie in the directory specified or an error will be
% generated.


% Parameters
max_topics = 30;
num_iterations = 15;


% Parameter (constant used to compute adjacency matrix for weight matrix)
sigm = .3;

% Parameter (number of search results to display)
num_display = 20;


if ~(exist('histmat','var') == 1)
    dictionary = [];
    
    % queries the user for the folder containing the txt files
    folder = uigetdir;
    dirListing = dir(folder);
    doc_count = 0;
   
    % constructs the dictionary
    for i = 1:length(dirListing)
        filename = fullfile(folder,dirListing(i).name);
        
        if ~isempty(regexp(filename,'\.(txt)', 'once'))
            A = parser(filename);
            
            dictionary = char(dictionary,A);
            dictionary = unique(dictionary,'rows');
            doc_count = doc_count + 1;
            
        end
        
    end
    
    
    load('stop_words.mat');
    
    % Remove stop words
    dictionary = setdiff(dictionary, stop_words,'rows');
    celldict = cellstr(dictionary);
    
    
    clear A stop_words;
    
    % Initialize histogram matrix
    histmat = zeros(length(dictionary), length(dirListing));
    
    
    % Construct histogram matrix (histmat)
    filenames = cell(1,length(dirListing));
    for i = 1:length(dirListing)
        filename = fullfile(folder,dirListing(i).name);
        filenames(i) = cellstr(dirListing(i).name);
        
        if ~isempty(regexp(filename,'\.(txt)', 'once'))
            A = sort(parser_cell(filename));
            
            
            index = 1;
            % cycle through each word in the dictionary
            for j = 1:length(dictionary)
                count = 0;
                for k = index:length(A);
                    if cellfun(@strcmp,celldict(j),A(k))
                        count = count + 1;
                        index = k+1;
                    end
                    
                    if lexcmp(celldict(j),A(k))==-1
                        break;
                    end
                    
                end
                histmat(j,i) = count;
            end
        end
    end
    
    clear i j k index A count filename dirListing;
    
    % Remove hidden documents containing no words from histogram
    filenames( :, ~any(histmat, 1) ) = [];
    histmat( :, ~any(histmat, 1) ) = [];
    
end

% % TF-IDF (Optional)
modified_histmat = tf_idf(histmat);
 

% Get the file from the user and make sure it exists in the specified
% directory
more_words = true;
while more_words
    
    need_valid_file = true;
    while(need_valid_file)
        % Queries the user to enter a file name
        prompt1 = 'Please enter the search document, (e.g. file1.txt):\n';
        search_file = input(prompt1,'s');
        
        
        % Check whether or not the file entered by the user lies in the specified
        % folder
        for i = 1:size(filenames,2)
            if strcmp(search_file,filenames(i))
                need_valid_file = false;
            end
        end
        
        % If not, display an error message
        if need_valid_file
            display('Error, file not found. Please enter the name of a file contained in the directory to be searched.');
        end
    end
    
    
    distances = zeros(1,doc_count);
    distances_gram = zeros(1,doc_count);
    distances_mod = zeros(1,doc_count);
    distances_eigs = zeros(1,doc_count);
    ind = strcmp(filenames,search_file);
    
    for dummyVar = 1:num_iterations
        for num_topics = 3:max_topics
            
            
            % % TF-IDF (Optional)
            % modified_histmat = tf_idf(histmat);
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % Apply nonnegative matrix factorization to the histogram matrix
            
            % [U V] = nnmf(histmat,num_topics);             % without tf-idf
            [U V] = nnmf(modified_histmat,num_topics);    % with tf-idf
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            
            
            % % Gram matrix modification
            gram = (U ./ norm(U))' * (U ./ norm(U));
            V_gram = gram * V;
            clear gram;
            
            
            % % Adjacency matrix modification
            A = squareform(pdist((U ./ norm(U))'));
            W = exp( - (A .^2) ./ sigm);
            V_mod = W * V;
            clear W A;
            
            
            % % Spectral clustering on V
            [~, e] = size(V);
            D1 = squareform(pdist(V'));
            sig = 1;
            A1 = exp( -(D1.^2)./sig^2 ) - eye(e);
            %    A1 = 1-D1;
            %A1 = bsxfun(@rdivide,A1,sum(A1,2));
            [V_eigs E1] = eig(A1);
            E1 = diag(E1);
            [E1, IX]= sort(E1, 'descend');
            
            V_eigs = V_eigs(1:end, IX);
            [~, num_evecs] = max(abs(gradient(E1(2:end))));
            V_eigs = V_eigs(:,2:num_evecs + 1);
            V_eigs = V_eigs';
            clear num_evecs D1 A1 E1 IX;
            
            
            
            % Compute differences between search file and other files in terms of their
            % representations in terms of the topics using Euclidean norm
            dif = bsxfun(@minus,V(:,ind),V) .^ 2;
            
            % % %     To run search based on histogram alone (without NMF), use the
            % % %     following line:
            %     dif = bsxfun(@minus,modified_histmat(:,ind),modified_histmat) .^ 2;
            
            distances = distances + sqrt(sum(dif,1));
            clear dif;
         
            
            
            % % Similarity using adjacency matrix altered topics
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % % %     To run search using V matrix from NMF (documents in terms of
            % % %     topics)
            dif_mod = bsxfun(@minus,V_mod(:,ind),V_mod) .^ 2;
            
            distances_mod = distances_mod + sqrt(sum(dif_mod,1));
            clear dif_mod;

            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            
            % % Similarity using gram matrix altered topics
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % % %     To run search using V matrix from NMF (documents in terms of
            % % %     topics)
            dif_gram = bsxfun(@minus,V_gram(:,ind),V_gram) .^ 2;
            
            distances_gram = distances_gram + sqrt(sum(dif_gram,1));
            clear dif_gram;
            
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            
            
            % % Similarity using eigenvectors from spectral clustering applied to V
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % % %     To run search using V matrix from NMF (documents in terms of
            % % %     topics)
            dif_eigs = bsxfun(@minus,V_eigs(:,ind),V_eigs) .^ 2;
            
            distances_eigs = distances_eigs + sqrt(sum(dif_eigs,1));
            clear dif_eigs;
            
        end
    end
    
    
%     Display search results using NMF
    [score, indx] = sort(distances,'ascend');
    display('Top 15 results (most to least similar, using topic vectors):');
    results = filenames;
    score = 1 - (score ./ num_iterations);
%     score = 1 - (score ./ max(score));
%     score = exp(- (score .^ 2) );
    results(2,indx) = num2cell(score);
    display(results(:,indx(1:num_display))');
    
    
%     Display search results from adjacency matrix altered topics
    [score, indx] = sort(distances_mod,'ascend');
    display('Top 15 results (most to least similar, using adjacency matrix):');
    results = filenames;
    score = 1 - (score ./ num_iterations);
    results(2,indx) = num2cell(score);
    display(results(:,indx(1:num_display))');
    
    
%     Display search results from gram matrix altered topics
    [score, indx] = sort(distances_gram,'ascend');
    display('Top 15 results (most to least similar, using gram matrix):');
    results = filenames;
    score = 1 - (score ./ num_iterations);
    results(2,indx) = num2cell(score);
    display(results(:,indx(1:num_display))');
    
    
%     Display search results from eigenvectors of spectral clustering
    [score, indx] = sort(distances_eigs,'ascend');
    display('Top 15 results (most to least similar, using eigenvectors of V):');
    results = filenames;
    score = 1 - (score ./ num_iterations);
    results(2,indx) = num2cell(score);
    display(results(:,indx(1:num_display))');
    
    

%     Display top words in each topic
    display('Top ten words in each topic:');
    display(celldict(top_words(U,10)));
    
    
    
    
    prompt2 = 'Would you like to search for another document in the same directory? (y/n)\n';
    response = input(prompt2,'s');
    
    if strcmp(response,'n') || strcmp(response,'N')
        more_words = false;
    end
    
end
clear distances indx results i ind more_words need_valid_file num_topics prompt1 prompt2 response result search_file