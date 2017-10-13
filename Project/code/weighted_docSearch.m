% Queries the user for a directory to search then for a series of files to search.
% The file entered must lie in the directory specified or an error will be
% generated.
% Compares different methods side-by-side (spectral clustering, nmf, gram
% matrix, weighted adjacency matrix)

close all
% % % % % PARAMETERS % % % % %

% Parameter (number of topics to be used in nonnegative matrix
% factorization)
num_topics = 10;

% Parameter (constant used to compute adjacency matrix for weight matrix)
sigm = .3;

% Parameter (number of search results to display)
num_display = 10;

% Parameter (indicates whether or not to use cosine similarity to constuct
% adjacency matrix for reweighting of V
use_cos_adj = false;

% Parameter (indicates whether or not to use tf-idf)
use_tfidf = true;

% Parameter (indicates which corpus to load)
% Choices are: {all, bio, biochem, chem, hist}
corpus = 'all';

% Parameter (indicates how many clusters are to be used by k-means)
clusters = 6;

% Parameter (indicates which similarity measure to use for clustering)
% Options:
%             sqeuclidean - Squared Euclidean distance
%             cityblock - L1 distance
%             cosine - 1 minus cosine of the angle betweent the points (treated as vectors)
%             correlation - 1 minus the correlation between points
dist_type = 'cosine';

% Parameter (kmeans will run num_trials times and return the clustering
% which gave the best residual)
num_trials = 10;

% Parameter (nnmf will run num_facts times and return the clustering
% which gave the best residual)
num_facts = 30;

% % % % % % % % % % % % % % % % % % % % % %

if exist('corpus','var')
    load([corpus '_hist.mat'])
    doc_count = length(filenames);
end

if ~exist('histmat','var') 
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
    
    
    % Remove stop words
    load('stop_words.mat');
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


% Apply nonnegative matrix factorization (and possibly tf-idf) to the histogram matrix
if use_tfidf
    modified_histmat = tf_idf(histmat);
    [U V] = nnmf(modified_histmat,num_topics,'algorithm','mult','replicates',num_facts);
else
    [U V] = nnmf(histmat,num_topics,'algorithm','mult','replicates',num_facts);
end


% % Gram matrix modification
gram = U.' * U / (norm(U)^2);
V_gram = gram * V;
clear gram;


% % Adjacency matrix modification
if use_cos_adj
    A = squareform(pdist((U ./ norm(U))','cosine'));
else
    A = squareform(pdist((U ./ norm(U))'));
end
W = exp( - (A .^2) ./ sigm);
V_adj = W * V;
clear W A;


% Display top words in each topic
display('Top ten words in each topic:');
display(celldict(top_words(U,10)));


% Display documents in terms of topics grayscale image

figure();

subplot(1,3,1);
imshow((1 - V'),'InitialMagnification','fit');
axis on;
colormap(gray);
title('Topic Vectors');
xlabel('Topics');
ylabel('Documents');

subplot(1,3,2);
imshow((1 - V_adj'),'InitialMagnification','fit');
axis on;
colormap(gray);
title('Adjacency Matrix Topic Vectors');
xlabel('Topics');
ylabel('Documents');

subplot(1,3,3);
imshow((1 - V_gram'),'InitialMagnification','fit');
axis on;
colormap(gray);
title('Gram Matrix Topic Vectors');
xlabel('Topics');
ylabel('Documents');


% % Horizontal plot for poster
% figure()
% imshow(1-V,'InitialMagnification','fit');
% axis on;
% colormap(gray);
% title('Topic Vectors','FontSize',24)
% xlabel('Documents','FontSize',18)
% ylabel('Topics','FontSize',18)

% % Gram vectors
figure()
imshow(1-V_gram,'InitialMagnification','fit');
axis on;
colormap(gray);
title('Gram Topic Vectors','FontSize',24)
xlabel('Documents','FontSize',18)
ylabel('Topics','FontSize',18)


% Cluster/classify the documents using k-means

% Topic vectors
idx = kmeans(V.',clusters,'Distance',dist_type,'Replicates',num_trials);

% Adjacency matrix reweighting
idx_adj = kmeans(V_adj.',clusters,'Distance',dist_type,'Replicates',num_trials);

% Gram reweighting
idx_gram = kmeans(V_gram.',clusters,'Distance',dist_type,'Replicates',num_trials);


% % Bar plots of clusterings
% figure()
% subplot(1,3,1)
% bar(idx)
% xlabel('Document')
% ylabel('Cluster')
% title('Topic Vectors')
% 
% subplot(1,3,2)
% bar(idx_adj)
% xlabel('Document')
% ylabel('Cluster')
% title('Adjacency Reweighting')
% 
% subplot(1,3,3)
% bar(idx_gram)
% xlabel('Document')
% ylabel('Cluster')
% title('Gram Reweighting')


% Print performance of clusterings
ground_truth = ones(size(idx));
for k=2:(doc_count/5)
    ground_truth((5*(k-1)+1):5*k) = k;
end


[px, ~] = purity(ground_truth, idx);
fprintf('Percentage of correct classifications (topic vectors): ');
fprintf([num2str(px) '\n\n'])

[px_adj, ~] = purity(ground_truth, idx_adj);
fprintf('Percentage of correct classifications (adjacency reweighted topic vectors): ');
fprintf([num2str(px_adj) '\n\n'])

[px_gram, ~] = purity(ground_truth, idx_gram);
fprintf('Percentage of correct classifications (gram reweighted topic vectors): ');
fprintf([num2str(px_gram) '\n\n'])




clear distances indx results i ind more_words need_valid_file num_topics prompt1 prompt2 response result search_file