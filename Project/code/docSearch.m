% Queries the user for a directory to search then for a series of files to search.
% The file entered must lie in the directory specified or an error will be
% generated.


% Parameter (number of topics to be used in nonnegative matrix
% factorization)
num_topics = 3;


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


% % TF-IDF (Optional)
modified_histmat = tf_idf(histmat);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Apply nonnegative matrix factorization to the histogram matrix
% [U V] = nnmf(histmat,num_topics);             % without tf-idf
[U V] = nnmf(modified_histmat,num_topics);    % with tf-idf
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


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
                break;
            end
        end
        
        % If not, display an error message
        if need_valid_file
            display('Error, file not found. Please enter the name of a file contained in the directory to be searched.');
        end
    end
    
    
    % Compute differences between search file and other files in terms of their
    % representations in terms of the topics using Euclidean norm
    ind = strcmp(filenames,search_file);

    
% Euclidean similarity
% % %     To run search using V matrix from NMF (documents in terms of
% % %     topics)
    dif = bsxfun(@minus,V(:,ind),V) .^ 2;
    
% % %     To run search based on histogram alone (without NMF)
%     dif = bsxfun(@minus,modified_histmat(:,ind),modified_histmat) .^ 2;
   
    distances = sum(dif,1);
clear dif;



% % Cosine similarity
% col_norm = norm(V(:,ind));
% distances = bsxfun(@dot,V,repmat(V(:,ind),1,doc_count));
% distances = distances ./ (sqrt(sum(V .^2, 1)) .* col_norm);
% distances = 1 - distances;
% clear col_norm;

    
    
%     Display results
    [score, indx] = sort(distances,'ascend');
    display('Top 15 results (most to least similar):');
    results = filenames;
    score = 1 - score;
%     score = 1 - (score ./ max(score));
%     score = exp(- (score .^ 2) );
    results(2,indx) = num2cell(score);
    display(results(:,indx(1:15))');
    

%     Display top words in each topic
    display('Top ten words in each topic:');
    display(celldict(top_words(U)));
    
%     Display documents in terms of topics grayscale image

    imshow((1 - V'),'InitialMagnification','fit');
    axis on;
    colormap(gray);
    title('Documents in terms of topics');
    xlabel('Topics');
    ylabel('Documents');
    
    
    prompt2 = 'Would you like to search for another document in the same directory? (y/n)\n';
    response = input(prompt2,'s');
    
    if strcmp(response,'n') || strcmp(response,'N')
        more_words = false;
    end
    
end
clear distances indx results i ind more_words need_valid_file num_topics prompt1 prompt2 response result search_file