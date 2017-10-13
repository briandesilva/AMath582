% This script computes the histogram matrix of a set of text documents then
% the nonnegative matrix factorization of the histogram. Next it outputs
% the top 10 words from each of the num_topics topics, where num_topics
% is a parameter.



% Instructions:
% 1. Set the parameter num_topics to the number of topics you would like to compute
% 2. Run the script
% 3. Select the folder containing the text documents to be analyzed from
%    the folder selection window (Note: The script will read ALL text files
%    from the folder you select)

num_topics = 3;




dictionary = [];
% queries the user for the folder containing the txt files
folder = uigetdir;
dirListing = dir(folder);
doc_count = 0;

fprintf('Constructing dictionary...')
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
disp('done.')


load('stop_words.mat');

% Remove stop words
dictionary = setdiff(dictionary, stop_words,'rows');
celldict = cellstr(dictionary);


clear A stop_words;

% Initialize histogram matrix
histmat = zeros(length(dictionary), length(dirListing));


fprintf('Constructing histogram matrix...')
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
disp('done.')

clear i j k index A count filename dirListing;

% Remove hidden documents containing no words from histogram
filenames( :, ~any(histmat, 1) ) = [];
histmat( :, ~any(histmat, 1) ) = [];


% Run nonnegative matrix factorization on the histogram matrix
[U, ~] = nnmf(histmat,num_topics);

% Find and display the top ten words from each topic
[~, IX] = sort(U,'descend');
top_indices = IX(1:10,:);

display('Top ten words from each topic:');
top_ten = celldict(top_indices);
display(top_ten);

