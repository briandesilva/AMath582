% creates the a dictionary file using the words found in a set of specified
% text documents
tic
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

% Remove stopwords from dictionary
dictionary = setdiff(dictionary, stop_words,'rows');
celldict = cellstr(dictionary);

clear A stop_words;

% Construct histogram matrix
histmat = zeros(length(dictionary), length(dirListing));


filenames = cell(1,length(dirListing));
for i = 1:length(dirListing)
    filename = fullfile(folder,dirListing(i).name);
    filenames(i) = cellstr(dirListing(i).name);
    
    if ~isempty(regexp(filename,'\.(txt)', 'once'))
        A = sort(parser_cell(filename));
        
        
        index = 1;
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
clear i j k index celldict A count;

filenames( :, ~any(histmat, 1) ) = [];
histmat( :, ~any(histmat,1) ) = [];

toc
