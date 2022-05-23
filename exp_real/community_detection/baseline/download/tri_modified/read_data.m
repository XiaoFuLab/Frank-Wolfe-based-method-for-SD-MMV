clear


%% read DBLP

adjacency = cell(7,1);
community = cell(7,1);
for d = 1:5
    
    name = ['DBLP',num2str(d),'_adjacency.txt'];
    file = fopen(name);
    cont = textscan( file, '%f %f' );
    
    row = cont{1}; col = cont{2};
    n = max([row;col]);
    adjacency{d} = sparse(row,col,1,n,n);
    
    name = ['DBLP',num2str(d),'_community.txt'];
    file = fopen(name);
    cont = textscan( file, '%f %f %f' );
    
    row = cont{1}; col = cont{2}; val = cont{3};
    n = max(row); k = max(col);
    community{d} = full(sparse(row,col,val));
    
end

%% read MAG

for d = 1:2
    
    name = ['MAG',num2str(d),'_adjacency.txt'];
    file = fopen(name);
    cont = textscan( file, '%f %f' );
    
    row = cont{1}; col = cont{2};
    n = max([row;col]);
    adjacency{d+5} = sparse(row,col,1,n,n);
    
    name = ['MAG',num2str(d),'_community.txt'];
    file = fopen(name);
    cont = textscan( file, '%f %f %f' );
    
    row = cont{1}; col = cont{2}; val = cont{3};
    n = max(row); k = max(col);
    community{d+5} = full(sparse(row,col,val));
    
end

save('networks.mat','adjacency','community')