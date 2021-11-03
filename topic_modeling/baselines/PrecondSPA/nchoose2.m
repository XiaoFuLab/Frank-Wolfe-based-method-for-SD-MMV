% Create an r-by-r*(r-1)/2 matrix H such that each column has exactly two 
% non-zero entries equal to one, and all columns are different. 

function H = nchoose2(r)

R = nchoosek(r,2);
H = zeros(r,R);

k = 1; 
for i = 1 : r
    for j = i+1:r
        H(i,k) = 1; H(j,k) = 1; 
        k = k + 1; 
    end
end