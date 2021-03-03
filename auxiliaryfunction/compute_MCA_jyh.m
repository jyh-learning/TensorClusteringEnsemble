function MCA_ML=compute_MCA_jyh(BP)
% input the base partition, where each column BP is a base partitions
n= size(BP,1); % number of samples
m=size(BP,2); % number of the BPs


MCA_ML=zeros(n);
for i=1:n
    for j=1:n
        if BP(i,:)==BP(j,:)
            MCA_ML(i,j)=1;
        end
    end
end

