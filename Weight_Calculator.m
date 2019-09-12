filename='payload_interaction_weights.xlsx';
n=4; % Number of weighted criteria

comp=xlsread(filename,'Decision Matrix');

%Matrix normalized with column sums (Shinpaugh Lecture 5)
normal=zeros(n);
for i=1:n
    normal(:,i)=comp(:,i)./sum(comp(:,i));
end

%Consistency indexes for random matrix (from Google)
ri=[0 0 0.58 0.9 1.12 1.24 1.32 1.41 1.45 1.49 1.51 1.48 1.56 1.57 1.59];
ri=ri(n);

%Calculating the weight for each criterion
weights=zeros(1,n);
for i=1:n
    weights(i)=mean(normal(i,:));
end

%Multiplying each column in normal by weights
comp2=zeros(n);
for i=1:n
    comp2(:,i)=comp(:,i).*weights(i);
end

%Vector of the sum of the rows of comp2
wsum=sum(comp2');

%Dividing wsum by weights
v=wsum./weights;

%Calculating lambda max
lambda=mean(v);

%Calculating consistency index
ci=(lambda-n)/(n-1);

%Consistency ratio
cr=ci/ri;
if cr <= 0.1
    fprintf("\nComparisons are Consistent\n")
else
    fprintf("\nComparisons are NOT Consistent\n")
end

%Exporting the weights to the excel file
xlswrite(filename,weights','Weights',['B2:B' num2str(n+1)])

