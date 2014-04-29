%% Script that will try to compute correlation matrix of huge dataset by dividing work into pieces

% So, we have this wholeThing array that has all the things we need.
% For, each column we should compute correlation to other elements, and
% then see what happens.

% We can create matrices of size 20823x200
for i=2:20823
    rowCorrelation = zeros(20824,1);
    rowCorrelation(i) = 1;
   for j=1:20823
      firstColumn = wholeThing(:,i);
      secondColumn = wholeThing(:,j);
      correlation = corrcoef(firstColumn,secondColumn).^8;
      % 0,0 -> i,i of the whole thing. Goes nowhere
      % 0,1 -> i,j. Only use one as the rest is useless.
      % 1,0 -> j,i
      % 1,1 -> j,j. Goes nowhere.
      if (j~=i)
        rowCorrelation(j,1) = correlation(1,2);
      end
   end
   dlmwrite('corResults.csv',transpose(rowCorrelation),'delimiter',',','-append');
   % Row correlation is the thing we have to write into file.
end