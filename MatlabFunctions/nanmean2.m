function means = nanmean2(data,dim)

if nargin == 1
    dim = 1;
end

if dim == 1 %(columns)
    dim_length = length(data(1,:));
    means = zeros(1,dim_length);
    for i = 1:dim_length
        col = data(:,i);
        m = mean(col(~isnan(col)));
        means(i) = m;
    end
elseif dim == 2 %(rows)
    dim_length = length(data(:,1));
    means = zeros(dim_length,1);
    for i = 1:dim_length
        col = data(i,:);
        m = mean(col(~isnan(col)));
        means(i) = m;
    end
else
    disp('Error')
end