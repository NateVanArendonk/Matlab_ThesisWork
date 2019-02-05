% For help see http://matlaboratory.blogspot.co.uk/2013/08/alternative-function-nanmean-multi.html

function m=nanmean3(data,dim)

if nargin == 1
    dim = 1;
end

dims=ndims(data);

%Permute data so target dimension is dimension 1
if dim>1
    dim_vec=1:dims;
    target_dims=dim:-1:1;
    dim_vec(1:length(target_dims))=target_dims;
    
    r_data=permute(data,[dim_vec]);
else 
    r_data=data;
end

%Perform nanmean
r_data_size=size(r_data);
r_data_lin=r_data(:);

step=length(r_data(:,1));
ind=[0-(step-1), 0];
range=[1:length(r_data_lin)/length(r_data(:,1))];
means=zeros(1,max(range));
for i= range
    ind=ind+step;
        
    col=r_data_lin(ind,1);
    m = mean(col(~isnan(col)));
    means(i) = m;
end
means=reshape(means,[1, r_data_size(2:end)]);

%Permute back so target dimension goes back where it came from
if dim>1
    m=permute(means,dim_vec);
else
    m=means;
end