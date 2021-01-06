function out = poly2d_projection(im,order,varargin)
inputs = inputParser;
inputs.addRequired('im',@isnumeric);
inputs.addRequired('order',@isnumeric);
inputs.addParameter('msk',true(size(im)),@islogical);
inputs.parse(im,order,varargin{:});

msk = inputs.Results.msk;

nmodes=((order+1)^2-(order+1))/2+(order+1);

%build the polynomial basis and the basis spatial covariance matrix
%over the image mask pixels

jnm=zeros(3,nmodes);
sz=size(im);
xyrt=getGridCoordinates(sz(2),sz(1),1);
basis=zeros(sz(1),sz(2),nmodes);
j_mode=0;
for n=0:order
    for k=0:n
        j_mode=j_mode+1;
        basis(:,:,j_mode) = (xyrt.x2D.^(n-k)) .* (xyrt.y2D.^k);
        jnm(:,j_mode) = [j_mode,n-k,k];
    end
end

Gij=zeros(nmodes,nmodes);
for k=1:nmodes
    for j=1:k
        Gij(k,j)=sum(sum(basis(:,:,k).*basis(:,:,j).*msk));
        Gij(j,k)=Gij(k,j);
    end
end


%get modes coefficients
  b=zeros(1,nmodes);
  for k=1:nmodes
      b(k)=sum(sum(im.*basis(:,:,k).*msk));
  end
  a=pinv(Gij)*b';


  %build the model from the modes and the coefficients
  model=0;
  for k=1:nmodes 
      model=model+a(k)*basis(:,:,k);
  end
  
  %results
  out.a = a;
  out.b = b;
  out.Gij = Gij;
  out.basis = basis;
  out.nmodes = nmodes;
  out.jnm = jnm;
  out.model = model;
  
    