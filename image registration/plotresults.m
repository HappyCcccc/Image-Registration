
h1=histogram(costofs,16);
%h1 = bar(costofs,16);
hold on
h2 = bar(1.04, 1 ,0.1,'r');
hold on
h3 = bar(5.98, 1 ,0.1,'y');
legend([h1 h2 h3]);
%applyhatch(gcf,'\+/');


function applyhatch(h,patterns,colorlist)

%  By Ben Hinkle, bhinkle@mathworks.com

%  This code is in the public domain. 

  
oldppmode = get(h,'paperpositionmode');
oldunits = get(h,'units');
set(h,'paperpositionmode','auto');
set(h,'units','pixels');
figsize = get(h,'position');
if nargin == 2
  colorlist = [];
end
bits = hardcopy(h,'-dzbuffer','-r0');
set(h,'paperpositionmode',oldppmode);

bwidth = size(bits,2);
bheight = size(bits,1);
bsize = bwidth * bheight;
if ~isempty(colorlist)
  colorlist = uint8(255*colorlist);
  [colors,colori] = nextnonbw(0,colorlist,bits);
else
  colors = (bits(:,:,1) ~= bits(:,:,2)) | ...
	   (bits(:,:,1) ~= bits(:,:,3));
end
pati = 1;
colorind = find(colors);
while ~isempty(colorind)
  colorval(1) = bits(colorind(1));
  colorval(2) = bits(colorind(1)+bsize);
  colorval(3) = bits(colorind(1)+2*bsize);
  if iscell(patterns)
    pattern = patterns{pati};
  elseif isa(patterns,'char')
    pattern = makehatch(patterns(pati));
  else
    pattern = patterns;
  end
  pattern = uint8(255*(1-pattern));
  pheight = size(pattern,2);
  pwidth = size(pattern,1);
  ratioh = ceil(bheight/pheight);
  ratiow = ceil(bwidth/pwidth);
  bigpattern = repmat(pattern,[ratioh ratiow]);
  if ratioh*pheight > bheight
    bigpattern(bheight+1:end,:) = [];
  end
  if ratiow*pwidth > bwidth
    bigpattern(:,bwidth+1:end) = [];
  end
  bigpattern = repmat(bigpattern,[1 1 3]);
  color = (bits(:,:,1) == colorval(1)) & ...
	  (bits(:,:,2) == colorval(2)) & ...
	  (bits(:,:,3) == colorval(3));
  color = repmat(color,[1 1 3]);
  bits(color) = bigpattern(color);
  if ~isempty(colorlist)
    [colors,colori] = nextnonbw(colori,colorlist,bits);
  else
    colors = (bits(:,:,1) ~= bits(:,:,2)) | ...
	     (bits(:,:,1) ~= bits(:,:,3));
  end
  colorind = find(colors);
  pati = (pati + 1);
  if pati > length(patterns)
    pati = 1;
  end
end

newfig = figure('units','pixels','visible','off');
imaxes = axes('parent',newfig,'units','pixels');
im = image(bits,'parent',imaxes);
fpos = get(newfig,'position');
set(newfig,'position',[fpos(1:2) figsize(3) figsize(4)+1]);
set(imaxes,'position',[0 0 figsize(3) figsize(4)+1],'visible','off');
set(newfig,'visible','on');
end

function [colors,out] = nextnonbw(ind,colorlist,bits)
out = ind+1;
colors = [];
while out <= size(colorlist,1)
  if isequal(colorlist(out,:),[255 255 255]) | ...
	isequal(colorlist(out,:),[0 0 0])
    out = out+1;
  else
    colors = (colorlist(out,1) == bits(:,:,1)) & ...
	     (colorlist(out,2) == bits(:,:,2)) & ...
	     (colorlist(out,3) == bits(:,:,3));
    return
  end
end
end

function A = makehatch(hatch)

%  By Ben Hinkle, bhinkle@mathworks.com

%  This code is in the public domain. 


n = 6;
A=zeros(n);
switch (hatch)
 case '/'
  A = fliplr(eye(n));
 case '\'
  A = eye(n);
 case '|'
  A(:,1) = 1;
 case '-'
  A(1,:) = 1;
 case '+'
  A(:,1) = 1;
  A(1,:) = 1;
 case 'x'
  A = eye(n) | fliplr(diag(ones(n-1,1),-1));
 case '.'
  A(1:2,1:2)=1;
 otherwise
  error(['Undefined hatch pattern "' hatch '".']);
end
end