function lms=ntr2lms(ntr, illum, cones, wavelengths)
% converts n-transition functions (one per row, defined by transition
% wavelengths) to lms using the given cone sensitivities and illumination
% type I: first wavelength smaller than last, else assume type II
finish = wavelengths(end);
start = wavelengths(1);
nWavelengths=size(wavelengths,1);

% create interpolation functions
nCones=size(cones,2);
illcones=cones.*repmat(illum,1,nCones);
cumcones=cumtrapz(illcones); % cumulative cones
cumcones=repmat(cumcones(end,:),nWavelengths,1)-cumcones;

ppcumcones=spline(wavelengths',cumcones'); % interpolated
getcum=@(WLs) ppval(ppcumcones, WLs); % get cumulative response
maxcum=getcum(start); % maximal response (white) needed for type 2

cum2lmsT1=@(v) sum([v(:,1:2:end) -v(:,2:2:end)],2);  % Type I
cum2lmsT2=@(v) maxcum-cum2lmsT1(v); % Type II = white - type I
lmsntrT1=@(WLs) cum2lmsT1(getcum(WLs))'; % Type I sensor response
lmsntrT2=@(WLs) cum2lmsT2(getcum(WLs))'; % Type II sensor response

% get lms...
lms=zeros(size(ntr,1),nCones);
for i=1:size(ntr,1)
    if ntr(i,1)<=ntr(i,end)
        lms(i,:)=lmsntrT1(ntr(i,:));
%         display([num2str(i) ' : Type 1']);
    else
        lms(i,:)=lmsntrT2(ntr(i,end:-1:1));
%         display([num2str(i) ' : Type 2']);
    end
end