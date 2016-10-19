

dname = 'U:\Narain\Microtubules\140204\'; % folder with images

fnames = dir([dname '*.ht3']); 

% irf = load('d:\MatlabCalc QY\Data\R6G for mix\irf.txt'); % irf
for j=3:numel(fnames)
    name = fnames(j).name;
    Process_file([dname name], 1);
    
%     FastFLIM([dname name], [1.5 3])    
end

for j=1:numel(fnames)
    name=fnames(j).name;
    ViewFLIM([dname name])
end

