function [press] = loadPressure(pressname)

disp('Loading file...');

fid = fopen(pressname);
if fopen(fid) == -1
    error('Pressure file not found.');
end

press = load(pressname,'-ascii');
disp('File loaded.');   

end

