function dlgboxInput1()
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
prompt = {'Enter matrix size:','Enter colormap name:'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'20','hsv'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

end

