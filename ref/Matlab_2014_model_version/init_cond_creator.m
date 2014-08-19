function init_cond_creator(filename, newfilename)
  data = load(filename);
  
  fid = fopen(newfilename, 'w');

style = '';
for n = 1:(29 +6+ 4*2) % functions, differential variables and time
    style = [style ' %6.6e'];
end

fprintf(fid, style, [data(end, 2:44)]); 

fclose(fid);
