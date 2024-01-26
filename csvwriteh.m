function csvwriteh(nameout,header,data)

    fout=openw(nameout,'wt');
    fprintf(fout,'%s,',header{1:end-1});
    fprintf(fout,'%s\n',header{end});
    fclose(fout);
    csvwrite(nameout,data,0,0,'-append');
    return
