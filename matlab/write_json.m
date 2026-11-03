function write_json(path,value)
f=fopen(path,'w');assert(f>=0);clean=onCleanup(@()fclose(f));fprintf(f,'%s',jsonencode(value,'PrettyPrint',true));
end
