function s=paper_file_sha256(path)
% Streaming SHA-256 for reference archives and extracted numerical files.
md=java.security.MessageDigest.getInstance('SHA-256');f=fopen(path,'rb');assert(f>=0);
closeFile=onCleanup(@()fclose(f)); %#ok<NASGU>
while ~feof(f),b=fread(f,8*1024*1024,'*uint8');md.update(b);end
s=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
