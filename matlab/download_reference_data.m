function report=download_reference_data(folder)
% Download the six full-reference archives from the fixed GitHub release.
repo=fileparts(fileparts(mfilename('fullpath')));
M=jsondecode(fileread(fullfile(repo,'data','reference_archives.json')));
if ~exist(folder,'dir'),mkdir(folder);end
downloads=fullfile(folder,'downloads');if ~exist(downloads,'dir'),mkdir(downloads);end
for k=1:numel(M.archives)
 A=M.archives(k);path=fullfile(downloads,A.name);
 if ~exist(path,'file')||~strcmp(paper_file_sha256(path),A.sha256)
  fprintf('Downloading %s (%.1f MB)\n',A.name,A.bytes/1e6);
  websave(path,A.url,weboptions('Timeout',3600));
 end
 assert(strcmp(paper_file_sha256(path),A.sha256),'Archive SHA-256 mismatch: %s',A.name);
 unzip(path,folder);
 for j=1:numel(A.files)
  p=fullfile(folder,A.files(j).path);assert(exist(p,'file')==2&&strcmp(paper_file_sha256(p),A.files(j).sha256));
 end
end
report=struct('tag',M.tag,'directory',folder,'archives',numel(M.archives),'allChecksumsPass',true);
fprintf('REFERENCE_DATA_DOWNLOAD_PASS\n');
end
