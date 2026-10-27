function report = diagnose_connectivity(source, outDir)
%DIAGNOSE_CONNECTIVITY Independent postprocessing of the saved physical graph.
%   D = diagnose_connectivity reads results/experiment_data.mat. It never
%   integrates, edits the state, changes dc, or thresholds positive weights.
%   Optional arguments select a MAT input and a separate output directory.
%   Run audit_connectivity.py for the EPS plot and 80-digit spot checks.
root = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(source), source = fullfile(root,'results','experiment_data.mat'); end
if nargin < 2 || isempty(outDir), outDir = fullfile(root,'connectivity_check'); end
if ~exist(outDir,'dir'), mkdir(outDir); end
tmp = load(source,'result'); R = tmp.result; P = R.P;
assert(P.N == 5 && P.dc == 3.2 && P.lambdaLower == .2);
report = struct('matlabVersion',version,'source',source, ...
    'scope','Saved samples only; no reintegration or continuous-time conclusion.', ...
    'edgeRule','r^2 < dc^2, diagonal excluded; no edge weight cutoff', ...
    'spectralWarningRule','100*N*eps*max(1,norm(L,2)); not a disconnection rule', ...
    'distanceWarningRule','100*eps*max(1,dc); never changes edges', ...
    'integrationSettings',R.integrationSettings);
keys = {'proposed','untightened'};
stems = {'tightened','untightened'};
series = struct();
for j = 1:2
    key = keys{j}; S = R.(key); t = S.t; nt = numel(t);
    lambda2=zeros(nt,1); tree=lambda2; components=lambda2; edgeCount=lambda2;
    minWeight=lambda2; maxDistance=lambda2; spectralTol=lambda2;
    underflowEdges=lambda2; nearRadiusPairs=lambda2;
    allWeights=zeros(nt,10); allDistances=zeros(nt,10); allSupport=false(nt,10);
    pairs=nchoosek(1:P.N,2); upper=triu(true(P.N),1);
    distanceTol=100*eps*max(1,P.dc);
    for k=1:nt
        Y=zeros(P.N,2);
        for i=1:P.N, cols=12*(i-1)+(1:2); Y(i,:)=S.z(k,cols); end
        ds=(Y(:,1)-Y(:,1).').^2+(Y(:,2)-Y(:,2).').^2;
        support=ds<P.dc^2; support(1:P.N+1:end)=false;
        W=zeros(P.N); W(support)=exp(-ds(support)./(P.dc^2-ds(support)));
        L=diag(sum(W,2))-W;
        ev=sort(eig(L)); lambda2(k)=ev(2); % retain raw eigenvalues
        tree(k)=P.N*det(L(1:P.N-1,1:P.N-1))/P.treeScale-P.lambdaLower;
        components(k)=numel(unique(conncomp(graph(support))));
        edgeCount(k)=nnz(support & upper);
        positive=W(support & W>0); if isempty(positive),minWeight(k)=NaN;else,minWeight(k)=min(positive);end
        maxDistance(k)=sqrt(max(ds(:)));
        spectralTol(k)=100*P.N*eps*max(1,norm(L,2));
        underflowEdges(k)=nnz(support & W==0 & upper);
        nearRadiusPairs(k)=nnz(abs(sqrt(ds)-P.dc)<=distanceTol & upper);
        for e=1:10
            a=pairs(e,1); b=pairs(e,2);
            allWeights(k,e)=W(a,b); allDistances(k,e)=sqrt(ds(a,b)); allSupport(k,e)=support(a,b);
        end
    end
    [minLambda,imin]=min(lambda2); [minTree,itree]=min(tree); negative=tree<0;
    h=struct('minLambda2',minLambda,'minLambda2Time',t(imin), ...
        'minSpectralMargin',minLambda-P.lambdaLower,'below02Samples',nnz(lambda2<P.lambdaLower), ...
        'below02Intervals',sampleIntervals(t,lambda2<P.lambdaLower), ...
        'disconnectedSamples',nnz(components>1),'disconnectedIntervals',sampleIntervals(t,components>1), ...
        'componentCounts',unique(components).','edgeCounts',unique(edgeCount).', ...
        'minPositiveWeight',min(minWeight),'maxPairDistance',max(maxDistance), ...
        'minDistanceToRadius',min(abs(allDistances-P.dc),[],'all'), ...
        'underflowedPositiveEdges',sum(underflowEdges),'nearRadiusPairSamples',sum(nearRadiusPairs), ...
        'distanceWarningTolerance',distanceTol,'spectralWarningToleranceRange',[min(spectralTol),max(spectralTol)], ...
        'nearZeroLambdaSamples',nnz(abs(lambda2)<=spectralTol), ...
        'minTree',minTree,'minTreeTime',t(itree),'lambda2AtMinTree',lambda2(itree), ...
        'negativeTreeIntervals',sampleIntervals(t,negative),'lambda2WhenTreeNegative',[]);
    if any(negative),h.lambda2WhenTreeNegative=[min(lambda2(negative)),max(lambda2(negative))];end
    report.(key)=h;
    T=table(t,lambda2,tree,lambda2-P.lambdaLower,components,edgeCount,minWeight,maxDistance,spectralTol, ...
        'VariableNames',{'time_s','lambda2','tree_certificate','lambda2_minus_02','components','edge_count', ...
        'min_positive_weight','max_pair_distance','spectral_warning_tolerance'});
    writetable(T,fullfile(outDir,[stems{j} '_connectivity_matlab.csv']));
    series.(key)=struct('table',T,'pairs',pairs,'weights',allWeights,'distances',allDistances,'support',allSupport);
    fprintf('%s: min lambda2=%.17g at t=%.17g, min weight=%.17g, disconnected samples=%d\n', ...
        key,minLambda,t(imin),min(minWeight),nnz(components>1));
end
fid=fopen(fullfile(outDir,'connectivity_matlab_summary.json'),'w'); assert(fid>=0);
c=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
save(fullfile(outDir,'connectivity_matlab.mat'),'report','series','-v7');
fprintf('MATLAB_CONNECTIVITY_PASS\n');
end

function rows=sampleIntervals(t,flag)
% These are observed sample endpoints, not exact continuous crossing times.
a=find(flag & ~[false;flag(1:end-1)]); b=find(flag & ~[flag(2:end);false]);
rows=struct('firstSample',{},'lastSample',{},'sampleCount',{},'entryBracket',{},'exitBracket',{});
for k=1:numel(a)
    entry=[]; departure=[];
    if a(k)>1,entry=[t(a(k)-1),t(a(k))];end
    if b(k)<numel(t),departure=[t(b(k)),t(b(k)+1)];end
    rows(k)=struct('firstSample',t(a(k)),'lastSample',t(b(k)), ...
        'sampleCount',b(k)-a(k)+1,'entryBracket',entry,'exitBracket',departure);
end
end
