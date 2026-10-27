"""Recompute physical-graph diagnostics from the archived MATLAB trajectory.

No integration, state interpolation, controller change, or edge-weight cutoff.
Run: python audit_connectivity.py --out connectivity_check
Dependencies: numpy scipy matplotlib mpmath.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
import numpy as np
import scipy
from scipy.io import loadmat
from scipy.linalg import helmert
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def intervals(t, flag):
    """Saved-sample intervals and brackets; no interpolated crossing claims."""
    flag = np.asarray(flag, dtype=bool)
    starts = np.flatnonzero(flag & ~np.r_[False, flag[:-1]])
    ends = np.flatnonzero(flag & ~np.r_[flag[1:], False])
    return [dict(first_sample=float(t[a]), last_sample=float(t[b]),
                 sample_count=int(b-a+1),
                 entry_bracket=None if a == 0 else [float(t[a-1]), float(t[a])],
                 exit_bracket=None if b+1 == len(t) else [float(t[b]), float(t[b+1])])
            for a, b in zip(starts, ends)]


def component_count(a):
    unseen = set(range(len(a)))
    count = 0
    while unseen:
        count += 1
        stack = [unseen.pop()]
        while stack:
            i = stack.pop()
            neighbors = set(np.flatnonzero(a[i])) & unseen
            unseen.difference_update(neighbors)
            stack.extend(neighbors)
    return count


def high_precision(y, dc):
    """80-digit recomputation from exact saved binary64 coordinates.

    This checks spectral arithmetic, not the error of the ODE trajectory.
    """
    with mp.workdps(80):
        n = len(y)
        yy = [[mp.mpf(float(v)) for v in row] for row in y]
        d = mp.mpf(float(dc))**2
        L = mp.zeros(n)
        for i in range(n):
            for j in range(i+1, n):
                s = sum((yy[i][k]-yy[j][k])**2 for k in range(y.shape[1]))
                w = mp.exp(-s/(d-s)) if s < d else mp.mpf(0)
                L[i,i] += w; L[j,j] += w
                L[i,j] = L[j,i] = -w
        ev = mp.eigsy(L, eigvals_only=True)
        return mp.nstr(ev[1], 70)


def analyze(r, key, saved_plot, out):
    s, P = r[key], r['P']
    t = np.asarray(s['t'])
    n = int(P['N']); dc = float(P['dc']); floor = float(P['lambdaLower'])
    assert n == 5 and dc == 3.2 and floor == 0.2
    # State order: each agent owns a length-10 decision estimate. Its own
    # two entries are its physical position: MATLAB columns 12*(i-1)+[1,2].
    y = np.stack([s['z'][:,12*i:12*i+2] for i in range(n)], axis=1)
    assert np.array_equal(y, saved_plot['Y'])
    assert np.array_equal(t, saved_plot['t'])
    delta = y[:,:,None,:]-y[:,None,:,:]
    d2 = np.sum(delta**2, axis=3)
    adjacency = d2 < dc**2  # exact support test, no distance tolerance
    adjacency[:,np.arange(n),np.arange(n)] = False
    w = np.zeros_like(d2)
    exponent = np.full_like(d2, -np.inf)
    exponent[adjacency] = -d2[adjacency]/(dc**2-d2[adjacency])
    w[adjacency] = np.exp(exponent[adjacency])
    L = -w.copy()
    L[:,np.arange(n),np.arange(n)] = w.sum(axis=2)
    ev, vec = np.linalg.eigh(L)  # symmetric solver; retain raw values
    lam = ev[:,1]
    q = helmert(n, full=False).T
    reduced = np.linalg.eigvalsh(q.T @ L @ q)[:,0]
    residual = np.max(np.linalg.norm(L@vec-vec*ev[:,None,:],axis=1),axis=1)
    eps = np.finfo(float).eps
    spectral_tol = 100*n*eps*np.maximum(1,np.max(np.abs(ev),axis=1))
    distance_tol = 100*eps*max(1,dc)  # warning band only, never changes edges
    pairs = list(zip(*np.triu_indices(n,1)))
    distances = np.stack([np.sqrt(d2[:,i,j]) for i,j in pairs],axis=1)
    weights = np.stack([w[:,i,j] for i,j in pairs],axis=1)
    support = np.stack([adjacency[:,i,j] for i,j in pairs],axis=1)
    components = np.array([component_count(a) for a in adjacency])
    represented_components = np.array([component_count(a) for a in w > 0])
    h = n*np.linalg.det(L[:,:n-1,:n-1])/float(P['treeScale'])-floor
    neg = h < 0
    imin = int(np.argmin(lam)); ih = int(np.argmin(h))
    special = {0, imin, ih, len(t)-1}
    if neg.any():
        ii = np.flatnonzero(neg)
        special.update([int(ii[0]), int(ii[-1]), int(ii[np.argmin(lam[neg])]),
                        int(ii[np.argmax(lam[neg])])])
    special.update(np.flatnonzero((np.abs(lam) <= spectral_tol) |
                                 (np.abs(lam-floor) <= spectral_tol)).tolist())
    hp=[]
    for i in sorted(special):
        value=high_precision(y[i],dc)
        hp.append(dict(time=float(t[i]), lambda2_80digit=value,
                       double_value=float(lam[i]), abs_difference=abs(float(value)-lam[i])))
    topology=[]
    change=np.flatnonzero(np.r_[True,np.any(support[1:]!=support[:-1],axis=1)])
    for a,b in zip(change,np.r_[change[1:]-1,len(t)-1]):
        topology.append(dict(first_sample=float(t[a]),last_sample=float(t[b]),
                             edges=[[int(i+1),int(j+1)] for k,(i,j) in enumerate(pairs) if support[a,k]],
                             component_count=int(components[a])))
    summary=dict(sample_count=len(t),time_range=[float(t[0]),float(t[-1])],
        min_lambda2=float(lam[imin]),min_lambda2_time=float(t[imin]),
        min_lambda2_minus_02=float(lam[imin]-floor),
        below_02_sample_count=int(np.sum(lam<floor)),
        below_02_intervals=intervals(t,lam<floor),
        disconnected_sample_count=int(np.sum(components>1)),
        disconnected_intervals=intervals(t,components>1),
        component_counts=np.unique(components).tolist(),
        edge_counts=np.unique(support.sum(axis=1)).tolist(),topology_intervals=topology,
        min_positive_weight=float(w[adjacency & (w>0)].min()),
        underflowed_positive_edge_count=int(np.sum(adjacency & (w==0))//2),
        support_vs_represented_components_mismatches=int(np.sum(components!=represented_components)),
        max_pair_distance=float(distances.max()),
        min_distance_to_radius=float(np.min(np.abs(distances-dc))),
        near_radius_pair_samples=int(np.sum(np.abs(distances-dc)<=distance_tol)),
        distance_warning_tolerance=distance_tol,
        spectral_warning_tolerance_range=[float(spectral_tol.min()),float(spectral_tol.max())],
        near_zero_lambda2_samples=int(np.sum(np.abs(lam)<=spectral_tol)),
        near_02_lambda2_samples=int(np.sum(np.abs(lam-floor)<=spectral_tol)),
        max_eigenpair_residual=float(residual.max()),
        max_reduced_laplacian_difference=float(np.max(np.abs(lam-reduced))),
        max_archived_lambda2_difference=float(np.max(np.abs(lam-saved_plot['lambda2']))),
        max_archived_tree_difference=float(np.max(np.abs(h-saved_plot['hTree']))),
        min_tree_certificate=float(h[ih]),min_tree_time=float(t[ih]),
        lambda2_at_min_tree=float(lam[ih]),
        negative_tree_intervals=intervals(t,neg),
        lambda2_when_tree_negative=None if not neg.any() else [float(lam[neg].min()),float(lam[neg].max())],
        high_precision_checks=hp,
        figure_window=[0,400],min_lambda2_in_figure=float(lam[t<=400].min()))
    assert summary['max_archived_lambda2_difference'] < 1e-12
    assert summary['max_archived_tree_difference'] < 1e-12
    stem='tightened' if key=='proposed' else 'untightened'
    data=np.column_stack((t,lam,h,lam-floor,components,support.sum(axis=1),
                          weights.min(axis=1),distances.max(axis=1),spectral_tol,residual,
                          weights,distances,support.astype(int)))
    headers=['time_s','lambda2','tree_certificate','lambda2_minus_02','components','edge_count',
             'min_pair_weight','max_pair_distance','spectral_warning_tolerance','eigenpair_residual']
    headers += [f'{kind}_{i+1}_{j+1}' for kind in ('weight','distance','edge') for i,j in pairs]
    np.savetxt(out/f'{stem}_connectivity.csv',data,delimiter=',',header=','.join(headers),comments='',fmt='%.17g')
    return summary, dict(t=t,lambda2=lam,hTree=h,Y=y)


def plot(data, out):
    # Existing manuscript encoding: tightened red solid; untightened blue dashed.
    plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','DejaVu Serif'],
        'mathtext.fontset':'stix','font.size':8,'axes.labelsize':9,'axes.linewidth':0.5,
        'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','path.simplify':False})
    fig,ax=plt.subplots(figsize=(3.5,1.72))
    fig.subplots_adjust(left=.14,right=.955,bottom=.24,top=.95)
    for key,color,style,label in [('proposed','red','-','Tightened'),('untightened','blue','--','Untightened')]:
        d=data[key];take=d['t']<=400
        ax.plot(d['t'][take],d['lambda2'][take],style,color=color,lw=.85,label=label)
    ax.axhline(.2,color='black',linestyle=':',lw=.7)
    ax.text(396,.26,'Reference: 0.2',ha='right',va='bottom',fontsize=8)
    ax.set(xlim=(0,400),ylim=(0,4.2),xticks=np.arange(0,401,100),yticks=np.arange(0,5),
           xlabel=r'Time $t$ (s)',ylabel=r'$\lambda_2(L_{\chi}(y))$')
    ax.grid(True,color='.88',linewidth=.35)
    ax.tick_params(direction='in',length=2,width=.45)
    ax.legend(loc='upper right',ncol=2,frameon=False,fontsize=7.5,
              borderaxespad=.2,columnspacing=.7,handlelength=1.4,handletextpad=.25)
    for ext in ('eps','pdf','svg','png'):
        fig.savefig(out/f'algebraic_connectivity.{ext}',dpi=300,facecolor='white')
    plt.close(fig)


def main():
    root=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,default=root/'results/experiment_data.mat')
    ap.add_argument('--plot-data',type=Path,default=root/'results/plot_data.mat')
    ap.add_argument('--out',type=Path,default=root/'connectivity_check')
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    r=loadmat(args.source,simplify_cells=True)['result']
    p=loadmat(args.plot_data,simplify_cells=True)
    summary=dict(source_sha256=sha256(args.source),plot_data_sha256=sha256(args.plot_data),
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        matplotlib=matplotlib.__version__,mpmath=mp.__version__,
        method='Physical positions; symmetric full Laplacian eigenvalues; distance support and independent BFS.',
        thresholds={'edge_support':'squared distance < dc^2; diagonal excluded; no weight cutoff',
                    'spectral_warning':'100*N*binary64_epsilon*max(1,||L||_2); diagnostic only',
                    'distance_warning':'100*binary64_epsilon*max(1,dc); diagnostic only',
                    'tree_negative':'hTree < 0, no tolerance',
                    'spectral_margin_violation':'lambda2 < 0.2, no tolerance',
                    'disconnected':'distance-support component_count > 1',
                    'archive_agreement_assertion':1e-12},
        integration_settings=r['integrationSettings'],
        scope='Saved samples only. No continuous-time claim or reintegration; high precision does not improve ODE accuracy.')
    data={}
    for key,pkey in [('proposed','T'),('untightened','U')]:
        summary[key],data[key]=analyze(r,key,p[pkey],args.out)
    plot(data,args.out)
    (args.out/'connectivity_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    for key in data:
        s=summary[key]
        print(key, 'min lambda2 =',format(s['min_lambda2'],'.17g'),
              'at',s['min_lambda2_time'],'components',s['component_counts'],
              'min weight',s['min_positive_weight'])
    print('CONNECTIVITY_AUDIT_PASS')


if __name__=='__main__':
    main()
