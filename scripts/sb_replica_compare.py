"""Replica comparison: distance overlay, occupancy, RMSF line plot, summary table."""
import json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9, 'figure.dpi': 150, 'font.family': 'sans-serif'})

OUT = r'D:\Researching\Peptide epitope\md pic'
MAIN = r'D:\Researching\Peptide epitope\md results\salt_bridge'
R2 = r'D:\Researching\Peptide epitope\md results\salt_bridge_replica2'
R3 = r'D:\Researching\Peptide epitope\md results\salt_bridge_replica3'

# ── Data ──
d_main = np.load(f'{MAIN}/md_dist.npz')
d_r2   = np.load(f'{R2}/md_dist_s42.npz')
d_r3   = np.load(f'{R3}/md_dist_s99.npz')

# Main uses 'mo' (min OE distance), replicas use 'minOE' — same metric
dist = {'Main (GBSA)': d_main['mo'], 'Rep2 (PME)': d_r2['minOE'], 'Rep3 (PME)': d_r3['minOE']}
sb   = {'Main (GBSA)': d_main['sb'], 'Rep2 (PME)': d_r2['sb'], 'Rep3 (PME)': d_r3['sb']}
t_ns = d_main['t_ns']
colors = {'Main (GBSA)': '#e74c3c', 'Rep2 (PME)': '#3498db', 'Rep3 (PME)': '#2ecc71'}
labels_full = {'Main (GBSA)': 'Main (GBSA)', 'Rep2 (PME)': 'Rep2 — seed=42 (PME)', 'Rep3 (PME)': 'Rep3 — seed=99 (PME)'}
labels_short = {'Main (GBSA)': 'Main', 'Rep2 (PME)': 'Rep2', 'Rep3 (PME)': 'Rep3'}

# RMSF
r_main = np.load(f'{MAIN}/rmsf.npz')
r2f, r2c = f'{R2}/rmsf.npz', f'{R3}/rmsf.npz'
try: r_r2 = np.load(r2f); r_r3 = np.load(r3f)
except: r_r2 = r_r3 = None

# ── Figure 1: 3-panel comparison ──
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
fig.suptitle('K-E63 Salt Bridge — Replica Comparison', fontsize=13, fontweight='bold', y=1.02)

# Panel A: histogram overlay
ax = axes[0]
for name, d in dist.items():
    ax.hist(d, bins=60, density=True, alpha=0.4, color=colors[name], label=labels_short[name])
    mu, std = d.mean(), d.std()
    xs = np.linspace(mu-4*std, mu+4*std, 200)
    kde = np.exp(-0.5*((xs-mu)/std)**2)/(std*np.sqrt(2*np.pi))
    ax.plot(xs, kde, color=colors[name], linewidth=1.5)
ax.axvline(4.0, color='gray', ls='--', lw=0.8)
ax.set_xlabel('CA-E63 min distance (Å)'); ax.set_ylabel('Density')
ax.set_title('A  Distance Distribution', loc='left', fontweight='bold')
ax.legend(fontsize=8)

# Panel B: distance time series (smoothed)
ax = axes[1]
window = 25
for name, d in dist.items():
    smooth = np.convolve(d, np.ones(window)/window, mode='valid')
    ax.plot(np.arange(len(smooth))/10, smooth, color=colors[name], lw=0.6, label=labels_short[name])
ax.axhline(4.0, color='gray', ls='--', lw=0.8)
ax.set_xlabel('Time (ns)'); ax.set_ylabel('CA-E63 min distance (Å)')
ax.set_title('B  Distance Time Series (25-frame MA)', loc='left', fontweight='bold')
ax.legend(fontsize=8)

# Panel C: occupancy bar
ax = axes[2]
names = list(sb.keys())
occ = [s.mean()*100 for s in sb.values()]
bars = ax.bar(range(len(names)), occ, color=[colors[n] for n in names], width=0.5)
ax.set_xticks(range(len(names))); ax.set_xticklabels([labels_short[n] for n in names])
ax.set_ylabel('Salt Bridge Occupancy (%)')
ax.set_ylim(0, 110)
ax.set_title('C  SB Occupancy (< 4.0 Å)', loc='left', fontweight='bold')
for bar, v in zip(bars, occ):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+2, f'{v:.1f}%', ha='center', fontsize=9)

plt.tight_layout()
plt.savefig(f'{OUT}/FigSX_replica_comparison.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/FigSX_replica_comparison.pdf', bbox_inches='tight')
plt.close()
print('1. FigSX_replica_comparison.png/pdf')

# ── Figure 2: RMSF 3-line overlay (if replica RMSF exists) ──
if r_r2 is not None and r_r3 is not None:
    fig, ax = plt.subplots(figsize=(12, 3.5))
    ax.plot(r_main['rmsf'], color=colors['Main (GBSA)'], lw=0.8, label='Main')
    ax.plot(r_r2['rmsf'], color=colors['Rep2 (PME)'], lw=0.8, alpha=0.7, label='Rep2')
    ax.plot(r_r3['rmsf'], color=colors['Rep3 (PME)'], lw=0.8, alpha=0.7, label='Rep3')
    ax.set_xlabel('Residue index'); ax.set_ylabel('RMSF (Å)')
    ax.set_title('Per-Residue RMSF — 3-Run Comparison')
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(f'{OUT}/FigSX_rmsf_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{OUT}/FigSX_rmsf_comparison.pdf', bbox_inches='tight')
    plt.close()
    print('2. FigSX_rmsf_comparison.png/pdf')
else:
    print('2. RMSF comparison skipped (replica RMSF not available)')

# ── Summary table ──
rows = []
for name in names:
    d = dist[name]
    rows.append(f'{name:20s} & {d.mean():.3f} & {d.std():.3f} & {d.min():.3f} & {d.max():.3f} & {sb[name].mean()*100:.1f}\\% \\\\')
summary = (
    'K-E63 Salt Bridge — 3-Run Comparison\n'
    + '='*70 + '\n'
    + f'{"Run":20s} {"Mean(Å)":>8s} {"SD(Å)":>6s} {"Min":>6s} {"Max":>6s} {"SB%":>6s}\n'
    + '-'*70 + '\n'
)
for name in names:
    d = dist[name]
    summary += f'{name:20s} {d.mean():8.3f} {d.std():6.3f} {d.min():6.3f} {d.max():6.3f} {sb[name].mean()*100:6.1f}\n'
summary += '='*70

with open(f'{MAIN}/replica_comparison.txt', 'w') as f:
    f.write(summary)
print(summary)
print('3. replica_comparison.txt saved')
