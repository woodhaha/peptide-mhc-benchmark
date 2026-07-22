#!/usr/bin/env python3
"""Fig2: Clustering — find representative conformations from 100ns trajectory.
Subsamples Ca every 10 residues for speed → KMeans 5 clusters → centroid PDBs."""
import numpy as np, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.cluster import KMeans
import mdtraj as md

DATA = Path("D:/Researching/Peptide epitope/md results/salt_bridge")
OUT  = Path("D:/Researching/Peptide epitope/md pic")

print("Loading trajectory for clustering...", flush=True)
top = md.load(str(DATA / "final.pdb"))
traj = md.load(str(DATA / "traj.dcd"), top=str(DATA / "final.pdb"))
n_frames = traj.n_frames
t_ns = np.arange(n_frames) * 0.1

# Subsample Ca for speed (every 5th Ca → ~90 Ca atoms)
ca = traj.topology.select('protein and name CA')
ca_sub = ca[::5]
traj_ca = traj.atom_slice(ca_sub)
xyz = traj_ca.xyz.reshape(n_frames, -1)
print(f"  {len(ca)} Ca → {len(ca_sub)} subsampled, {xyz.shape[1]} dims", flush=True)

# KMeans
n_clusters = 5
km = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
labels = km.fit_predict(xyz)
counts = np.bincount(labels)

# Centroid frames
centroid_idx = []
for i in range(n_clusters):
    mask = labels == i
    if mask.sum() < 3:
        centroid_idx.append(0); continue
    dists = np.linalg.norm(xyz[mask] - km.cluster_centers_[i], axis=1)
    centroid_idx.append(int(np.where(mask)[0][np.argmin(dists)]))

# Save centroids
for i, fi in enumerate(centroid_idx):
    traj[fi].save(str(DATA / f"centroid_{i}_f{fi}.pdb"))
print(f"  Centroids: {[f'{fi} ({fi*0.1:.1f}ns)' for fi in centroid_idx]}", flush=True)
print(f"  Cluster sizes: {dict(zip(range(n_clusters), counts))}", flush=True)

np.savez(str(DATA / "clusters.npz"), labels=labels, centroid_idx=np.array(centroid_idx), t_ns=t_ns)

# NOTE: try PCA projection for visualization
ca_full = traj.atom_slice(ca)
xyz_full = ca_full.xyz.reshape(n_frames, -1)
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pc = pca.fit_transform(xyz_full / 10)  # scale to nm for numerical stability
var_exp = pca.explained_variance_ratio_ * 100

np.savez(str(DATA / "pca_proj.npz"), pc1=pc[:,0], pc2=pc[:,1], var_exp=var_exp)

# ── Figure: PCA projection colored by cluster + centroid markers ──
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle('K-E63 Trajectory Clustering — PCA Projection with Representative Conformations',
             fontsize=12, fontweight='bold')

colors = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6']
cluster_colors = [colors[l] for l in labels]

# 2a. PCA scatter colored by cluster
ax = axes[0]
sc = ax.scatter(pc[:,0], pc[:,1], c=[colors[l] for l in labels], s=8, alpha=0.6, edgecolors='none')
for i, fi in enumerate(centroid_idx):
    ax.scatter(pc[fi,0], pc[fi,1], c=colors[i], s=80, marker='D', edgecolors='black', linewidth=1,
               zorder=5, label=f'C{i+1} (t={fi*0.1:.0f}ns, n={counts[i]})')
# Start/end markers
ax.scatter(pc[0,0], pc[0,1], c='green', s=60, marker='*', edgecolors='black', zorder=5, label='Start')
ax.scatter(pc[-1,0], pc[-1,1], c='red', s=60, marker='v', edgecolors='black', zorder=5, label='End')
ax.set_xlabel(f'PC1 ({var_exp[0]:.1f}%)')
ax.set_ylabel(f'PC2 ({var_exp[1]:.1f}%)')
ax.set_title('Cluster Assignment (KMeans=5)')
ax.legend(fontsize=6, loc='upper right')
ax.grid(alpha=0.15)
# Population bar
for i in range(n_clusters):
    ax.text(pc[centroid_idx[i],0], pc[centroid_idx[i],1]-2,
            f'n={counts[i]}', fontsize=7, ha='center', color=colors[i], fontweight='bold')

# 2b. Time series of cluster membership
ax = axes[1]
for i in range(n_clusters):
    mask = labels == i
    ax.fill_between(t_ns, i-0.4, i+0.4, where=mask, color=colors[i], alpha=0.5, step='mid')
    ax.text(t_ns[len(t_ns)//2], i, f'C{i+1} ({counts[i]}fr)', fontsize=8, ha='center', va='center')
ax.set_xlim(0, 100)
ax.set_ylim(-0.5, n_clusters-0.5)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Cluster')
ax.set_yticks(range(n_clusters))
ax.set_yticklabels([f'C{i+1}' for i in range(n_clusters)])
ax.set_title('Cluster Occupancy Over Time')
ax.grid(axis='x', alpha=0.15)

# Print per-cluster time ranges
for i in range(n_clusters):
    cl_t = np.where(labels == i)[0] * 0.1
    if len(cl_t) > 0:
        print(f"  Cluster {i+1}: t={cl_t[0]:.0f}–{cl_t[-1]:.0f}ns, {len(cl_t)} frames ({counts[i]/n_frames*100:.1f}%)")

plt.tight_layout()
plt.savefig(OUT / 'Fig2_clustering.png', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'Fig2_clustering.pdf', bbox_inches='tight')
plt.close()
print('Fig2_clustering.png/pdf done')
