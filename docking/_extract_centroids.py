"""Extract representative conformations (cluster centroids) for all compounds."""
import numpy as np, os, json
import MDAnalysis as mda, warnings; warnings.filterwarnings("ignore")
from sklearn.cluster import KMeans

COMPS = ["NAFTAZONE","BIOTIN","ESKETAMINE","FUROSEMIDE",
         "GABAPENTIN_ENACARBIL","HYDROMORPHONE","PHENOBARBITAL","RIBOFLAVIN"]
BASE = "/root/autodl-tmp/SMVT_AutoDL_Package/trajectories"

summary = []
for c in COMPS:
    d = f"{BASE}/{c}"
    print(f"\n=== {c} ===", flush=True)
    u = mda.Universe(f"{d}/{c}_final.pdb", f"{d}/{c}_100ns.dcd")
    n = len(u.trajectory)
    ca = u.select_atoms("name CA")
    n_ca = ca.n_atoms

    # Load PCA projections
    pca_d = np.load(f"{d}/{c}_pca_full.npz")
    X = np.column_stack([pca_d["pc1"], pca_d["pc2"]])

    # Find optimal k (same method as before)
    from sklearn.metrics import silhouette_score
    best_k, best_s = 3, 0
    for k in range(2, 9):
        km = KMeans(n_clusters=k, n_init=10, random_state=42)
        lbl = km.fit_predict(X)
        s = silhouette_score(X, lbl)
        if s > best_s: best_s, best_k = s, k

    km = KMeans(n_clusters=best_k, n_init=10, random_state=42)
    labels = km.fit_predict(X)
    sizes = np.bincount(labels)

    # For each cluster: find the frame closest to centroid
    cluster_info = []
    for cl in range(best_k):
        mask = labels == cl
        idx = np.where(mask)[0]
        clust_pts = X[mask]
        centroid = km.cluster_centers_[cl]
        # Frame closest to centroid in PC1-PC2 space
        dists = np.sum((clust_pts - centroid)**2, axis=1)
        nearest_idx = idx[np.argmin(dists)]

        # Also find lowest RMSD frame in this cluster
        cluster_info.append({
            "cluster": int(cl),
            "size": int(sizes[cl]),
            "pct": float(sizes[cl]/n*100),
            "nearest_frame": int(nearest_idx),
            "nearest_time_ns": round(nearest_idx * 100.0 / 1000, 1),
        })

    # Save per-frame cluster labels + centroid info
    np.savez_compressed(f"{d}/{c}_cluster_labels.npz",
                       labels=labels, n_clusters=best_k,
                       silhouette=best_s, sizes=sizes,
                       centroids=km.cluster_centers_)

    # The "best" conformation = biggest cluster centroid frame
    top = max(cluster_info, key=lambda x: x["pct"])
    print(f"  Best cluster #{top[cluster]}: {top[pct]:.0f}% of frames", flush=True)
    print(f"  Centroid frame #{top[nearest_frame]} @ {top[nearest_time_ns]}ns", flush=True)

    # Save centroid frame as PDB
    u.trajectory[top["nearest_frame"]]
    centroid_pdb = f"{d}/{c}_centroid.pdb"
    # Select only protein + ligand (no waters)
    sel = u.select_atoms("protein or resname UNK")
    sel.write(centroid_pdb)
    print(f"  Saved {centroid_pdb} ({sel.n_atoms} atoms)", flush=True)

    summary.append({
        "compound": c, "n_clusters": best_k, "silhouette": round(best_s, 3),
        "top_cluster_pct": round(top["pct"], 1),
        "centroid_frame": top["nearest_frame"],
        "centroid_time_ns": top["nearest_time_ns"],
        "centroid_pdb": centroid_pdb,
        "cluster_sizes": [int(s) for s in sizes],
    })

json.dump(summary, open(f"{BASE}/centroid_summary.json", "w"), indent=2)
print(f"\n{=*55}")
for s in summary:
    print(f"  {s[compound]:25s}  cluster#{s[n_clusters]}  "
          f"top={s[top_cluster_pct]:.0f}%  "
          f"frame={s[centroid_frame]} @ {s[centroid_time_ns]}ns")
print('='*55)
print("ALL DONE", flush=True)
