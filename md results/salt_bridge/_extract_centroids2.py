#!/usr/bin/env python3
"""Extract representative conformations (cluster centroids from PCA) for all compounds."""
import numpy as np, os, json, warnings
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

COMPS = ["NAFTAZONE","BIOTIN","ESKETAMINE","FUROSEMIDE",
         "GABAPENTIN_ENACARBIL","HYDROMORPHONE","PHENOBARBITAL","RIBOFLAVIN"]
BASE = "/root/autodl-tmp/SMVT_AutoDL_Package/trajectories"
SEP = "=" * 55

summary = []
for c in COMPS:
    d = BASE + "/" + c
    print("\n=== " + c + " ===", flush=True)
    u = mda.Universe(d + "/" + c + "_final.pdb", d + "/" + c + "_100ns.dcd")
    n = len(u.trajectory)
    ca = u.select_atoms("name CA")

    pca_d = np.load(d + "/" + c + "_pca_full.npz")
    X = np.column_stack([pca_d["pc1"], pca_d["pc2"]])

    best_k, best_s = 3, 0
    for k in range(2, 9):
        km = KMeans(n_clusters=k, n_init=10, random_state=42)
        lbl = km.fit_predict(X)
        s = silhouette_score(X, lbl)
        if s > best_s:
            best_s, best_k = s, k

    km = KMeans(n_clusters=best_k, n_init=10, random_state=42)
    labels = km.fit_predict(X)
    sizes = np.bincount(labels)

    cluster_info = []
    for cl in range(best_k):
        mask = labels == cl
        idx = np.where(mask)[0]
        pts = X[mask]
        centroid = km.cluster_centers_[cl]
        dists = np.sum((pts - centroid) ** 2, axis=1)
        nearest_idx = idx[np.argmin(dists)]
        cluster_info.append({
            "cluster": int(cl), "size": int(sizes[cl]),
            "pct": float(sizes[cl] / n * 100),
            "nearest_frame": int(nearest_idx),
            "nearest_time_ns": round(nearest_idx * 100.0 / 1000, 1),
        })

    np.savez_compressed(d + "/" + c + "_cluster_labels.npz",
                       labels=labels, n_clusters=best_k,
                       silhouette=best_s, sizes=sizes,
                       centroids=km.cluster_centers_)

    top = max(cluster_info, key=lambda x: x["pct"])
    print("  Best cluster #{}: {:.0f}% of frames".format(top["cluster"], top["pct"]), flush=True)
    print("  Centroid frame #{} @ {}ns".format(top["nearest_frame"], top["nearest_time_ns"]), flush=True)

    u.trajectory[top["nearest_frame"]]
    centroid_pdb = d + "/" + c + "_centroid.pdb"
    sel = u.select_atoms("protein or resname UNK")
    sel.write(centroid_pdb)
    print("  Saved {} ({} atoms)".format(centroid_pdb, sel.n_atoms), flush=True)

    summary.append({
        "compound": c, "n_clusters": best_k, "silhouette": round(best_s, 3),
        "top_cluster_pct": round(top["pct"], 1),
        "centroid_frame": top["nearest_frame"],
        "centroid_time_ns": top["nearest_time_ns"],
        "centroid_pdb": centroid_pdb,
        "cluster_sizes": [int(s) for s in sizes],
    })

json.dump(summary, open(BASE + "/centroid_summary.json", "w"), indent=2)
print("\n" + SEP)
for s in summary:
    print("  {:25s}  cluster#{}  top={:.0f}%  frame={} @ {}ns".format(
        s["compound"], s["n_clusters"], s["top_cluster_pct"],
        s["centroid_frame"], s["centroid_time_ns"]))
print(SEP)
print("ALL DONE", flush=True)
