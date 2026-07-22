import numpy as np, MDAnalysis as mda, warnings, os, json
warnings.filterwarnings("ignore")

COMPS = ["NAFTAZONE","BIOTIN","ESKETAMINE","FUROSEMIDE",
         "GABAPENTIN_ENACARBIL","HYDROMORPHONE","PHENOBARBITAL","RIBOFLAVIN"]
BASE = "/root/autodl-tmp/SMVT_AutoDL_Package/trajectories"

for c in COMPS:
    d = BASE + "/" + c
    labels_f = d + "/" + c + "_cluster_labels.npz"
    if not os.path.exists(labels_f):
        print("SKIP " + c + ": no cluster_labels"); continue
    
    lbl = np.load(labels_f)
    labels = lbl["labels"]
    sizes = np.bincount(labels)
    top_cl = np.argmax(sizes)
    
    u = mda.Universe(d + "/" + c + "_final.pdb", d + "/" + c + "_100ns.dcd")
    n = len(u.trajectory)
    
    # Get PC1-PC2 to find centroid frame
    pca = np.load(d + "/" + c + "_pca_full.npz")
    X = np.column_stack([pca["pc1"], pca["pc2"]])
    mask_top = labels == top_cl
    idx_top = np.where(mask_top)[0]
    pts = X[mask_top]
    centroid = pts.mean(axis=0)
    dists = np.sum((pts - centroid)**2, axis=1)
    best_idx = idx_top[np.argmin(dists)]
    
    u.trajectory[best_idx]
    # Fix: select protein OR LIG (not resname UNK)
    sel = u.select_atoms("protein or resname LIG")
    out = d + "/" + c + "_centroid.pdb"
    sel.write(out)
    
    n_lig = len(u.select_atoms("resname LIG"))
    print(c + ": frame=" + str(best_idx) + " prot=" + str(sel.n_atoms) + " lig_atoms=" + str(n_lig) + " -> " + out, flush=True)

print("ALL FIXED", flush=True)
