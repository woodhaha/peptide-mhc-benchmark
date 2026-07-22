"""Extract PCA+eigenvectors + DCCM for all compounds. Tuned for 226K atoms."""
import numpy as np, os
import MDAnalysis as mda, warnings
warnings.filterwarnings("ignore")

COMPS = ["NAFTAZONE","BIOTIN","ESKETAMINE","FUROSEMIDE",
         "GABAPENTIN_ENACARBIL","HYDROMORPHONE","PHENOBARBITAL","RIBOFLAVIN"]
BASE = "/root/autodl-tmp/SMVT_AutoDL_Package/trajectories"

for c in COMPS:
    d = f"{BASE}/{c}"
    print(f"\n=== {c} ===", flush=True)
    u = mda.Universe(f"{d}/{c}_final.pdb", f"{d}/{c}_100ns.dcd")
    n = len(u.trajectory)
    ca = u.select_atoms("name CA")
    nc = ca.n_atoms

    # ── PCA (save eigenvectors for porcupine) ──
    out = f"{d}/{c}_pca_full.npz"
    if not os.path.exists(out):
        XYZ = np.empty((n, nc*3))
        for i, ts in enumerate(u.trajectory):
            XYZ[i] = ca.positions.flatten()
            if (i+1)%300==0: print(f"  PCA {i+1}/{n}", flush=True)
        m = XYZ.mean(0); Xc = XYZ - m
        U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        proj = U * S
        ve = (S**2)/(S**2).sum()
        np.savez_compressed(out, pc1=proj[:,0], pc2=proj[:,1],
                           proj=proj[:,:10], eigvec=Vt[:10,:],
                           mean=m, var_exp=ve[:20], eigvals=S[:20]**2)
        print(f"  PCA: PC1={ve[0]*100:.1f}% PC2={ve[1]*100:.1f}%", flush=True)

    # ── DCCM (sampled to ~200 frames) ──
    out2 = f"{d}/{c}_dccm.npy"
    if not os.path.exists(out2):
        step = max(1, n//200)
        # For 226K atoms, full Ca-Ca is ~600x600 = fine
        pos_list = []
        for i, ts in enumerate(u.trajectory):
            if i%step==0: pos_list.append(ca.positions.copy())
        XYZ = np.array(pos_list).reshape(len(pos_list), nc*3)
        corr = np.corrcoef(XYZ.T)
        np.save(out2, corr)
        print(f"  DCCM: {corr.shape}", flush=True)

    # ── Cosine content ──
    out3 = f"{d}/{c}_cosine.npz"
    if not os.path.exists(out3):
        Tn = np.arange(n)
        cos_vals = np.empty(5)
        for k in range(5):
            cosk = np.cos(k*np.pi*Tn/(n-1))
            num = (proj[:n,k]*cosk).sum()**2
            den = (proj[:n,k]**2).sum()
            cos_vals[k] = 2*num/((n-1)*den) if den>0 else 0
        np.savez_compressed(out3, cosine_content=cos_vals)
        print(f"  Cosine: {cos_vals}", flush=True)

    print(f"  DONE", flush=True)
print("\nALL DONE", flush=True)
