"""Continue from saved model_comparison.csv — run PSSM, CV, IEDB."""
import os, time, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
rng = np.random.default_rng(42)

AA = list("ARNDCQEGHILKMFPSTWYV")
AA2IDX = {aa: i for i, aa in enumerate(AA)}
BLOSUM62 = np.array([[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],[-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],[-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],[-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],[0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],[-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],[-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],[0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],[-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],[-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],[-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],[-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2],[-2,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1],[-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1],[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2],[1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],[0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0],[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3],[-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1],[0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4]], dtype=np.float32)
row_min = BLOSUM62.min(1, keepdims=True); row_max = BLOSUM62.max(1, keepdims=True)
BLOSUM_NORM = (BLOSUM62 - row_min) / (row_max - row_min + 1e-8)

import keras; from keras import layers
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import StratifiedKFold

df = pd.read_csv("02_Data/raw/real_peptides.csv")
peps = df["peptide"].values; n = len(peps)
X = np.zeros((n, 180), dtype=np.float32)
for j in range(9):
    idx = [AA2IDX[p[j]] for p in peps]
    X[:, j*20:(j+1)*20] = BLOSUM_NORM[idx]
tm = (df.data_type == "train").values
Xtr, Xte = X[tm], X[~tm]
ytr = df.label_num.values[tm].astype(np.int32)
yte = df.label_num.values[~tm].astype(np.int32)

def bdffn():
    return keras.Sequential([layers.Input((180,)),layers.Dense(360,activation="relu"),layers.BatchNormalization(),layers.Dropout(0.5),layers.Dense(180,activation="relu"),layers.BatchNormalization(),layers.Dropout(0.4),layers.Dense(90,activation="relu"),layers.Dropout(0.3),layers.Dense(45,activation="relu"),layers.Dense(3,activation="softmax")],name="Deep_FFN")
def bffn():
    return keras.Sequential([layers.Input((180,)),layers.Dense(180,activation="relu"),layers.Dropout(0.4),layers.Dense(90,activation="relu"),layers.Dropout(0.3),layers.Dense(3,activation="softmax")],name="FFN_Jessen")
def bcnn():
    return keras.Sequential([layers.Input((9,20,1)),layers.Conv2D(32,(3,3),activation="relu",padding="same"),layers.Dropout(0.25),layers.Flatten(),layers.Dense(180,activation="relu"),layers.Dropout(0.4),layers.Dense(90,activation="relu"),layers.Dropout(0.3),layers.Dense(3,activation="softmax")],name="CNN")
def blstm():
    return keras.Sequential([layers.Input((9,20)),layers.LSTM(64),layers.Dropout(0.3),layers.Dense(32,activation="relu"),layers.Dense(3,activation="softmax")],name="LSTM")
def bres():
    inp = layers.Input((180,)); x = layers.Dense(128,activation="relu")(inp); skip = x; x = layers.Dense(128,activation="relu")(x); x = layers.BatchNormalization()(x); x = layers.add([x,skip]); x = layers.Dense(64,activation="relu")(x); skip2 = x; x = layers.Dense(64,activation="relu")(x); x = layers.BatchNormalization()(x); x = layers.add([x,skip2]); x = layers.Dropout(0.3)(x); x = layers.Dense(3,activation="softmax")(x); return keras.Model(inp,x,name="ResNet")

t0 = time.perf_counter()

# ===== PSSM =====
print("=== PSSM Labeling Comparison ===")
pssm = {f"p{i}": d for i, d in enumerate([
    dict(A=0.1,R=0.0,N=0.0,D=0.0,C=0.0,Q=0.0,E=0.0,G=0.2,H=0.0,I=0.2,L=0.1,K=0.0,M=0.1,F=0.4,P=0.0,S=0.2,T=0.2,W=0.3,Y=0.5,V=0.1),
    dict(A=0.5,R=-1,N=-1,D=-1,C=-1,Q=0.3,E=-1,G=0,H=-1,I=0.7,L=1,K=-1,M=0.9,F=0,P=-1,S=0,T=0.3,W=-1,Y=0,V=0.7),
    dict(A=0,R=0,N=0,D=0.3,E=0.3,Q=0,G=0,H=0,I=0,L=0,K=0,M=0,F=0,P=0,S=0,T=0,C=-0.2,W=0,Y=0,V=0),
    dict(A=0.1,R=0.1,N=0,D=0,E=0,Q=0.1,G=0,H=0.1,I=0.1,L=0.1,K=0.1,M=0,F=0,P=-0.2,S=0.1,T=0.1,C=-0.1,W=0,Y=0,V=0.1),
    dict(A=0,R=0.1,N=0,D=0,E=0,Q=0,G=0,H=0.1,I=0,L=0.1,K=0,M=0,F=0,P=-0.1,S=0,T=0,C=-0.1,W=0,Y=0,V=0),
    dict(A=0,R=0,N=0,D=0,E=0.1,Q=0,G=0,H=0,I=0.1,L=0.1,K=0,M=0,F=0,P=-0.1,S=0,T=0,C=0,W=0,Y=0,V=0),
    dict(A=0,R=0,N=0,D=0,E=0,Q=0,G=0,H=0,I=0,L=0.1,K=0,M=0,F=0.1,P=-0.1,S=0,T=0,C=0,W=0,Y=0,V=0),
    dict(A=0,R=0,N=0,D=0,E=0,Q=0,G=0,H=0,I=0.1,L=0.2,K=0,M=0,F=0,P=-0.1,S=0,T=0,C=0,W=0,Y=0,V=0.1),
    dict(A=0.4,R=-1,N=-1,D=-1,C=-1,Q=0.2,E=-1,G=0,H=-1,I=0.7,L=0.8,K=-1,M=0.5,F=0,P=-1,S=0.1,T=0.3,W=-1,Y=0,V=1.0),
], 1)}

nc = 15000; gc = rng.integers(0, 20, (nc, 9))
gp = ["".join(AA[i] for i in r) for r in gc]
sc = np.zeros(nc)
for pos in range(9):
    pssm_d = pssm[f"p{pos+1}"]
    for i, p in enumerate(gp): sc[i] += pssm_d.get(p[pos], 0)
sbq = np.quantile(sc, 0.98); wbq = np.quantile(sc, 0.93)
pl = np.full(nc, "NB", dtype="<U2"); pl[sc >= wbq] = "WB"; pl[sc >= sbq] = "SB"
pln = np.zeros(nc, dtype=np.int32); pln[pl == "WB"] = 1; pln[pl == "SB"] = 2
mc = min((pln == i).sum() for i in range(3))
sel = np.concatenate([rng.choice(np.where(pln == i)[0], mc, replace=False) for i in range(3)])
sp = [gp[s] for s in sel]; sl = pln[sel]
Xp = np.zeros((len(sp), 180), dtype=np.float32)
for j in range(9):
    idx2 = [AA2IDX[p[j]] for p in sp]
    Xp[:, j*20:(j+1)*20] = BLOSUM_NORM[idx2]
nt = len(sel) // 10; ip = rng.permutation(len(sel))
Xpt, Xpe = Xp[ip[nt:]], Xp[ip[:nt]]; ypt, ype = sl[ip[nt:]], sl[ip[:nt]]
print(f"PSSM data: {len(sel):,} balanced ({nt} test)")

mr = pd.read_csv("02_Data/cleaned/model_comparison.csv", index_col=0)
cfg = [("Deep FFN", bdffn, "ffn"), ("FFN (Jessen)", bffn, "ffn"), ("CNN", bcnn, "cnn"),
       ("LSTM", blstm, "lstm"), ("ResNet", bres, "ffn"), ("Random Forest", None, "ffn")]
print(f"{'Model':<20} {'MHCflurry':>10} {'PSSM':>10} {'Gap':>8}")
print("-" * 50)
psr = {}
for name, bldr, mode in cfg:
    if name == "Random Forest":
        rf2 = RandomForestClassifier(100, random_state=42, n_jobs=-1)
        rf2.fit(Xpt, ypt); ypp = rf2.predict(Xpe)
    else:
        try: keras.backend.clear_session()
        except: pass
        if mode == "cnn": m = bldr(); Xt2, Xe2 = Xpt.reshape(-1,9,20,1), Xpe.reshape(-1,9,20,1)
        elif mode == "lstm": m = bldr(); Xt2, Xe2 = Xpt.reshape(-1,9,20), Xpe.reshape(-1,9,20)
        else: m = bldr(); Xt2, Xe2 = Xpt, Xpe
        m.compile(loss="sparse_categorical_crossentropy", optimizer=keras.optimizers.Adam(0.001), metrics=["accuracy"])
        m.fit(Xt2, ypt, epochs=150, batch_size=50, validation_split=0.2,
              callbacks=[keras.callbacks.EarlyStopping(monitor="val_loss", patience=10, restore_best_weights=True)], verbose=0)
        ypp = np.argmax(m.predict(Xe2, verbose=0), axis=1)
    pacc = accuracy_score(ype, ypp) * 100
    macc = float(mr.loc[name, "accuracy"]) * 100 if name in mr.index else 0
    f1p = f1_score(ype, ypp, average="macro")
    print(f"  {name:<20} {macc:>9.1f}% {pacc:>9.1f}% {pacc-macc:>+7.1f}pp")
    psr[name] = {"accuracy": pacc/100, "macro_f1": f1p}
pd.DataFrame(psr).T.round(4).to_csv("02_Data/cleaned/pssm_comparison.csv")
print("  Saved: pssm_comparison.csv")

# ===== 5-Fold CV =====
print("\n=== 5-Fold CV (Deep FFN) ===")
skf = StratifiedKFold(5, shuffle=True, random_state=42)
cva = []
for f, (ti, vi) in enumerate(skf.split(Xtr, ytr)):
    try: keras.backend.clear_session()
    except: pass
    m = bdffn()
    m.compile(loss="sparse_categorical_crossentropy", optimizer=keras.optimizers.Adam(0.001), metrics=["accuracy"])
    m.fit(Xtr[ti], ytr[ti], epochs=150, batch_size=50, validation_data=(Xtr[vi], ytr[vi]),
          callbacks=[keras.callbacks.EarlyStopping(monitor="val_loss", patience=10, restore_best_weights=True)], verbose=0)
    _, acc = m.evaluate(Xte, yte, verbose=0); cva.append(acc)
    print(f"  Fold {f+1}: {acc*100:.1f}%")
print(f"  CV: {np.mean(cva)*100:.1f}% +/- {np.std(cva)*100:.1f}%")
pd.DataFrame({"fold": range(1,6), "accuracy": [a*100 for a in cva]}).to_csv("02_Data/cleaned/cv_summary.csv", index=False)
print("  Saved: cv_summary.csv")

# ===== IEDB =====
print("\n=== IEDB Benchmark ===")
iedb = pd.read_csv("02_Data/cleaned/iedb_benchmark_results.csv")
iep = iedb["peptide"].values
Xi = np.zeros((len(iep), 180), dtype=np.float32)
for j in range(9):
    iidx = [AA2IDX.get(p[j], 0) if len(p) >= 9 else 0 for p in iep]
    Xi[:, j*20:(j+1)*20] = BLOSUM_NORM[iidx]
try: keras.backend.clear_session()
except: pass
best = bdffn()
best.compile(loss="sparse_categorical_crossentropy", optimizer=keras.optimizers.Adam(0.001), metrics=["accuracy"])
best.fit(Xtr, ytr, epochs=150, batch_size=50, validation_split=0.2,
         callbacks=[keras.callbacks.EarlyStopping(monitor="val_loss", patience=10, restore_best_weights=True)], verbose=0)
ipr = best.predict(Xi, verbose=0); ipc = np.argmax(ipr, axis=1)
tb = (iedb["true_label"] == "POS").values.astype(int)
pb = (ipc != 0).astype(int)
sens = (tb & pb).sum() / max(tb.sum(), 1)
spec = ((1-tb) & (1-pb)).sum() / max((1-tb).sum(), 1)
print(f"  Sensitivity: {sens*100:.1f}%  Specificity: {spec*100:.1f}%")
print(f"  TP={(tb&pb).sum()} FN={(tb&(~pb)).sum()} TN={((1-tb)&(1-pb)).sum()} FP={((1-tb)&pb).sum()}")

io = iedb.copy()
io["pred_SB"] = ipr[:,2]; io["pred_WB"] = ipr[:,1]; io["pred_NB"] = ipr[:,0]
io["pred_class"] = [["NB","WB","SB"][p] for p in ipc]; io["pred_binary"] = pb
io["binding_score"] = ipr[:,2] + ipr[:,1] * 0.5
io.to_csv("02_Data/cleaned/iedb_benchmark_results_v2.csv", index=False)
print("  Saved: iedb_benchmark_results_v2.csv")

# ===== Protein Scan =====
print("\n=== Protein Epitope Scanning ===")
proteins = {
    "MART-1": "MPREDAHFIYGYPKKGHGHSYTTAEEAAGIGILTVILGVLLLIGCWYCRRRNGYRALMDKSLHVGTQCALTRRCPQEGFDHRDSKVSLQEKNCEPVVPNAPPAYEKLSAEQSPPPYSP",
    "gp100/PMEL": "MDLVLKRCLLHLAVIGALLAVGATKVPRNQDWLGVSRQLRTKAWNRQLYPEWTEAQRLDCWRGGQVSLKVSNDGPTLIGANASFSIALNFPGSQKVLPDGQVIWVNNTIINGSQVWGGQPVYPQETDDACIFPDGGPCPSGSWSQKRSFVYVWKTWGQYWQVLGGPVSGLSIGTGRAMLGTHTMEVTVYHRRGSRSYVPLAHSSSAFTITDQVPFSVSVSQLRALDGGNKHFLRNQPLTFALQLHDPSGYLAEADLSYTWDFGDSSGTLISRALVVTHTYLEPGPVTAQVVLQAAIPLTSCGSSPVPGTTDGHRPTAEAPNTTAGQVPTTEVVGTTPGQAPTAEPSGTTSVQVPTTEVISTAPVQMPTAESTGMTPEKVPVSEVMGTTLAEMSTPEATGMTPAEVSIVVLSGTTAAQVTTTEWVETTARELPIPEPEGPDASSIMSTESITGSLGPLLDGTATLRLVKRQVPLDCVLYRYGSFSVTLDIVQGIESAEILQAVPSGEGDAFELTVSCQGGLPKEACMEISSPGCQPPAQRLCQPVLPSPACQLVLHQILKGGSGTYCLNVSLADTNSLAVVSTQLIMPGQEAGLGQVPLIVGILLVLMAVVLASLIYRRRLMKQDFSVPQLPHSSSHWLRLPRIFCSCPIGENSPLLSGQQV",
    "Tyrosinase": "MLLAVLYCLLWSFQTSAGHFPRACVSSKNLMEKECCPPWSGDRSPCGQLSGRGSCQNILLSNAPLGPQFPFTGVDDRESWPSVFYNRTCQCSGNFMGFNCGNCKFGFWGPNCTERRLLVRRNIFDLSAPEKDKFFAYLTLAKHTISSDYVIPIGTYGQMKNGSTPMFNDINIYDLFVWMHYYVSMDALLGGSEIWRDIDFAHEAPAFLPWHRLFLLRWEQEIQKLTGDENFTIPYWDWRDAEKCDICTDEYMGGQHPTNPNLLSPASFFSSWQIVCSRLEEYNSHQSLCNGTPEGPLRRNPGNHDKSRTPRLPSSADVEFCLSLTQYESGSMDKAANFSFRNTLEGFASPLTGIADASQSSMHNALHIYMNGTMSQVQGSANDPIFLLHHAFVDSIFEQWLRRHRPLQEVYPEANAPIGHNRESYMVPFIPLYRNGDFFISSKDLGYDYSYLQDSDPDSFQDYIKSYLEQASRIWSWLLGAAMVGAVLTALLAGLVSLLCRHKRKQLPEEKQPLLMEKEDYHSLYQSHL",
    "NY-ESO-1": "MQAEGRGTGGSTGDADGPGGPGIPDGPGGNAGGPGEAGATGGRGPRGAGAARASGPGGGAPRGPHGGAASGLNGCCRCGARGPESRLLEFYLAMPFATPMEAELARRSLAQDAPPLPVPGVLLKEFTVSGNILTIRLTAADHRQLQLSISSCLQQLSLLMWITQCFLPVFLAQPPSGQRR",
    "WT1": "MGSDVRDLNALLPAVPSLGGGGGCALPVSGAAQWAPVLDFAPPGASAYGSLGGPAPPPAPPPPPPPPPHSFIKQEPSWGGAEPHEEQCLSAFTVHFSGQFTGTAGACRYGPFGPPPPSQASSGQARMFPNAPYLPSCLESQPAIRNQGYSTVTFDGTPSYGHTPSHHAAQFPNHSFKHEDPMGQQGSLGEQQYSVPPPVYGCHTPTDSCTGSQALLLRTPYSSDNLYQMTSQLECMTWNQMNLGATLKGVAAGSSSSVKWTEGQSNHSTGYESDNHTTPILCGAQYRIHTHGVFRGIQDVRRVPGVAPTLVRSASETSEKRPFMCAYPGCNKRYFKLSHLQMHSRKHTGEKPYQCDFKDCERRFSRSDQLKRHQRRHTGVKPFQCKTCQRKFSRSDHLKTHTRTHTGKTSEKPFSCRWPSCQKKFARSDELVRHHNMHQRNMTKLQLAL",
    "p53": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD",
    "KRAS": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEK",
    "M1": "MASQGTKRSYEQMETDGERQNATEIRASVGKMIGGIGRFYIQMCTELKLSDYEGRLIQNSLTIENMVLWVLSKEERRYFPRTKGLRIAYNEVGHTNEDEYTQEDVKKSRFEIQKKLREMATKDDVKKFQSKLEKEELDLVKKMLKDAIGQDKKL",
    "CMV pp65": "MESRGRRCPEMISVLGPISGHVLKAVFSRGDTPVLPHETRLLQTGIHVRVSQPSLILVSQYTPDSTPCHRGDNQLQVQHTYFTGSEVENVSVNVHNPTGRSICPSQEPMSIYVYALPLKMLNIPSINVHHYPSAAERKHRHLPVADAVIHASGKQMWQARLTVSGLAWTRQQNQWKEPDVYFTSAFVFPTKDVALRHVVCAHELVCSMENTRATKMQVIGDQYVKVYLESFCEDVPSGKLFMHVTLGSDVEEDLTMTRNPQPFMRPHERNGFTVLCPKNMIIKPGKISHIMLDVAFTSHEHFGLLCPKSIPGLSISGNLLMNGQQIFLEVQAIRETVELRQYDPVAALFFFDIDLLLQRGPQYSEHPTFTSQYRIQGKLEYRHTWDRHDEGAAQGDDDVWTSGSDSDEELVTTERKTPRVTGGGAMAGASTSAGRKRKSASSATACTAGVMTRGRLKAESTVAPEEDTDEDSDNEIHNPAVFTWPPWQAGILARNLVPMVATVQGQNLKYQEFFWDANDIYRIFAELEGVWQPAAQPKRRRHRQDALPGPCIASTPKKHRG",
    "Spike RBD": "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF",
}
print(f"  Scanning {len(proteins)} proteins...")
scan_rows = []
for name, seq in proteins.items():
    for start in range(len(seq) - 8):
        pep = seq[start:start+9]
        if len(pep) != 9: continue
        p2 = pep[1] if len(pep) > 1 else ""
        p9 = pep[8] if len(pep) > 8 else ""
        scan_rows.append({"protein": name, "start_pos": start+1, "end_pos": start+9,
                          "peptide": pep, "p2": p2, "p9": p9})

scan_df = pd.DataFrame(scan_rows)
scan_peps = scan_df["peptide"].values
Xs = np.zeros((len(scan_peps), 180), dtype=np.float32)
for j in range(9):
    sidx = [AA2IDX[p[j]] for p in scan_peps]
    Xs[:, j*20:(j+1)*20] = BLOSUM_NORM[sidx]
sprobs = best.predict(Xs, verbose=0); spreds = np.argmax(sprobs, axis=1)
scan_df["prob_SB"] = sprobs[:,2]; scan_df["prob_WB"] = sprobs[:,1]; scan_df["prob_NB"] = sprobs[:,0]
scan_df["pred_class"] = [["NB","WB","SB"][p] for p in spreds]
scan_df["binding_score"] = sprobs[:,2] + sprobs[:,1]*0.5
scan_df["anchors"] = scan_df["p2"].isin(["L","M","I","V"]) & scan_df["p9"].isin(["V","L","I"])
scan_df.to_csv("02_Data/cleaned/protein_epitope_scan_extended.csv", index=False)
print(f"  Scanned {len(scan_df):,} 9-mers across {len(proteins)} proteins")
print(f"  SB: {(spreds==2).sum()}  WB: {(spreds==1).sum()}  NB: {(spreds==0).sum()}")
print("  Saved: protein_epitope_scan_extended.csv")

# ===== Mutation Scan =====
print("\n=== Mutation Scanning ===")
mutations = [
    ("p53", 175, "R", "H"), ("p53", 220, "Y", "C"), ("p53", 245, "G", "S"),
    ("p53", 248, "R", "W"), ("p53", 249, "R", "S"), ("p53", 273, "R", "H"),
    ("p53", 282, "R", "W"),
    ("KRAS", 12, "G", "D"), ("KRAS", 12, "G", "V"), ("KRAS", 12, "G", "C"),
    ("KRAS", 12, "G", "R"), ("KRAS", 13, "G", "D"), ("KRAS", 61, "Q", "H"),
    ("KRAS", 61, "Q", "L"), ("KRAS", 61, "Q", "R"), ("KRAS", 146, "A", "T"),
]
mut_rows = []
for prot, pos, wt, mut in mutations:
    seq = proteins[prot]
    pos0 = pos - 1
    for w in range(-4, 5):
        start = pos0 - 4 + w
        end = start + 9
        if start < 0 or end > len(seq): continue
        wt_pep = seq[start:end]
        mut_pep = list(wt_pep)
        mut_idx = pos0 - start
        mut_pep[mut_idx] = mut
        mut_pep = "".join(mut_pep)
        mut_rows.append({"protein": prot, "mutation": f"{wt}{pos}{mut} -- {prot}",
                         "position": pos, "wt_aa": wt, "mut_aa": mut,
                         "window_start": start+1, "wt_peptide": wt_pep, "mut_peptide": mut_pep})
mut_df = pd.DataFrame(mut_rows)
mut_df = mut_df.drop_duplicates(subset=["protein","mutation","mut_peptide"])

mw = mut_df["wt_peptide"].values; mm = mut_df["mut_peptide"].values
Xw = np.zeros((len(mw), 180), dtype=np.float32)
Xm = np.zeros((len(mm), 180), dtype=np.float32)
for j in range(9):
    Xw[:, j*20:(j+1)*20] = BLOSUM_NORM[[AA2IDX[p[j]] for p in mw]]
    Xm[:, j*20:(j+1)*20] = BLOSUM_NORM[[AA2IDX[p[j]] for p in mm]]
wp = best.predict(Xw, verbose=0); mp = best.predict(Xm, verbose=0)
wscore = wp[:,2] + wp[:,1]*0.5; mscore = mp[:,2] + mp[:,1]*0.5
delta = mscore - wscore
mut_df["wt_score"] = wscore; mut_df["mut_score"] = mscore
mut_df["delta_score"] = delta
mut_df["wt_class"] = [["NB","WB","SB"][p] for p in np.argmax(wp, axis=1)]
mut_df["mut_class"] = [["NB","WB","SB"][p] for p in np.argmax(mp, axis=1)]
mut_df["effect"] = "unchanged"
for i in range(len(mut_df)):
    if mut_df["wt_class"].iloc[i] in ["NB"] and mut_df["mut_class"].iloc[i] in ["SB","WB"]:
        mut_df.loc[mut_df.index[i], "effect"] = "CREATED (neoepitope)"
    elif mut_df["delta_score"].iloc[i] > 0.1:
        mut_df.loc[mut_df.index[i], "effect"] = "ENHANCED"
    elif mut_df["wt_class"].iloc[i] in ["SB","WB"] and mut_df["mut_class"].iloc[i] == "NB":
        mut_df.loc[mut_df.index[i], "effect"] = "DESTROYED"
    elif mut_df["delta_score"].iloc[i] < -0.05:
        mut_df.loc[mut_df.index[i], "effect"] = "WEAKENED"
mut_df.to_csv("02_Data/cleaned/mutation_scan_results.csv", index=False)
print(f"  Scanned {len(mut_df)} mutation windows")
eff = mut_df[mut_df["effect"] != "unchanged"]
print(f"  Epitope-altering: {len(eff)}")
for _, r in eff.iterrows():
    print(f"    {r['mutation']:30s} {r['effect']:25s} Δ={r['delta_score']:+.3f}")
print("  Saved: mutation_scan_results.csv")

# ===== Top 20 Epitopes =====
print("\n=== Top 20 Epitope Candidates ===")
top20 = scan_df.nlargest(20, "binding_score")[
    ["protein","start_pos","peptide","p2","p9","pred_class","binding_score"]
].reset_index(drop=True)
top20.to_csv("02_Data/cleaned/top20_epitope_candidates.csv", index=False)
for i, r in top20.head(10).iterrows():
    print(f"  {i+1:2d}. {r['peptide']}  {r['protein']}  {r['binding_score']:.3f}")
print("  Saved: top20_epitope_candidates.csv")

print(f"\n{'='*50}")
print(f"ALL DONE — {time.perf_counter()-t0:.0f}s")
