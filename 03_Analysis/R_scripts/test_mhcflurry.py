import warnings; warnings.filterwarnings('ignore')
from mhcflurry import Class1AffinityPredictor
p = Class1AffinityPredictor.load()
print('MHCflurry loaded OK')

df = p.predict_to_dataframe(
    peptides=['LLTDAQRIV','LMAFYLYEV','VMSPITLPT','SLHLTNCFV'],
    allele='HLA-A*02:01',
    include_percentile_ranks=True
)
print('Columns:', df.columns.tolist())
print(df.head().to_string())
