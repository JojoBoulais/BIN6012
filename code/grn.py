from grnndata import GRNAnnData, read_h5ad
from grnndata import utils
import numpy as np
import anndata

grn = np.random.rand(10,10)
adata = anndata.AnnData(X=np.random.rand(10,10))

grn = GRNAnnData(adata, grn=grn)

print(grn) #shows the number of elements
grn.varp['GRN'] or grn.grn #shows the GRN
# subgrn = grn.get(['gene1', 'gene2']) #only gets some elements from the GRN
subgrn = grn.get(['gene1', 'gene2']) #only gets some elements from the GRN
subgrn.targets #shows the target connections
subgrn.plot() # displays the network

subgrn.write_h5ad('grn.h5ad') #writes it
read_h5ad('grn.h5ad') #reads it