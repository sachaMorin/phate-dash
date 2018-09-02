import numpy as np
import phate
import os
import scprep
import magic
import pickle

data_name = "PPA"
data = scprep.io.load_10X_HDF5(os.path.join(
    os.path.abspath(os.sep),
    "data", "lab", "DataSets",
    "Krause_2018_primary_parathyroid_adenoma",
    "ParaY9_HHT_cellranger",
    "filtered_gene_bc_matrices_h5.h5",  # raw_gene_bc_matrices
), gene_labels='both', allow_duplicates=True)

data = scprep.filter.remove_rare_genes(data, min_cells=3)
data = scprep.normalize.library_size_normalize(data)
data = scprep.transform.sqrt(data)

ph = phate.PHATE(n_components=2)
phate_data = ph.fit_transform(data)
np.save("{}Phate2d.npy".format(data_name), phate_data)

ph.set_params(n_components=3)
phate3_data = ph.transform()
np.save("{}Phate3d.npy".format(data_name), phate3_data)

mg = magic.MAGIC()
_ = mg.fit_transform(data)

# reduce memory footprint
del mg.graph.data
del mg.graph.data_nu
del mg.graph._kernel
del mg.graph._diff_op
del mg.graph._knn_tree

with open('magic.pickle', 'wb') as handle:
    pickle.dump(mg, handle, protocol=pickle.HIGHEST_PROTOCOL)
