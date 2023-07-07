import anndata2ri
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.conversion import localconverter
import pandas as pd

# Load R packages
readRDS = ro.r['readRDS']
pandas2ri.activate()

organ_list = ['Blood', 'BM', 'Brain', 'Liver', 'LN', 'SmallInt', 'Spleen']
disease_list = ['WT', 'AD']
ms_path = '/Users/kataokayuunosuke/MS'
image_path = f'{ms_path}/image'

organ = 'Blood'
organ_path = "/Users/kataokayuunosuke/MS/rds/Blood_all.integrated_doublet.removed_annotated.rds"

organ_rdata = readRDS(organ_path)
ro.r("organ_rdata <- readRDS('/Users/kataokayuunosuke/MS/rds/Blood_all.integrated_doublet.removed_annotated.rds')")
organ_xdata = ro.r("organ_rdata@assays[['RNA']]@counts")

