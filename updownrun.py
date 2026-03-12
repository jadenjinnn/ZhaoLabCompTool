from databases import DrugBank
from utils import get_genes_from_directory_updown
from l2s2_updown import enrich_l2s2_up_down

disease_genes_celltype = get_genes_from_directory_updown("./data/sources/Huntington/HL_puHD_E_CTRL")

db = DrugBank()

for name, genes in disease_genes_celltype.items(): 
    enrich_l2s2_up_down(list(genes[0].keys()),list(genes[1].keys()), db, name, "E_l2s2_updown")
