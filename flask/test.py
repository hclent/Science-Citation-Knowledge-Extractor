from content_management import *
from nes import *

nes_categories= ['BioProcess', 'CellLine', 'Cellular_component', 'Family', 'Gene_or_gene_product', 'Organ', 'Simple_chemical', 'Site', 'Species', 'TissueType']
w_number = 150
nes_list =  pickle.load(open("/home/hclent/data/nes/nes_18952863+18269575.pickle", "rb")) #pre-processed already
data_samples =  pickle.load(open("/home/hclent/data/data_samples/data_samples_18952863+18269575.pickle", "rb")) #pre-processed already
query = '18952863+18269575'
print("making the x, y, z")
x, y, z = vis_heatmap(data_samples, nes_list, nes_categories, w_number)
print("making the seadata")
seaData = make_seaborn_data(x, y, z)
print("saving the seamap")
blah = makeClusterMap(seaData, query)
print(blah)
print("done")
