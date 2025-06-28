#lists
'''
gene_list = ['tp53', 'pca12', 'brca1']
gene_list.append('sum25')
print(gene_list)
gene_list.insert(2,'gene5')
print(gene_list)
#gene_list.remove('gene5')
#print(gene_list) 
#gene_list.pop(2)
#print(gene_list)
#To get index of the element:
print(gene_list.index('gene5'))
num_list = [1,45,65,34,23,35,66]
num_list.sort()
print(num_list)

#TUPLES:
gene_1 = ((110,200),(500,600))
for start , end in gene_1:
    print('Start position of gene1 is', start ,'End postion of gene1 is', end)

#Sets
bases = set('ATGCGCGCTGCU')
print('Unique bases are:', bases )

#Dctionaries:
genes = {
    'brca1':1,
    'tp53':2,
    'gene3':3
}
print(genes)
for name, num in genes.items():
    print('The name of gene is: ', name,'The number is: ',num)
'''
