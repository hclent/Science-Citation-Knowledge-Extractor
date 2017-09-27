import os

'''
jointxts is for creating a massive text file containing all of our corpus.
with this we are able to train fasttext.
Once you have the jointxt file, you can cd into the fasttext repo and
./fasttext skipgram -input path/to/jointxt.txt -output model 
And this will produce your new model  
'''

pmcids_folder = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', 'data', 'pmcids'))


folders = os.listdir(pmcids_folder)
subfolders = [  os.path.abspath(os.path.join(pmcids_folder, f))  for f in folders ]


train_files = []

for s in subfolders:
    try:
        intrafolders = os.listdir(s)
        for i in intrafolders:
            look_folders = os.path.abspath(os.path.join(s, i))
            files = os.listdir(look_folders)
            for texts in files:
                if '.txt' in texts:
                    fullname = os.path.abspath(os.path.join(look_folders, texts))
                    train_files.append(fullname)
    except Exception as e:
        print(str(s) + " is not right ...")

# print("BIO TEXTS: ")
# print(len(train_files)) #17,108
#
#
# with open('/home/hclent/data/nns/cyverse/cyverseModel.txt', 'w') as outfile:
#     for fname in train_files:
#         with open(fname) as infile:
#             for line in infile:
#                 outfile.write(line)
#
# print("done! :)")