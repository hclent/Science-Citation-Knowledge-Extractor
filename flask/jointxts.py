import os


pmcids_folder = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', 'data', 'pmcids'))


folders = os.listdir(pmcids_folder)
subfolders = [  os.path.abspath(os.path.join(pmcids_folder, f))  for f in folders ]


train_files = []

for s in subfolders:
    intrafolders = os.listdir(s)
    for i in intrafolders:
        look_folders = os.path.abspath(os.path.join(s, i))
        files = os.listdir(look_folders)
        for texts in files:
            if '.txt' in texts:
                fullname = os.path.abspath(os.path.join(look_folders, texts))
                train_files.append(fullname)

print(len(train_files)) #16,721



with open('/home/hclent/data/nns/16k/16721.txt', 'w') as outfile:
    for fname in train_files:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

print("done! :)")