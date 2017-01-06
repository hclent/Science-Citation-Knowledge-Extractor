import sys, os.path
try:
    from urllib.request import urlopen
except ImportError:
    from urllib import urlopen

url = "http://www.gutenberg.org/files/21325/21325-0.txt"
response = urlopen(url)
raw = response.read().decode('utf8')
print(type(raw))
print(len(raw))
print(raw[:75])

#print to txt
prefix = '/home/hclent/data/corpora'
name = 'grecoroman_med'
completeName = os.path.join(prefix, (str(name)+'.txt'))  #pmcid.txt #save to suffix path
sys.stdout = open(completeName, "w")
print(raw)