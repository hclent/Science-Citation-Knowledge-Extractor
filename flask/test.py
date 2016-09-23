from content_management import *
import pickle


main_info = [] #main_info are formatted citations from citations table in db
target_urls = []


citations_18269575, urls_18269575 = new_citations_from_db(18269575)


citations_18952863, urls_18952863 = new_citations_from_db(18952863)

for mi in citations_18269575:
    main_info.append(mi)

print("first doc to main")

for mi in citations_18952863:
    main_info.append(mi)

print("second doc to main")

for url in urls_18269575:
    target_urls.append(url)

print("first doc to url")

for url in urls_18952863:
    target_urls.append(url)


print("second doc to url")


citations_with_links = list(zip(main_info, target_urls))
print(citations_with_links)

save_path = "/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/"
data_completeName = os.path.join(save_path, ('coge_citations'+'.pickle'))  #with the query for a name
pickle.dump( citations_with_links, open( data_completeName, "wb" ) )
print("dumped to pickle")

