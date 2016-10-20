import pickle
from content_management import new_citations_from_db

main_info = []
target_urls = []



main1, db_urls1 = new_citations_from_db(18952863)
for mi in main1:
    main_info.append(mi)
for url in db_urls1:
    target_urls.append(url)

main2, db_urls2 = new_citations_from_db(18269575)
for mi in main2:
    main_info.append(mi)
for url in db_urls2:
    target_urls.append(url)


citations_with_links = list(zip(main_info, target_urls))



pickle.dump(citations_with_links, open( '/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_citations.pickle', "wb" ) )
print("dumped to pickle")