from content_management import *
from Entrez_IR import *

#code for updating tables in pmid_info.db

def recurateInputPapers(user_input):
    #check for user_input's pmcid. if pmcid exists, scrape it
    self_pmcid, abstract, article = getSelfText(user_input)
    #put pmcid, abstract check, and article check into db
    conn, c = connection()
    c.execute("UPDATE inputPapers SET abstract=?, whole_article=?, pmcid =? WHERE pmid=?",(abstract, article, self_pmcid, user_input))  # put user pmid into db
    conn.commit()
    #re-scrape citing papers with new file saving parameters
    pmc_ids = getCitationIDs(user_input)
    all_abstract_check, all_article_check = getContentPMC(pmc_ids)


def testLock(testid):
    conn, c = connection()
    new = 'bbb'
    c.execute("UPDATE annotations SET pmcid=? WHERE pmcid=?",(new, testid))  # put user pmid into db
    conn.commit()


# def recurateAll(user_input):
#   self_info, new_info, target_journals, target_dates, num_citations = run_IR_not_db(user_input)
#   conn, c = connection()
#
#   #update inputPapers for num_citations
#   for tup in self_info:
#       conn, c = connection()
#       c.execute("UPDATE inputPapers SET num_citations = ? WHERE pmid=?", (num_citations, user_input)) #put user pmid into db
#       conn.commit()
#   logging.info("UPDATING self_info to inputPapers db")
#
#   data_samples, nes_list, total_sentences, sum_tokens = do_SOME_multi_preprocessing(user_input)
#
#   logging.info(total_sentences)
#   logging.info(sum_tokens)
#
#
#   #if there are no new citations:
#   if num_citations == len(total_sentences):
#       i = 0
#       #update citations for abstract, whole_article, sents, tokens
#       for tup in new_info:
#           logging.info("TUP IN MAIN: ")
#           logging.info(tup)
#           pmcid = tup[0]
#           title = tup[1]
#           s = ', '
#           author = str(s.join(tup[2]))
#           journal = tup[3]
#           pubdate = tup[4]
#           url = tup[5]
#           abstract = tup[6]
#           whole = tup[7]
#           unix = time.time()
#           date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
#
#           logging.info("sents, tokens: ")
#           sents = total_sentences[i]
#           logging.info(sents)
#           tokens = sum_tokens[i]
#           logging.info(tokens)
#
#           # conn, c = connection()
#           # c.execute("UPDATE citations SET abstract=?, whole_article=?, sents =?, tokens=? WHERE pmcid=?", (abstract, whole, sents, tokens, pmcid) )
#           # conn.commit()
#
#
#           #c.execute("SELECT * FROM citations WHERE pmcid = (?)", (pmcid, ))
#           entry = c.fetchone()
#
#         #if the entry doesn't exist for some reason, create it (documents scraped and annotated earlier)
#           if entry is None:
#               logging.info("No entry found")
#               c.execute("INSERT INTO citations (datestamp, pmcid, title, author, journal, pubdate, citesPmid, url, abstract, whole_article, sents, tokens) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (date, pmcid, title, author, journal, pubdate, user_input, url, abstract, whole, sents, tokens)) #put user pmid into db
#               conn.commit()
#         #if the entry already exists, just up date it
#           else:
#               logging.info("Entry found")
#               conn, c = connection()
#               c.execute("UPDATE citations SET abstract=?, whole_article=?, sents =?, tokens=? WHERE pmcid=?", (abstract, whole, sents, tokens, pmcid) )
#               conn.commit()
#
#
#           i += 1
#       logging.info("UPDATING new_info to citations db")
#
#   #if there are more citations than annotated documents
#   if num_citations != len(total_sentences):
#       logging.info("gotta re-annotated the documents again :( ")
#       data_samples, nes_list, total_sentences, sum_tokens = do_ALL_multi_preprocessing(user_input)
#       i = 0
#       #update citations for abstract, whole_article, sents, tokens
#       for tup in new_info:
#           logging.info("TUP IN MAIN: ")
#           logging.info(tup)
#           pmcid = tup[0]
#           title = tup[1]
#           s = ', '
#           author = str(s.join(tup[2]))
#           journal = tup[3]
#           pubdate = tup[4]
#           url = tup[5]
#           abstract = tup[6]
#           whole = tup[7]
#           unix = time.time()
#           date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
#
#           logging.info("sents, tokens: ")
#           sents = total_sentences[i]
#           logging.info(sents)
#           tokens = sum_tokens[i]
#           logging.info(tokens)
#
#           c.execute("SELECT * FROM citations WHERE pmcid = (?) AND citesPmid = (?)", (pmcid, user_input, ))
#           entry = c.fetchone()
#
#         #if its a new entry, add everything
#           if entry is None:
#               logging.info("No entry found")
#               c.execute("INSERT INTO citations (datestamp, pmcid, title, author, journal, pubdate, citesPmid, url, abstract, whole_article, sents, tokens) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (date, pmcid, title, author, journal, pubdate, user_input, url, abstract, whole, sents, tokens)) #put user pmid into db
#               conn.commit()
#
#         #otherwise just update it
#           else:
#               logging.info("Entry found")
#               conn, c = connection()
#               c.execute("UPDATE citations SET abstract=?, whole_article=?, sents =?, tokens=? WHERE pmcid=?", (abstract, whole, sents, tokens, pmcid) )
#               conn.commit()
#
#           i += 1
#
#
#       logging.info("UPDATING new_info to citations db")











