from content_management import *

def recurate(user_input):
  self_info, new_info, target_journals, target_dates, num_citations = run_IR_not_db(user_input)
  conn, c = connection()

  #update inputPapers for num_citations
  for tup in self_info:
      conn, c = connection()
      c.execute("UPDATE inputPapers SET num_citations = ? WHERE pmid=?", (num_citations, user_input)) #put user pmid into db
      conn.commit()
  logging.info("UPDATING self_info to inputPapers db")

  data_samples, nes_list, total_sentences, sum_tokens = do_SOME_multi_preprocessing(user_input)

  logging.info(total_sentences)
  logging.info(sum_tokens)


  i = 0
  #update citations for abstract, whole_article, sents, tokens
  for tup in new_info:
      logging.info("TUP IN MAIN: ")
      logging.info(tup)
      pmcid = tup[0]
      abstract = tup[6]
      whole = tup[7]

      logging.info("sents, tokens: ")
      sents = total_sentences[i]
      logging.info(sents)
      tokens = sum_tokens[i]
      logging.info(tokens)

      conn, c = connection()
      c.execute("UPDATE citations SET abstract=?, whole_article=?, sents =?, tokens=? WHERE pmcid=?", (abstract, whole, sents, tokens, pmcid) )
      conn.commit()
      i += 1
  logging.info("UPDATING new_info to citations db")


recurate("25884402")

