from database_management import *


total_sentences = [51, 123, 80, 186, 27, 18]
sum_tokens = [2383, 4622, 3045, 5434, 755, 535]
new_info = [('3659593', 'Korean Red Ginseng Saponin Fraction Downregulates Proinflammatory Mediators in LPS Stimulated RAW264.7 Cells and Protects Mice against Endotoxic Shock', ['Yayeh T', 'Jung KH', 'Jeong HY', 'Park JH', 'Song YB', 'Kwak YS', 'Kang HS', 'Cho JY', 'Oh JW', 'Kim SK', 'Rhee MH'], 'Journal of Ginseng Research', '2012 Jul', 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3659593', 'yes', 'yes'),
('3645680', 'Effect of Opuntia humifusa Supplementation and Acute Exercise on Insulin Sensitivity and Associations with PPAR-γ and PGC-1α Protein Expression in Skeletal Muscle of Rats', ['Kang J', 'Lee J', 'Kwon D', 'Song Y'], 'International Journal of Molecular Sciences', '2013 Mar 28', 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3645680', 'yes', 'yes'),
('3397493', 'Opuntia humifusa Supplementation Increased Bone Density by Regulating Parathyroid Hormone and Osteocalcin in Male Growing Rats', ['Kang J', 'Park J', 'Choi SH', 'Igawa S', 'Song Y'], 'International Journal of Molecular Sciences', '2012 Jun 4', 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3397493', 'yes', 'yes'),
('3140049', 'The Role of Th17 in Neuroimmune Disorders: A Target for CAM Therapy. Part III', ['Vojdani A', 'Lambert J', 'Kellermann G'], 'Evidence-based Complementary and Alternative Medicine : eCAM', '2011 Jun 16', 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3140049', 'yes', 'yes'),
('2596873', 'Quercetin transiently increases energy expenditure but persistently decreases circulating markers of inflammation in C57BL/6J mice fed a high-fat diet', ['Stewart LK', 'Soileau JL', 'Ribnicky D', 'Wang ZQ', 'Raskin I', 'Poulev A', 'Majewski M', 'Cefalu WT', 'Gettys TW'], 'Metabolism: clinical and experimental', '2008 Jul', 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2596873', 'yes', 'yes'),
('2322858', 'Risks and Benefits of Commonly used Herbal Medicines in México', ['Rodriguez-Fragoso L', 'Reyes-Esparza J', 'Burchiel S', 'Herrera-Ruiz D', 'Torres E'], 'Toxicology and applied pharmacology', '2007 Oct 12', 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2322858', 'yes', 'yes')]


print(len(total_sentences))
print(len(sum_tokens))
print(len(new_info))

i = 0
for tup in new_info:
    print(tup)
    pmcid = tup[0]
    print(pmcid)
    sents = total_sentences[i]
    print(sents)
    tokens = sum_tokens[i]
    print(tokens)
    conn, c = connection()
    c.execute("UPDATE citations SET sents=?, tokens=? WHERE pmcid=?", (sents, tokens, pmcid))
    conn.commit()
    i += 1



