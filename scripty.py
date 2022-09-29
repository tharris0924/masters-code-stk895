#%%


from nltk.corpus import brown
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#%%

def initial_clean(text):
    """
    Function to clean text of websites, email addressess and any punctuation
    We also lower case the text
    Args:
        text: raw corpus

    Returns: tokenized corpus

    """
    text = re.sub(r"((\S+)?(http(s)?)(\S+))|((\S+)?(www)(\S+))|((\S+)?(\@)(\S+)?)", " ", text)
    text = re.sub(r"[^a-zA-Z ]", "", text)
    text = text.lower()
    #text = strip_short(text, minsize=3)
    #text = nltk.word_tokenize(text)
    if text.strip():
        return text

#join sentences into docs
docs = []
for fileid in brown.fileids():
    sents = []
    for sent in brown.sents(fileids=fileid):
        sents.extend([initial_clean(s) for s in sent]) #do basic pre-processing

    docs.append(sents)


#%%

#make dataframe to write results to
df = pd.DataFrame()
df_norm = pd.DataFrame()

said_count = []
showed_count = []
kennedy_count = []

said_norm = []
showed_norm = []
kennedy_norm = []

d_len = []
#three_word_len = []
N = 10
for d in docs:
    said_count.append(d.count("said"))
    showed_count.append(d.count("showed"))
    kennedy_count.append(d.count("kennedy"))
    d_len.append(len(d)) #document length
    three_word_len = d.count("said") + d.count("showed") + d.count("kennedy")
    if three_word_len > 0:
        said_norm.append(np.round((N * d.count("said"))/three_word_len))
        showed_norm.append(np.round((N * d.count("showed"))/three_word_len))
        kennedy_norm.append(np.round((N * d.count("kennedy"))/three_word_len))

df.loc[:,'N'] = d_len
df.loc[:,'said'] = said_count
df.loc[:,'showed'] = showed_count
df.loc[:,'kennedy'] = kennedy_count
df_norm.loc[:,'said_norm'] = said_norm
df_norm.loc[:,'showed_norm'] = showed_norm
df_norm.loc[:,'kennedy_norm'] = kennedy_norm
df_norm["N"] = df_norm.sum(axis=1)

#%%

print(df_norm)

#%%

df.to_csv('brown_word_counts.csv')
df_norm.to_csv('brown_normalized_counts.csv')

#%% md

## Make new normalized corpus

#%%

print(d)

#%%

docs_3words= []
l=['said', 'showed', 'kennedy']
doc_lens = []
said_count = []
showed_count = []
kennedy_count = []
for d in docs:
    words = [x for i, x in enumerate(d) if x in l]
    doc_lens.append(len(words))
    said_count.append(words.count("said"))
    showed_count.append(words.count("showed"))
    kennedy_count.append(d.count("kennedy"))
    #docs_3words.append([(x,words.count(x)) for x in set(words)])

#%%

print(docs_3words[0], doc_lens[0])

#%%

N = 10
docs_normalized = []
said_count = []
showed_count = []
kennedy_count = []

for d in range(len(docs_3words)):
    new_doc = []
    for word in docs_3words[d]:
        word_count = word[1]
        new_wordcount = (word_count * N)/doc_lens[d]
        if word[0] == 'said':
            new_word = (word[0], np.round(new_wordcount))
        new_doc.append(new_word)
    if new_doc: docs_normalized.append(new_doc)


#%%

print(docs_normalized)
