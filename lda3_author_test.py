#!/usr/bin/env python
"""
This script performes author Latent Dirichlet Allocation (LDA) analysis
on a set of titles and abstract of selected papers.
Kazuyoshi Tatsumi transferred the ipynb script written by Dr. Yoshifumi Amano
to this script.
The input documents should be prepared as a csv file by a separated script.
Kazuyoshi TATSUMI 2026/05/17
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import math
import ast
import os
import time
import requests
import nltk
from sklearn.metrics.pairwise import cosine_similarity
from collections import defaultdict, Counter
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
from gensim.models import LdaModel, AuthorTopicModel, CoherenceModel
from gensim.corpora import Dictionary
from wordcloud import WordCloud

nltk.download('stopwords')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']


def preprocess(textstring):
    stops = set(stopwords.words('english'))
    textstring = str(textstring).replace('-', ' ')
    tokens = word_tokenize(textstring)

    return [token.lower() for token in tokens if token.isalpha() and token
            not in stops]


def preprocess2(textstring):
    remove_word = ['c', 'elsevier', 'science', 'all', 'rights', 'reserved',
                   'a', 'the', 'it', 'we', 'paper', 'space', 'astronomy',
                   'astronomical', 'astrophysics', 'astrophysical']
    return [textstring[j] for j in range(len(textstring)) if textstring[j]
            not in remove_word]


def do_initial_setup():
    nltk.download('punkt_tab')
    df = pd.read_csv('Neutron.csv')
    df.Abstract = df.Abstract.apply(preprocess).apply(preprocess2)
    df.Title = df.Title.apply(preprocess).apply(preprocess2)
    df.to_csv('df_Neutron.csv')


def plot_number_of_papers_per_year(
    df_pre,
):
    N_paper_year = []
    for year in range(1991, 2027):
        N_paper_year.append([year, len(df_pre[df_pre.Year == year])])
    df_N_paper_year = pd.DataFrame(N_paper_year)
    #print(df_N_paper_year.sum(axis=0))
    #print(df_N_paper_year)
    plt.subplots(figsize=(5, 5))
    plt.plot(df_N_paper_year.iloc[:-1, 0], df_N_paper_year.iloc[:-1, 1],
             marker='o', color='black', markersize=8)
    # plt.legend(bbox_to_anchor=(1, 1), fontsize=14)
    plt.xlabel('Year', fontsize=20)
    plt.ylabel('Number of papers', fontsize=20)
    plt.tick_params(which='minor', width=1, length=4)
    plt.tick_params(which='major', labelsize=18, width=1.5, length=6)
    plt.yticks([0, 1000, 2000, 3000])
    plt.minorticks_on()
    plt.tight_layout()
    plt.show()


def prepare_necesarry_quantities_for_lda_from_data(
    data,
):
    frequency = defaultdict(int)

    # count the number of occurrences of the word
    for text in data:
        for token in ast.literal_eval(text):
            frequency[token] += 1

    # build only words above 30 into an array
    texts = [[token for token in ast.literal_eval(text)
              if frequency[token] > 30] for text in data]

    #print(texts[0])
    dictionary = Dictionary(texts)
    dictionary.filter_extremes(no_below=5, no_above=0.8)
    # For fast computation, I modified to use texts instead of data.values.
    # corpus = [dictionary.doc2bow(ast.literal_eval(summary))
    #           for summary in data.values]
    corpus = [dictionary.doc2bow(text) for text in texts]
    #print(corpus)
    id2word = dict(dictionary.items())
    #print(id2word)
    return corpus, id2word, texts, dictionary,


def count_number_of_words(
    corpus,
    target: str,
):
    count_word_num = 0
    for i in range(len(corpus)):
        count_word_num += len(corpus[i])
    print("# of words in "+target+":", count_word_num)


def perpare_quantities_for_authors(
    df_pre,
):
    frequency = defaultdict(int)

    # count the number of occurrences of the word
    author_list = []
    for text in df_pre.Author:
        #print(text)
        if not pd.isna(text):
            text = text.replace('.,', ';').replace(' and', ';').replace(
                    ' ', '').split(';')
            author_list.append(text)
            for token in text:
                frequency[token] += 1

    # build only words above 5 into an array
    authors = [[token for token in text if frequency[token] >= 5]
               for text in author_list]

    # print(frequency)

    author2doc = {}
    for docid, author in enumerate(authors):
        for each_author in author:
            if each_author in author2doc.keys():
                author2doc[each_author].append(docid)
            else:
                author2doc[each_author] = [docid]
    return author2doc


def get_num_of_topics_dependence_of_cp(
    corpus,
    id2word,
    texts,
    dictionary,
    cpfile,
):
    # https://qiita.com/Spooky_Maskman/items/0d03ea499b88abf56819
    coherence_perplexity = []
    coherence_perplexity.append(['n_topic', 'perplexity', 'coherence'])
    for n_topic in range(1, 21):
        lda_model = LdaModel(
            corpus=corpus, id2word=id2word, num_topics=n_topic, random_state=0)
        coherence_model_lda = CoherenceModel(
            model=lda_model, texts=texts,
            dictionary=dictionary, coherence='c_v')
        print(n_topic, np.exp2(-lda_model.log_perplexity(corpus)),
              coherence_model_lda.get_coherence())
        coherence_perplexity.append(
            [n_topic, np.exp2(-lda_model.log_perplexity(corpus)),
             coherence_model_lda.get_coherence()])
    pd.DataFrame(coherence_perplexity).to_csv(cpfile, header=None)


def plot_coherence_and_preplrexity(
    df_coherence_perplexity,
):
    fig = plt.figure()

    ax1 = fig.subplots()
    ax2 = ax1.twinx()

    l1, = ax1.plot(df_coherence_perplexity['n_topic'],
                   df_coherence_perplexity['perplexity'],
                   c='blue', label='perplexity')
    ax1.set_ylabel('perplexity')
    l2, = ax2.plot(df_coherence_perplexity['n_topic'],
                   df_coherence_perplexity['coherence'],
                   c='red', label='coherence')
    ax2.set_ylabel('coherence')
    plt.xticks(np.arange(1, 21, 2))
    plt.legend(handles=[l1, l2])
    plt.show()


def get_model(
    corpus,
    id2word,
    author2doc,
    savefile,
    num_topics=12,
):
    if not os.path.exists(savefile):
        model = AuthorTopicModel(corpus=corpus, id2word=id2word,
                                 author2doc=author2doc, num_topics=num_topics,
                                 random_state=0)
        pickle.dump(model, open(savefile, 'wb'))
    model = pickle.load(open(savefile, 'rb'))
    # top_topics = list(model.top_topics(corpus))
    # print(top_topics)
    return model


def color_func(word, font_size, position, orientation, random_state,
               font_path):
    return 'darkturquoise'


def plot_word_cloud(
    model,
    target: str,
):
    num_topics = model.num_topics
    ncols = 4
    nrows = math.ceil(num_topics/ncols)
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(10, 8))
    axs = axs.flatten()

    for i, t in enumerate(range(model.num_topics)):
        x = dict(model.show_topic(t, 30))
        im = WordCloud(
            background_color='black',
            color_func=color_func,
            max_words=4000,
            width=300, height=300,
            random_state=0
        ).generate_from_frequencies(x)
        axs[i].imshow(im.recolor(colormap='Paired_r', random_state=244),
                      alpha=0.98)
        axs[i].axis('off')
        axs[i].set_title('Topic '+str(t+1), size=25)
    fig.suptitle('LDA, '+target, size=35)
    #fig.tight_layout()
    plt.subplots_adjust(wspace=-0.2, hspace=0.25)
    plt.show()


def get_topics_sum(
        df_pre,
        model,
        corpus,
):
    topics_sum = []
    for year in range(1991, 2024):
        topics_year = np.zeros(model.num_topics)
        count = 0
        for index in df_pre[df_pre['Year'] == year].index:
            count += 1
            topics = model.get_document_topics(corpus[index])
            for i in range(len(topics)):
                topics_year[topics[i][0]] += topics[i][1]
        if count != 0:
            topics_sum.append(np.hstack([year, topics_year/count]))
        #print(year, np.sum(topics_year/count))
    return topics_sum


def plot_year(
        topics_sum,
        model,
):
    cm = plt.get_cmap('tab20')
    plt.subplots(figsize=(8, 3))
    for i in range(1, model.num_topics+1):
        plt.plot(pd.DataFrame(topics_sum)[0],
                 pd.DataFrame(topics_sum)[i], linewidth=5,
                 label='topic'+str(i), c=cm(i))
    plt.legend(bbox_to_anchor=(1, 1), fontsize=14, labelspacing=0.1)
    plt.xlabel('Year', fontsize=20)
    plt.ylabel('Probability', fontsize=20)
    plt.tick_params(which='minor', width=1, length=4)
    plt.tick_params(which='major', labelsize=18, width=1.5, length=6)
    plt.yticks([0.06, 0.1, 0.14])
    plt.xlim(1993.5, 2023.5)
    plt.minorticks_on()
    #plt.title('Title', fontsize=40)
    plt.tight_layout()
    plt.show()


def show_top_phis(model, topn=10):
    for i in range(model.num_topics):
        print(f"\n{'='*15} Topic {i+1} {'='*15}")
        selected_phis = model.show_topic(i, topn=topn)
        df_phi = pd.DataFrame(selected_phis, columns=['word', 'probability'])
        print(df_phi.set_index('word'))


def get_theta_matrix(
        df_pre,
        model,
        corpus,
        ):
    thetas = np.zeros((len(df_pre), model.num_topics))
    for d_idx, crp in enumerate(corpus):
        theta_d = model.get_document_topics(crp, minimum_probability=0)
        for k_idx, prob in theta_d:
            thetas[d_idx][k_idx] = prob
    return thetas


def show_title_of_samples_with_10_largest_thetas_on_each_topic(
        df_pre,
        thetas,
        model,
        ):
    for k_idx in range(model.num_topics):
        print(f"\n{'='*20} Topic {k_idx+1} {'='*20}")
        top10_didxs = np.argsort(thetas[:, k_idx])[::-1][:10]
        for d_idx in top10_didxs:
            print(f"[{d_idx:d}][{thetas[d_idx, k_idx]:.4f}]"
                  + f" {df_pre.iloc[d_idx].Title}")


def prepare_journal_list_file(
        df_pre,
        ):
    issn_dic = Counter(df_pre.ISSN)
    #print(issn_dic)
    df_issn_dic = pd.DataFrame(issn_dic.items(), columns=['ISSN', 'count'])
    df_issn_dic = df_issn_dic.sort_values('count', ascending=False)
    #print(df_issn_dic)
    Journal_list = []
    for i in range(0, len(df_issn_dic)):
        num_paper = df_issn_dic.iloc[i]['count']
        if num_paper < 10:
            break
        issn = df_issn_dic.iloc[i]['ISSN']
        result = requests.get("https://api.openalex.org/sources",
                              params={'filter': f'issn:{issn}'})
        time.sleep(0.2)
        if result.json()['results'] == []:
            continue
        if i < 11 or i % 100 == 0:
            print(result.json()['results'][0]['display_name'],
                  result.json()['results'][0]['country_code'], issn, num_paper)
        Journal_list.append([result.json()['results'][0]['display_name'],
                             result.json()['results'][0]['country_code'],
                             issn, num_paper])
    pd.DataFrame(Journal_list).to_csv('Journal_list.csv')


def plot_countries(
        df_country_year,
        df_merge,
        ):
    cm = plt.get_cmap('tab20')
    plt.subplots(figsize=(10, 5))
    for i in range(10):
        plt.plot(df_country_year.index,
                 df_country_year[df_merge.country.value_counts().keys()[i]],
                 marker='o', color=cm(i), markersize=8,
                 label=df_merge.country.value_counts().keys()[i])
    #plt.legend(bbox_to_anchor=(1, 1), fontsize=14)
    plt.xlabel('Year', fontsize=20)
    plt.ylabel('Ratio of papers', fontsize=20)
    plt.tick_params(which='minor', width=1, length=4)
    plt.tick_params(which='major', labelsize=18, width=1.5, length=6)
    plt.legend(bbox_to_anchor=(1, 1), title='Publishing country', fontsize=14,
               title_fontsize=16)
    #plt.yscale('log')
    plt.yticks([0, 0.1, 0.2, 0.3, 0.4])
    plt.minorticks_on()
    plt.tight_layout()
    plt.show()


def show_journal_info(
        df_pre,
        jlistfile: str = 'Journal_list.csv',
        debug: bool = False,
):
    if not os.path.exists(jlistfile):
        prepare_journal_list_file(df_pre,)
    df_journal = pd.read_csv(jlistfile, index_col=0)
    df_journal.columns = ['Journal_from_openarex', 'country', 'ISSN', 'num_paper']
    df_merge = pd.merge(df_pre, df_journal[['country', 'ISSN']], on='ISSN',
                        how='left')
    df_country_year = pd.DataFrame()
    for year in range(1991, 2024):
        df_country_year = pd.concat(
            [df_country_year,
             df_merge[df_merge.Year == year].country.value_counts() /
             np.sum(df_merge[df_merge.Year == year].country.value_counts())
             ], axis=1)
    df_country_year = df_country_year.T
    df_country_year.index = np.arange(1991, 2024)
    if debug:
        print("print(df_journal)")
        print(df_journal)
        print("print(df_merge)")
        print(df_merge)
        print("print(df_contry_year)")
        print(df_country_year)
        country_list = Counter(df_journal.country)
        print("print(country_list)")
        print(country_list)
        df_country = pd.DataFrame(country_list, index=['abc']
                                  ).T.sort_values('abc', ascending=False)
        print("print(df_country)")
        print(df_country)
        print("print(df_country.index[:8])")
        print(df_country.index[:8])
        print("print(df_merge[df_merge.country == 'JP']" +
              ".Journal.value_counts().head(10))")
        print(df_merge[df_merge.country == 'JP'
                       ].Journal.value_counts().head(10))
    return df_country_year, df_merge


def get_df_pre(
    dffile: str,
    debug: bool = True,
):
    if not os.path.exists(dffile):
        do_initial_setup()
    df_pre = pd.read_csv(dffile)
    if debug:
        plot_number_of_papers_per_year(df_pre)
    return df_pre


def run_author_LDA(
    df_pre,
    data,
    #cpfile: str,
    savfile: str,
    target: str,
    num_topics: int = 12,
):
    corpus, id2word, texts, dictionary =\
        prepare_necesarry_quantities_for_lda_from_data(data)
    count_number_of_words(corpus, target)
    author2doc = perpare_quantities_for_authors(df_pre)
    #if not os.path.exists(cpfile):
    #    get_num_of_topics_dependence_of_cp(
    #        corpus, id2word, texts, dictionary, cpfile)
    #df_coherence_perplexity = pd.read_csv(cpfile)
    #plot_coherence_and_preplrexity(df_coherence_perplexity)
    model = get_model(corpus, id2word, author2doc, savfile,
                      num_topics=num_topics,)
    plot_word_cloud(model, target)
    #topics_sum = get_topics_sum(df_pre, model, corpus)
    #plot_year(topics_sum, model)
    author_vecs = [model.get_author_topics(author) for author
                   in model.id2author.values()]
    print(len(author_vecs))
    person = 'Takahara,A.'
    if target == "Title":
        for i in author2doc[person]:
            print(i, df_pre.Title[i])
        No_author = list(author2doc.keys()).index(person)
        print(author_vecs[No_author])
    if target == "Abstract":
        author_vecs_list = np.zeros((len(author_vecs), num_topics))
        for i in range(len(author_vecs)):
            for j in range(len(author_vecs[i])):
                author_vecs_list[i][author_vecs[i][j][0]] =\
                        author_vecs[i][j][1]
        No_author = list(author2doc.keys()).index(person)
        for i in range(len(author_vecs[No_author])):
            print(author_vecs[No_author][i])
        for i in author2doc[person]:
            print(i, df_pre.Title[i])
        cos_sim = np.sum(author_vecs_list[No_author]*author_vecs_list, axis=1)/(np.linalg.norm(author_vecs_list[No_author])*np.linalg.norm(author_vecs_list, axis=1))
        author_list = list(author2doc.keys())
        for i in range(1, 31):
            print(i, np.argsort(cos_sim)[-1*i], round(cos_sim[np.argsort(cos_sim)[-1*i]], 3),
                  author_list[np.argsort(cos_sim)[-1*i]], len(author2doc[author_list[np.argsort(cos_sim)[-1*i]]]))
    return corpus, model


def run(
    dffile: str,
    #cptitlefile: str,
    #cpabstractfile: str,
    num_topics: int = 12,
    savtitlefile: str = 'lda_author_title.sav',
    savabstractfile: str = 'lda_author_abstract.sav',
    debug: bool = True,
):
    df_pre = get_df_pre(dffile, debug=debug)
    corpus_title, model_title = run_author_LDA(
        df_pre, df_pre.Title, savtitlefile, 'Title', num_topics,)

    corpus_abstract, model_abstract = run_author_LDA(
        df_pre, df_pre.Abstract, savabstractfile, 'Abstract',
        num_topics=num_topics)
    show_top_phis(model_abstract)
    """
    thetas = get_theta_matrix(df_pre, model_abstract, corpus_abstract,)
    show_title_of_samples_with_10_largest_thetas_on_each_topic(
        df_pre, thetas, model_abstract,)
    df_country_year, df_merge = show_journal_info(df_pre, debug=debug)
    plot_countries(df_country_year, df_merge)
    """
