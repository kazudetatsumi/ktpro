#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import copy
import math
import ast
from collections import defaultdict, Counter
import seaborn as sns
import re
import requests

import nltk
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
nltk.download('stopwords')

from gensim.models import LdaModel, CoherenceModel, TfidfModel, AuthorTopicModel
from gensim.corpora import Dictionary
from pprint import pprint

from wordcloud import WordCloud
import gensim
import os

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
    # 最初に実行する。時間がかかるので、それ以降は保存したデータから読み込む。
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
    print(df_N_paper_year.sum(axis=0))
    print(df_N_paper_year)
    plt.subplots(figsize=(5, 5))
    plt.plot(df_N_paper_year.iloc[:-1, 0], df_N_paper_year.iloc[:-1, 1],
             marker='o', color='black', markersize=8)
    plt.legend(bbox_to_anchor=(1, 1), fontsize=14)
    plt.xlabel('Year', fontsize=20)
    plt.ylabel('Number of papers', fontsize=20)
    plt.tick_params(which='minor', width=1, length=4)
    plt.tick_params(which='major', labelsize=18, width=1.5, length=6)
    plt.yticks([0, 1000, 2000, 3000])
    plt.minorticks_on()
    plt.tight_layout()
    plt.show()


def perpare_necessary_quantities_for_lda_on_titles(
    df_pre,
):
    # dfの時はast.literal_evalが不要 (pandasのcsvファイル経由による処理)
    frequency_title = defaultdict(int)

    # count the number of occurrences of the word
    for text in df_pre.Title:
        for token in ast.literal_eval(text):
            frequency_title[token] += 1

    # build only words above 30 into an array
    texts_title = [[token for token in ast.literal_eval(text)
                   if frequency_title[token] > 30] for text in df_pre.Title]
    #print(texts_title)
    dictionary_title = Dictionary(texts_title)
    dictionary_title.filter_extremes(no_below=5, no_above=0.8)
    corpus_title = [dictionary_title.doc2bow(ast.literal_eval(summary))
                    for summary in df_pre.Title.values]
    #print(corpus_title)
    id2word_title = dict(dictionary_title.items())
    #print(id2word_title)
    #corpusにtfidfという指標もあるが、今回は使わない
    #tfidf_title = TfidfModel(corpus_title)
    #corpus_tfidf_title = tfidf_title[corpus_title]
    #corpus_tfidf_title

    #dfの時はast.literal_evalが不要 (pandasのcsvファイル経由による処理)
    return corpus_title, id2word_title, texts_title, dictionary_title,


def perpare_necesarry_quantities_for_lda_on_abstracts(
    df_pre,
):
    frequency = defaultdict(int)

    # count the number of occurrences of the word
    for text in df_pre.Abstract:
        for token in ast.literal_eval(text):
            frequency[token] += 1

    # build only words above 30 into an array
    texts_abstract = [[token for token in ast.literal_eval(text)
                      if frequency[token] > 30] for text in df_pre.Abstract]

    #print(texts_abstract[0])
    dictionary_abstract = Dictionary(texts_abstract)
    dictionary_abstract.filter_extremes(no_below=5, no_above=0.8)
    corpus_abstract = [dictionary_abstract.doc2bow(ast.literal_eval(summary))
                       for summary in df_pre.Abstract.values]
    #print(corpus_abstract)
    id2word_abstract = dict(dictionary_abstract.items())
    #print(id2word_abstract)
    return corpus_abstract, id2word_abstract, texts_abstract,\
        dictionary_abstract,


def count_number_of_words(
    corpus_title,
    corpus_abstract,
):
    #タイトルにおける単語数
    count_word_num_title = 0
    for i in range(len(corpus_title)):
        count_word_num_title += len(corpus_title[i])
    print("# of words in titles:", count_word_num_title)
    #アブストラクトにおける単語数
    count_word_num = 0
    for i in range(len(corpus_abstract)):
        count_word_num += len(corpus_abstract[i])
    print("# of words in abstracts", count_word_num)
    #dfの時はast.literal_evalが不要 (pandasのcsvファイル経由による処理)


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


def get_num_of_topics_dependence_of_coherence_and_preplexity_for_lda_title(
    corpus_title,
    id2word_title,
    texts_title,
    dictionary_title,
    cptitlefile,
):
    # https://qiita.com/Spooky_Maskman/items/0d03ea499b88abf56819
    coherence_perplexity_title = []
    coherence_perplexity_title.append(['n_topic', 'perplexity', 'coherence'])

    for n_topic in range(1, 21):
        lda_model = LdaModel(corpus=corpus_title, id2word=id2word_title,
                             num_topics=n_topic, random_state=0)
        coherence_model_lda = CoherenceModel(
                model=lda_model, texts=texts_title,
                dictionary=dictionary_title, coherence='c_v')
        print(n_topic, np.exp2(-lda_model.log_perplexity(corpus_title)),
              coherence_model_lda.get_coherence())
        coherence_perplexity_title.append(
                [n_topic, np.exp2(-lda_model.log_perplexity(corpus_title)),
                    coherence_model_lda.get_coherence()])
    pd.DataFrame(coherence_perplexity_title).to_csv(cptitlefile)


def get_coherence_perplexity_title_from_file(cptitlefile):
    # 1. 保存時についた余計な列(0,1,2..)をスキップして読み込む
    df_read = pd.read_csv(cptitlefile, index_col=0, skiprows=[0], header=None)
    # 2. 1行目（見出しの文字列）をリストとして取り出す
    header = df_read.iloc[0].tolist()
    # 3. 2行目以降（本来のデータ）を「数値（float）」に変換してリスト化
    data = df_read.iloc[1:].astype(float).values.tolist()
    # 4. くっつけると、計算直後と「全く同じリスト」が完成します
    coherence_perplexity_title = [header] + data
    return coherence_perplexity_title


def plot_coherence_and_preplrexity_of_lda_model(
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
    savefile,
    num_topics=12,
):
    model = LdaModel(corpus=corpus, id2word=id2word,
                     iterations=400, num_topics=num_topics, random_state=0)
    pickle.dump(model, open(savefile, 'wb'))
    # top_topics = list(model.top_topics(corpus))
    # print(top_topics_abstract)
    return model


def plot_word_cloud_on_model_title(
    model_title
):
    fig, axs = plt.subplots(ncols=3, nrows=math.ceil(model_title.num_topics/3),
                            figsize=(9.5, 15))
    axs = axs.flatten()

    def color_func(
            word, font_size, position, orientation, random_state, font_path):
        return 'darkturquoise'

    for i, t in enumerate(range(model_title.num_topics)):
        x = dict(model_title.show_topic(t, 30))
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
        axs[i].set_title('Topic '+str(t+1), size=17)

    # vis
    fig.suptitle('LDA, Title', size=20)
    fig.tight_layout()
    #fig.tight_layout(rect=[0, 0, 1, 0.93])
    plt.subplots_adjust(wspace=-0.2, hspace=0.25)
    plt.show()

    n_components = 12
    fig, axs = plt.subplots(ncols=4, nrows=math.ceil(n_components/4),
                            figsize=(10, 8))
    axs = axs.flatten()

    for i, t in enumerate(range(model_title.num_topics)):
        x = dict(model_title.show_topic(t, 30))
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
    fig.suptitle('LDA, Title', size=35)
    #fig.tight_layout()
    plt.subplots_adjust(wspace=-0.2, hspace=0.25)
    plt.show()


def plot_word_cloud_on_model_abstract(
    model_abstract,
):
    fig, axs = plt.subplots(
        ncols=3, nrows=math.ceil(model_abstract.num_topics/3),
        figsize=(9.5, 15))
    axs = axs.flatten()

    def color_func(word, font_size, position, orientation, random_state,
                   font_path):
        return 'darkturquoise'

    for i, t in enumerate(range(model_abstract.num_topics)):
        x = dict(model_abstract.show_topic(t, 30))
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
        axs[i].set_title('Topic '+str(t+1), size=35)

    fig.suptitle('LDA, Abstract', size=40)
    fig.tight_layout()
    plt.show()

    n_components = 12
    # WordCloud
    fig, axs = plt.subplots(ncols=4, nrows=math.ceil(n_components/4), figsize=(10, 8))
    axs = axs.flatten()

    #def color_func(word, font_size, position, orientation, random_state, font_path):
    #    return 'darkturquoise'

    for i, t in enumerate(range(model_abstract.num_topics)):
        x = dict(model_abstract.show_topic(t, 30))
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

    fig.suptitle('LDA, Abstract', size=35)
    #fig.tight_layout()
    plt.subplots_adjust(wspace=-0.2, hspace=0.25)
    plt.show()


def get_topics_sum_title(
        df_pre,
        model_title,
        corpus_title,
        num_topics_title=12,
):
    topics_sum_title = []
    for year in range(1991, 2024):
        topics_year = np.zeros(num_topics_title)
        count = 0
        for index in df_pre[df_pre['Year'] == year].index:
            count += 1
            topics = model_title.get_document_topics(corpus_title[index])
            for i in range(len(topics)):
                topics_year[topics[i][0]] += topics[i][1]
        topics_sum_title.append(np.hstack([year, topics_year/count]))
        print(year, np.sum(topics_year/count))
    return topics_sum_title


def get_topics_sum_abstract(
        df_pre,
        model,
        corpus_abstract,
        num_topics=12,
):
    topics_sum = []
    for year in range(1994, 2026):
        topics_year = np.zeros(num_topics)
        count = 0
        for index in df_pre[df_pre['Year'] == year].index:
            count += 1
            topics = model.get_document_topics(corpus_abstract[index])
            for i in range(len(topics)):
                topics_year[topics[i][0]] += topics[i][1]
        topics_sum.append(np.hstack([year, topics_year/count]))
        print(year, np.sum(topics_year/count))
    return topics_sum


def plot_year(
        topics_sum_title,
        num_topics_title=12,
):
    cm = plt.get_cmap('tab20')
    plt.subplots(figsize=(8, 3))
    for i in range(1, num_topics_title+1):
        plt.plot(pd.DataFrame(topics_sum_title)[0],
                 pd.DataFrame(topics_sum_title)[i], linewidth=5,
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
    plt.show()


def get_num_of_topics_dependence_of_cp_for_lda_abstract(
    corpus_abstract,
    id2word_abstract,
    texts_abstract,
    dictionary_abstract,
    cpabstractfile,
):
    # https://qiita.com/Spooky_Maskman/items/0d03ea499b88abf56819
    coherence_perplexity_abstract = []
    coherence_perplexity_abstract.append(['n_topic', 'perplexity', 'coherence']
                                         )
    for n_topic in range(1, 21):
        lda_model = LdaModel(
            corpus=corpus_abstract, id2word=id2word_abstract,
            num_topics=n_topic, random_state=0)
        coherence_model_lda = CoherenceModel(
            model=lda_model, texts=texts_abstract,
            dictionary=dictionary_abstract, coherence='c_v')
        print(n_topic, np.exp2(-lda_model.log_perplexity(corpus_abstract)),
              coherence_model_lda.get_coherence())
        coherence_perplexity_abstract.append(
            [n_topic, np.exp2(-lda_model.log_perplexity(corpus_abstract)),
             coherence_model_lda.get_coherence()])
    pd.DataFrame(coherence_perplexity_abstract).to_csv(cpabstractfile,
                                                       header=None)


def show_largest_10_phis(model):
    for i, t in enumerate(range(model.num_topics)):
        print("=========topic", str(i+1), "==========")
        x = dict(model.show_topic(t, 10))
        print(pd.DataFrame.from_dict(dict(x), orient='index'))


def do_unknown_things(
        df_pre, model,
        corpus_abstract,
        num_topics=12,
):
    topics_prev = np.zeros((len(df_pre), num_topics))
    for i in range(len(corpus_abstract)):
        topics = model.get_document_topics(corpus_abstract[i])
        for j in range(len(topics)):
            topics_prev[i][topics[j][0]] += topics[j][1]

    for topic_No in range(num_topics):
        print('---------------------'+str(topic_No+1)+'---------------------')
        for i in range(10):
            No_doc = np.arange(len(topics_prev))[
                    topics_prev[:, topic_No] ==
                    np.sort(topics_prev[:, topic_No])[-(i+1)]]
            if len(No_doc) != 1:
                print('There are more than two documents', i, j)
            print(No_doc[0], np.sort(topics_prev[:, topic_No])[-(i+1)],
                  df_pre.loc[No_doc[0]].Title)

    issn_dic = Counter(df_pre.ISSN)
    print(issn_dic)
    df_issn_dic = pd.DataFrame(issn_dic, index=['abc']).T.sort_values(
        'abc', ascending=False)
    print(df_issn_dic)
    Journal_list = []
    for i in range(0, len(df_issn_dic)):
        num_paper = df_issn_dic.iloc[i].values[0]
        issn = df_issn_dic.index[i]
        if num_paper < 10:
            break
        r = requests.get(f"https://api.openalex.org/sources", params={'filter': f'issn:{issn}'})
        if r.json()['results'] == []:
            continue
        if i < 11 or i % 100 == 0:
            print(r.json()['results'][0]['display_name'], r.json()['results'][0]['country_code'], issn, num_paper)
        Journal_list.append([r.json()['results'][0]['display_name'], r.json()['results'][0]['country_code'], issn, num_paper])
    pd.DataFrame(Journal_list).to_csv('Journal_list.csv')
    issn = '1538-4357'
    requests.get(f"https://api.openalex.org/sources", params={'filter': f'issn:{issn}'}).json()['results'][0]['country_code']
    df_journal = pd.read_csv('Journal_list.csv', index_col=0)
    df_journal.columns = ['journal', 'country', 'ISSN', 'num_paper']
    print(df_journal)
    df_merge = pd.merge(df_pre, df_journal[['country', 'ISSN']], on='ISSN', how='left')
    print(df_merge)
    df_country_year=pd.DataFrame()
    for year in range(1991, 2024):
        df_country_year = pd.concat(
            [df_country_year,
             df_merge[df_merge.Year == year].country.value_counts() /
             np.sum(df_merge[df_merge.Year == year].country.value_counts())
             ], axis=1)
    df_country_year = df_country_year.T
    df_country_year.index = np.arange(1991, 2024)
    print(df_country_year)
    country_list = Counter(df_journal.country)
    print(country_list)
    df_country = pd.DataFrame(country_list, index=['abc']
                              ).T.sort_values('abc', ascending=False)
    print(df_country)
    print(df_country.index[:8])
    country_name_list = ['United States of America',
                         'United Kingdom',
                         'Netherlands',
                         'Germany',
                         'Japan',
                         'Switzerland',
                         "People's Republic of China",
                         'Russian Federation', ]
    cm = plt.get_cmap('tab20')
    plt.subplots(figsize=(5, 6))
    for i in range(8):
        plt.plot(df_country_year.index,
                 df_country_year[df_merge.country.value_counts().keys()[i]],
                 marker='o', color=cm(i), markersize=8,
                 label=country_name_list[i])
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
    plt.show()


def run(
    dffile,
    cptitlefile,
    cpabstractfile,
):
    if not os.path.exists(dffile):
        do_initial_setup()
    df_pre = pd.read_csv(dffile)
    plot_number_of_papers_per_year(df_pre)
    corpus_title, id2word_title, texts_title, dictionary_title =\
        perpare_necessary_quantities_for_lda_on_titles(df_pre)
    corpus_abstract, id2word_abstract, texts_abstract, dictionary_abstract =\
        perpare_necesarry_quantities_for_lda_on_abstracts(df_pre)
    count_number_of_words(corpus_title, corpus_abstract,)
    if not os.path.exists(cptitlefile):
        get_num_of_topics_dependence_of_coherence_and_preplexity_for_lda_title(
            corpus_title, id2word_title, texts_title, dictionary_title,
            cptitlefile)
    df_coherence_perplexity_title = pd.read_csv(cptitlefile)
    plot_coherence_and_preplrexity_of_lda_model(df_coherence_perplexity_title)
    model_title = get_model(corpus_title, id2word_title, 'lda_title.sav',
                            num_topics=12,)
    plot_word_cloud_on_model_title(model_title)
    topics_sum_title = get_topics_sum_title(
            df_pre, model_title, corpus_title, num_topics_title=12,)
    plot_year(topics_sum_title, num_topics_title=12,)
    if not os.path.exists(cpabstractfile):
        get_num_of_topics_dependence_of_cp_for_lda_abstract(
            corpus_abstract, id2word_abstract, texts_abstract,
            dictionary_abstract, cpabstractfile)
    df_coherence_perplexity_abst = pd.read_csv(cpabstractfile)
    print(df_coherence_perplexity_abst)
    plot_coherence_and_preplrexity_of_lda_model(df_coherence_perplexity_abst)
    model_abstract = get_model(corpus_abstract, id2word_abstract,
                               'lda_abstract.sav', num_topics=12,)
    plot_word_cloud_on_model_abstract(model_abstract)
    show_largest_10_phis(model_abstract)
    topics_sum_abstract = get_topics_sum_abstract(
        df_pre, model_abstract, corpus_abstract, num_topics=12)
    plot_year(topics_sum_abstract, num_topics_title=12,)
    do_unknown_things(df_pre, model_abstract, corpus_abstract, num_topics=12)
