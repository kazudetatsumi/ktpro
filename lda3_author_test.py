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
import pandas as pd
from sklearn.cluster import KMeans
import nltk
from sklearn.metrics.pairwise import cosine_similarity
from collections import defaultdict, Counter
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
from gensim.models import LdaModel, AuthorTopicModel, CoherenceModel
from gensim.corpora import Dictionary
from wordcloud import WordCloud

#nltk.download('stopwords')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']


def _preprocess(textstring):
    stops = set(stopwords.words('english'))
    textstring = str(textstring).replace('-', ' ')
    tokens = word_tokenize(textstring)

    # Since nltk stopwords is only lowercase, token should be forced to be
    # lowercase.
    # return [token.lower() for token in tokens if token.isalpha() and token
    #         not in stops]

    return [token.lower() for token in tokens if token.isalpha() and
            token.lower() not in stops]


def preprocess(textstring):
    stops = set(stopwords.words('english'))
    # 1. on the whole sentences, words in idioms are coonected with underscore.
    text = str(textstring).lower()
    text = text.replace('monte carlo', 'monte_carlo')
    text = text.replace('monte-carlo', 'monte_carlo')
    text = text.replace('small angle', 'small_angle')
    text = text.replace('powder diffraction', 'powder_diffraction')
    text = text.replace('single crystal', 'single_crystal')
    text = text.replace('single-crystal', 'single_crystal')
    text = text.replace('pair distribution', 'pair_distribution')
    text = text.replace('unit cell', 'unit_cell')
    text = text.replace('unit-cell', 'unit_cell')
    text = text.replace('lattice parameter', 'lattice_parameter')
    text = text.replace('neutron diffraction', 'neutron_diffraction')
    text = text.replace('neutron scattering', 'neutron_scattering')
    text = text.replace('crystal structure', 'crystal_structure')
    text = text.replace('hydrogen storage', 'hydrogen_storage')
    text = text.replace('hydrogen absorption', 'hydrogen_absorption')
    text = text.replace('hydrogen atom', 'hydrogen_atom')
    text = text.replace('time of flight', 'time_of_flight')
    text = text.replace('time-of-flight', 'time_of_flight')
    text = text.replace('cross section', 'cross_section')
    text = text.replace('cross-section', 'cross_section')
    text = text.replace('crosssection', 'cross_section')
    text = text.replace('cross sections', 'cross_sections')
    text = text.replace('cross-sections', 'cross_sections')
    text = text.replace('crosssections', 'cross_sections')
    text = text.replace('prompt gamma', 'prompt_gamma')
    text = text.replace('prompt-gamma', 'prompt_gamma')
    text = text.replace('space group', 'space_group')
    text = text.replace('space-group', 'space_group')
    text = text.replace('phase transition', 'phase_transition')
    text = text.replace('in situ', 'in_situ')
    text = text.replace('in-situ', 'in_situ')
    text = text.replace('x ray', 'x_ray')
    text = text.replace('x-ray', 'x_ray')
    # 2. remaining '-' is replaced with space and divided the text into words
    #    (tokeniation).
    text = text.replace('-', ' ')
    tokens = word_tokenize(text)
    # 3. filtering
    cleaned_tokens = []
    for token in tokens:
        # keep words including underscores or not included in stopss.
        if (token.isalpha() or '_' in token) and token not in stops:
            # words whose number of letters <= 2 are excluded.
            if len(token) > 2:
                cleaned_tokens.append(token)
    return cleaned_tokens


def preprocess2(textstring):
    """
    remove_word = ['c', 'elsevier', 'science', 'all', 'rights', 'reserved',
                   'a', 'the', 'it', 'we', 'paper', 'space', 'astronomy',
                   'astronomical', 'astrophysics', 'astrophysical', 'using',
                   'results', 'similar', 'based', 'used', 'method',
                   'analysis', 'study', ]
    """
    remove_word = ['c', 'elsevier', 'science', 'all', 'rights', 'reserved',
                   'a', 'the', 'it', 'we', 'paper', 'using', 'results',
                   'similar', 'based', 'used', 'method', 'analysis', 'study',
                   'high', 'also', 'data', 'energy', 'measured', 'degrees',
                   'two', 'model', 'experimental', 'find', 'samples',
                   'different', 'measurements', 'observed', 'calculations',
                   'system', 'range', 'found', 'obtained', 'studied', 'one',
                   'determined', 'investigated', 'experiments', 'new', 'well',
                   'treatment', 'studies', 'code', 'developed',
                   'design', 'research', 'equation', 'observations', 'models',
                   'properties', 'show', 'may', 'effect', 'calculated',
                   'compared', 'values', 'methods', 'one', 'numerical',
                   'induced', 'measurement', 'performance', 'current',
                   'neutron', 'neutrons', 'sample', 'nan', 'mathrsfs',
                   'amsfonts', 'amsmath', 'document', 'time', 'due', 'first',
                   'within', 'present', 'order', 'type', 'process', 'rate',
                   'level', 'problems', 'low', 'small', 'fast', 'minimal',
                   'increase', 'materials', 'material', 'elements',
                   'calculation', 'parameters', 'distribution', 'effective',
                   'performed', 'presented', 'proposed', 'approach', 'work',
                   'experiment', 'however', 'three', 'total', 'evaluated',
                   'along', 'general', 'large', 'higher', 'accuracy',
                   'technique', 'techniques', 'facility', 'function',
                   'effects', 'systems', 'structural', 'structures', 'rates',
                   'temperatures', 'sources', 'states', 'site', 'sites',
                   'wasysym', 'amssymb', 'upgreek', 'amsbsy', 'dot', 'circle',
                   'nss', 'use', 'significant', 'mean', 'beam', 'source',
                   'sources', 'beams', 'reactor', 'target', 'pulse',
                   'detector', 'showed', 'density', 'content', 'flux',
                   'condition', 'conditions', 'contents', 'region',
                   ]
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


def prepare_quantities_for_authors(
    df_pre,
):
    frequency = defaultdict(int)

    # count the number of occurrences of the word
    author_list = []
    for text in df_pre.Author:
        #print(text)
        if not pd.isna(text):
            # 修正1：ピリオドが消えないように '.,' を '.;' に置換
            #text = text.replace('.,', ';').replace(' and', ';').replace(
            #        ' ', '').split(';')
            text = text.replace('.,', '.;').replace(' and', ';').replace(
                    ' ', '').split(';')
            author_list.append(text)
            for token in text:
                frequency[token] += 1
        # 修正2：著者が空欄でも、行ズレを防ぐために空のリストを追加する！（最重要）
        else:
            author_list.append([])


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
    print(savefile)
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
):
    num_topics = model.num_topics
    #ncols = 5
    #nrows = math.ceil(num_topics/ncols)
    nrows = 5
    ncols = math.ceil(num_topics/nrows)
    #fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(10, 8))
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows,)
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
        axs[i].set_title('#'+str(t))
    fig.tight_layout()
    #plt.subplots_adjust(wspace=-0.2, hspace=0.25)
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
        print(f"\n{'='*15} Topic {i} {'='*15}")
        selected_phis = model.show_topic(i, topn=topn)
        df_phi = pd.DataFrame(selected_phis, columns=['word', 'probability'])
        print(df_phi.set_index('word'))


def get_theta_matrix(
        df_pre,
        model,
        corpus,
        ):
    """
    calculate the whole set of theta_d and return it as thetas.
    """
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


def get_info_on_specific_author(
        authname,
        author2doc,
        model,
        theta_a_list,
        df_pre,
        theta_a,
):
    # AIの直言
    # author2docではなく、model由来の author_names からインデックスを探す
    # authidx = list(author2doc.keys()).index(authname)
    author_names = [model.id2author[i] for i in range(len(model.id2author))]
    authidxg = author_names.index(authname)
    for i in range(len(theta_a_list[authidxg])):
        print(theta_a_list[authidxg][i])
    #for i in range(len(theta_a_list[authidx])):
    #    print(theta_a_list[authidx][i])
    for i in author2doc[authname]:
        print(i, df_pre.Title[i])
    cos_sim = np.sum(theta_a[authidxg]*theta_a, axis=1)\
        / (np.linalg.norm(theta_a[authidxg]) *
           np.linalg.norm(theta_a, axis=1))
    deg = np.degrees(np.arccos(np.clip(cos_sim, -1, 1)))
    # AIの直言
    # author2doc.keys() を使わず、model の ID順に名前を並べる
    # author_list = list(author2doc.keys())
    author_list = [model.id2author[i] for i in range(len(model.id2author))]
    for i in range(30):
        print(i, np.argsort(deg)[i],
              deg[np.argsort(deg)[i]],
              author_list[np.argsort(deg)[i]],
              len(author2doc[author_list[np.argsort(deg)[i]]]))


def get_2D_map(theta_a, model, df_cluster):
    X = theta_a
    # list of 22,862 authors in correct order
    names_list = [model.id2author[i] for i in range(len(model.id2author))]
    # PCA (SVD)
    # Centering data
    X_mean = np.mean(X, axis=0)
    X_centered = X - X_mean
    # SVD
    U, s, Vh = np.linalg.svd(X_centered, full_matrices=False)
    # 2D coordination of (PC1, PC2)
    coords = U[:, :2] * s[:2]
    pc1 = coords[:, 0]
    pc2 = coords[:, 1]
    # contribution of PC1 and PC2
    var_exp = (s**2) / np.sum(s**2) * 100
    # Visualization
    plt.figure(figsize=(12, 10))
    plt.rcParams['font.family'] = 'Arial'
    # Plot all authors with thin gray
    plt.scatter(pc1, pc2, s=5, alpha=0.15, c='gray', edgecolors='none',
                label='All Researchers')
    #plt.scatter(pc1, pc2, s=5, alpha=0.15, c=df_cluster['Cluster'], cmap='Set1', edgecolors='none',
    #            label='All Researchers')
    # Plot of knwon authors
    targets = {
        'Kawakita,Y.': 'red',
        'Kitaguchi,M.': 'blue',
        'Kanaya,T.': 'green',
    }
    for name, color in targets.items():
        if name in names_list:
            idx = names_list.index(name)
            plt.scatter(pc1[idx], pc2[idx], c=color, s=100, marker='*',
                        edgecolors='black', linewidth=1, label=name)
    # Mean position of each topics shown as an arrow
    # indicates what topic you will be close to as you proceed along a specifc
    # direction.
    scale_factor = np.max(np.abs(coords)) * 1.2
    for i in range(5):
        # Vh[0, i] is contirubtion to PC1, and Vh[1, i] is contribution to PC2.
        plt.arrow(0, 0, Vh[0, i] * scale_factor, Vh[1, i] * scale_factor,
                  color='k', alpha=0.8, width=scale_factor*0.003,
                  head_width=scale_factor*0.01)
        plt.text(Vh[0, i] * scale_factor * 1.1, Vh[1, i] * scale_factor * 1.1,
                 f"Topic {i}", color='k', fontsize=12,
                 fontweight='bold')
    plt.xlabel(f'PC1 ({var_exp[0]:.1f}%)')
    plt.ylabel(f'PC2 ({var_exp[1]:.1f}%)')
    plt.title('Researcher Map (Author Topic based PCA and kmeans)')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.axhline(0, color='black', alpha=0.2)
    plt.axvline(0, color='black', alpha=0.2)
    plt.show()


def cluster_authors_by_theta_a(theta_a, author_names, k=6, random_state=42):
    """
    theta_a (トピック分布行列) を用いてK-meansクラスタリングを行い、
    全著者にラベルを付与して DataFrame として返す関数。
    """
    print(f"Running K-means clustering with k={k}...")

    # 1. K-meansの実行 (n_init=10 はアルゴリズムの安定化のため)
    kmeans = KMeans(n_clusters=k, random_state=random_state, n_init=10)
    cluster_labels = kmeans.fit_predict(theta_a)

    # 2. 結果を Pandas DataFrame に整理する
    # トピック確率も一緒に保存しておくと後で検索が超便利になります
    num_topics = theta_a.shape[1]
    df_authors = pd.DataFrame(theta_a, columns=[f"Topic_{i}" for i in range(num_topics)])
    df_authors['Author'] = author_names
    df_authors['Cluster'] = cluster_labels

    # 列の並び順を見やすく変更
    cols = ['Author', 'Cluster'] + [f"Topic_{i}" for i in range(num_topics)]
    df_authors = df_authors[cols]

    # 3. クラスタごとのサマリー（人口と特徴）を表示
    print("\n--- Cluster Summary ---")
    cluster_centers = kmeans.cluster_centers_
    for i in range(k):
        # その村の人数
        count = np.sum(cluster_labels == i)
        # その村の「重心（平均的な住人）」が最も強く持っているトピック
        dominant_topic = np.argmax(cluster_centers[i])
        max_prob = cluster_centers[i][dominant_topic]
        print(f"Cluster {i}: {count:5d} authors | Dominant Feature: Topic {dominant_topic} ({max_prob*100:.1f}%)")

    return df_authors, kmeans


def run_kmeans(theta_a, model):
    # ==========================================
    # 実行部分
    # ==========================================
    # k=6 で実行（必要に応じて k=5 や k=8 などに変更して試してください）
    # ※ names_list は前回の PCA の時に作った [model.id2author[i] for i in ...] のリストです
    names_list = [model.id2author[i] for i in range(len(model.id2author))]
    k_num = 6
    df_cluster, kmeans_model = cluster_authors_by_theta_a(
            theta_a, names_list, k=k_num)

    # --- 答え合わせ（3大スターの所属村を確認） ---
    print("\n--- Target Authors Clusters ---")
    targets = ['Kawakita,Y.', 'Kitaguchi,M.', 'Kanaya,T.']

    for target in targets:
        if target in names_list:
            # その人の行を抽出
            author_info = df_cluster[df_cluster['Author'] == target].iloc[0]
            c_label = author_info['Cluster']

            # わかりやすく表示
            print(f"⭐ {target}")
            print(f"   Assigned to : Cluster {c_label}")
            # トピック分布の数値を小数点3桁で表示
            probs = [round(p, 3) for p in author_info.iloc[2:].values]
            print(f"   Theta_a     : {probs}")
    return df_cluster


def run_author_LDA_on_abstract(
    df_pre,
    data,
    savfile: str,
    num_topics: int = 12,
):
    corpus, id2word, texts, dictionary =\
        prepare_necesarry_quantities_for_lda_from_data(data)
    count_number_of_words(corpus, "Abstract")
    author2doc = prepare_quantities_for_authors(df_pre)
    print("CHK1", len(author2doc))

    model = get_model(corpus, id2word, author2doc, savfile,
                      num_topics=num_topics,)
    print("CHK", type(model.id2author))
    print(len(model.id2author))
    print(len(model.id2author.values()))
    print("Kawakita,Y.", model.get_author_topics('Kawakita,Y.', minimum_probability=0))
    print("Shibayama,N.", model.get_author_topics('Shibayama,N.', minimum_probability=0))
    print("Kitaguchi,M.", model.get_author_topics('Kitaguchi,M.', minimum_probability=0))
    print("Yamaguchi,N.", model.get_author_topics('Yamaguchi,N.', minimum_probability=0))
    theta_a_list = [model.get_author_topics(author, minimum_probability=0) for author
                    in model.id2author.values()]
    theta_a = np.zeros((len(theta_a_list), num_topics))
    for i in range(len(theta_a_list)):
        for j in range(len(theta_a_list[i])):
            theta_a[i][theta_a_list[i][j][0]] = theta_a_list[i][j][1]
    return corpus, model, author2doc, theta_a_list, theta_a


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
    author2doc = prepare_quantities_for_authors(df_pre)
    model = get_model(corpus, id2word, author2doc, savfile,
                      num_topics=num_topics,)
    plot_word_cloud(model)
    author_vecs = [model.get_author_topics(author) for author
                   in model.id2author.values()]
    person = 'Kawakita,Y.'
    if target == "Title":
        for i in author2doc[person]:
            print(i, df_pre.Title[i])
        No_author = list(author2doc.keys()).index(person)
        print(author_vecs[No_author])
    if target == "Abstract":
        author_vecs_list = np.zeros((len(author_vecs), num_topics))
        print("author_vecs_list shape:", author_vecs_list.shape)
        for i in range(len(author_vecs)):
            for j in range(len(author_vecs[i])):
                author_vecs_list[i][author_vecs[i][j][0]] =\
                        author_vecs[i][j][1]
        No_author = list(author2doc.keys()).index(person)
        for i in range(len(author_vecs[No_author])):
            print(author_vecs[No_author][i])
        for i in author2doc[person]:
            print(i, df_pre.Title[i])
        cos_sim = np.sum(author_vecs_list[No_author]*author_vecs_list, axis=1)\
            / (np.linalg.norm(author_vecs_list[No_author]) *
               np.linalg.norm(author_vecs_list, axis=1))
        author_list = list(author2doc.keys())
        for i in range(1, 31):
            print(i, np.argsort(cos_sim)[-1*i],
                  round(cos_sim[np.argsort(cos_sim)[-1*i]], 3),
                  author_list[np.argsort(cos_sim)[-1*i]],
                  len(author2doc[author_list[np.argsort(cos_sim)[-1*i]]]))
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
    #corpus_title, model_title = run_author_LDA(
    #    df_pre, df_pre.Title, savtitlefile, 'Title', num_topics,)
    #corpus_abstract, model_abstract = run_author_LDA(
    #    df_pre, df_pre.Abstract, savabstractfile, 'Abstract',
    #    num_topics=num_topics)

    #corpus, model, author2doc, theta_a_list, theta_a =\
    corpus, model, author2doc, theta_a_list, theta_a =\
        run_author_LDA_on_abstract(
            df_pre, df_pre.Abstract, savabstractfile, num_topics=num_topics)
    plot_word_cloud(model)
    ## 辞書のキー（著者名）の中から 'Takahara' を含むものをすべて探す
    #takahara_list = [name for name in author2doc.keys() if 'Kawakita' in name]
#
#    print("--- Found Takahara variations ---")
#    for name in takahara_list:
#        # 名前と、その人が何件の論文を持っているかを表示
#        print(f"'{name}' : {len(author2doc[name])} docs")
    #get_info_on_specific_author(
    #    'Kawakita,Y.',
    #    author2doc,
    #    model,
    #    theta_a_list,
    #    df_pre,
    #    theta_a,
    #    )
    #show_top_phis(model, topn=30)
    df_cluster = run_kmeans(theta_a, model)
    get_2D_map(theta_a, model, df_cluster)
    #corpus_abstract, model_abstract = run_author_LDA(
    #    df_pre, df_pre.Abstract, savabstractfile, "Abstract", num_topics=num_topics)
    #show_top_phis(model_abstract)
    """
    thetas = get_theta_matrix(df_pre, model_abstract, corpus_abstract,)
    show_title_of_samples_with_10_largest_thetas_on_each_topic(
        df_pre, thetas, model_abstract,)
    df_country_year, df_merge = show_journal_info(df_pre, debug=debug)
    plot_countries(df_country_year, df_merge)
    """
