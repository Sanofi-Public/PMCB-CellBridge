import numpy as np
import scipy
import time
from datetime import datetime
from sklearn.decomposition import PCA,TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from fa2 import ForceAtlas2
import networkx as nx
import random


def tot_counts_norm_sparse(E, exclude_dominant_frac = 1, included = [], target_mean = 0):
    '''total counts normalization for sparse matrix'''
    E = E.tocsc()
    ncell = E.shape[0]
    if len(included) == 0:
        if exclude_dominant_frac == 1:
            tots_use = E.sum(axis=1)
        else:
            tots = E.sum(axis=1)
            wtmp = scipy.sparse.lil_matrix((ncell, ncell))
            wtmp.setdiag(1. / tots)
            included = np.asarray(~(((wtmp * E) > exclude_dominant_frac).sum(axis=0) > 0))[0,:]
            tots_use = E[:,included].sum(axis = 1)
            print(('Excluded %i genes from normalization' %(np.sum(~included))))
    else:
        tots_use = E[:,included].sum(axis = 1)

    if target_mean == 0:
        target_mean = np.mean(tots_use)

    w = scipy.sparse.lil_matrix((ncell, ncell))
    w.setdiag(float(target_mean) / tots_use)
    Enorm = w * E

    return Enorm.tocsc(), target_mean, included


def sparse_var(E, axis=0):
    '''helper function for the pca'''
    mean_gene = E.mean(axis=axis).A.squeeze()
    tmp = E.copy()
    tmp.data **= 2
    return tmp.mean(axis=axis).A.squeeze() - mean_gene ** 2


def sparse_multiply(E, a):
    '''helper function for the pca'''
    nrow = E.shape[0]
    w = scipy.sparse.lil_matrix((nrow, nrow))
    w.setdiag(a)
    return w * E


def runningquantile(x, y, p, nBins):
    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]

    dx = (x[-1] - x[0]) / nBins
    xOut = np.linspace(x[0]+dx/2, x[-1]-dx/2, nBins)

    yOut = np.zeros(xOut.shape)

    for i in range(len(xOut)):
        ind = np.nonzero((x >= xOut[i]-dx/2) & (x < xOut[i]+dx/2))[0]
        if len(ind) > 0:
            yOut[i] = np.percentile(y[ind], p)
        else:
            if i > 0:
                yOut[i] = yOut[i-1]
            else:
                yOut[i] = np.nan

    return xOut, yOut


def get_PCA_sparseInput(E, numpc=50, method='sparse', 
                        normalize=True, seed=None):
    # Takes sparse counts matrix as input.
    # If method == 'sparse', gene-level normalization maintains sparsity
    #     (no centering) and TruncatedSVD is used instead of normal PCA.

    base_ix = np.arange(E.shape[0])

    print({'numpc': numpc, 'method' : method, 'normalize' : normalize})

    if method == 'sparse':
        if normalize:
            zstd = np.sqrt(sparse_var(E[base_ix,:]))
            Z = sparse_multiply(E.T, 1 / zstd).T
        else:
            Z = E

        pca = TruncatedSVD(n_components=numpc, random_state=seed)
        pca.fit(Z[base_ix,:])
        return pca.transform(Z)

    else:
        if normalize:
            zmean = E[base_ix,:].mean(0)
            zstd = np.sqrt(sparse_var(E[base_ix,:]))
            Z = sparse_multiply((E - zmean).T, 1/zstd).T
        else:
            Z = E
        pca = PCA(n_components=numpc, random_state=seed)
        pca.fit(Z[base_ix,:])
        return pca.transform(Z)


def get_vscores_sparse(E, min_mean=0, nBins=50, fit_percentile=0.1, error_wt=1):
    '''gets v-scores, which are smoothened dispersion estimates'''

    mu_gene = E.mean(axis=0).A.squeeze()
    gene_ix = np.nonzero(mu_gene > min_mean)[0]
    mu_gene = mu_gene[gene_ix]

    tmp = E[:,gene_ix]
    tmp.data **= 2
    var_gene = tmp.mean(axis=0).A.squeeze() - mu_gene ** 2
    del tmp
    FF_gene = var_gene / mu_gene

    data_x = np.log(mu_gene)
    data_y = np.log(FF_gene / mu_gene)

    x, y = runningquantile(data_x, data_y, fit_percentile, nBins)
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]

    gLog = lambda input: np.log(input[1] * np.exp(-input[0]) + input[2])
    h,b = np.histogram(np.log(FF_gene[mu_gene>0]), bins=200)
    b = b[:-1] + np.diff(b)/2
    max_ix = np.argmax(h)
    c = np.max((np.exp(b[max_ix]), 1))
    errFun = lambda b2: np.sum(abs(gLog([x,c,b2])-y) ** error_wt)
    b0 = 0.1
    b = scipy.optimize.fmin(func = errFun, x0=[b0], disp=False)
    a = c / (1+b) - 1

    v_scores = FF_gene / ((1+a)*(1+b) + b * mu_gene);
    CV_eff = np.sqrt((1+a)*(1+b) - 1);
    CV_input = np.sqrt(b);

    return v_scores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b


def get_knn_graph2(X, k=5, dist_metric='euclidean',approx=False):
    '''builds knn graph for spring plot'''
    t00 = time.time()
    if approx:
        try:
            from annoy import AnnoyIndex
        except:
            approx = False
            print ('Could not find library "annoy" for approx. nearest neighbor search.')
            print ('Using sklearn instead.')
    if approx:
        print ('Using approximate nearest neighbor search.')
        # Approximate KNN using Annoy
        if dist_metric == 'cosine':
            dist_metric = 'angular'
        npc = X.shape[1]
        ncell = X.shape[0]
        annoy_index = AnnoyIndex(npc, metric=dist_metric)
        t0 = time.time()
        for i in range(ncell):
            annoy_index.add_item(i, list(X[i,:]))
        annoy_index.build(10) # 10 trees
        t1 = time.time() - t0
        #print 'Annoy: index built in %.5f sec' %t1

        t0 = time.time()
        knn = []
        for iCell in range(ncell):
            knn.append(annoy_index.get_nns_by_item(iCell, k + 1)[1:])
        knn = np.array(knn, dtype=int)
        t1 = time.time() - t0
        #print 'kNN built in %.5f sec' %(t1)
    else:
        t0 = time.time()
        if dist_metric == 'cosine':
            nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric, algorithm='brute').fit(X)
        else:
            nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric).fit(X)
        knn = nbrs.kneighbors(return_distance=False)
        t1 = time.time() - t0
        print('kNN built in %.5f sec' %(t1))

    links = set([])
    for i in range(knn.shape[0]):
        for j in knn[i,:]:
            links.add(tuple(sorted((i,j))))

    t11 = time.time() - t00
    print(('kNN built in %.5f sec' %(t11)))
    return list(links), knn


def make_spring_plot(E,
                     gene_list,
                     normalize = True,
                     exclude_dominant_frac = 1.0,
                     min_counts = 3,
                     min_cells = 3,
                     min_vscore_pctl = 90,
                     num_pc = 20,
                     precomputed_pca = [],
                     pca_norm = True,
                     pca_method = '',
                     k_neigh = 4,
                     num_force_iter = 1000,
                     dist_metric = 'euclidean',
                     use_approxnn = False,
                     use_seed=42):

    '''this function outputs spring coordinates and edges'''

    E = E.tocsc()

    if normalize:
        print ('Normalizing')
        E = tot_counts_norm_sparse(E, exclude_dominant_frac = exclude_dominant_frac)[0]
    
    if len(precomputed_pca) == 0:

        # Get gene stats (above Poisson noise, i.e. V-scores)
        print ('Filtering genes')
        Vscores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b = get_vscores_sparse(E)
    
        ix2 = Vscores>0
        Vscores = Vscores[ix2]
        gene_ix = gene_ix[ix2]
        mu_gene = mu_gene[ix2]
        FF_gene = FF_gene[ix2]
    
        # Filter genes: minimum V-score percentile and at least min_counts in at least min_cells
        min_log_vscore = np.percentile(np.log(Vscores), min_vscore_pctl)
    
        ix = (((E[:,gene_ix] >= min_counts).sum(0).A.squeeze() >= min_cells) & (np.log(Vscores) >= min_log_vscore))
        gene_filter = gene_ix[ix]
        n_filtered_genes = len(gene_filter)
        print(('Using %i genes' %(n_filtered_genes)))
    
        # RUN PCA
        # if method == 'sparse': normalizes by stdev
        # if method == anything else: z-score normalizes
        print ('Running PCA')
        Epca = get_PCA_sparseInput(E[:,gene_filter], numpc=num_pc,
                                   method=pca_method, normalize=pca_norm,
                                   seed=use_seed)
     
    else:
        print('Using a provided PCA')
        n_filtered_genes = 'None: PCA pre-computed'
        Epca = precomputed_pca[:, 0:num_pc]

    print ('Building kNN graph')

    links, knn_graph = get_knn_graph2(Epca, k=k_neigh, dist_metric = dist_metric, approx=use_approxnn)

    # Calculate Euclidean distances in the PC space (will be used to build knn graph)
    #print 'Getting distance matrix'
    #D = get_distance_matrix(Epca)
    #D = scipy.spatial.distance.squareform(pdist(Epca, dist_metric))
    # print(links[:10])

    print ('Running ForceAtlas2')

    random.seed(use_seed)

    G = nx.Graph()
    G.add_nodes_from(list(range(Epca.shape[0])))
    G.add_edges_from(list(links))

    forceatlas2 = ForceAtlas2(
                  # Behavior alternatives
                  outboundAttractionDistribution=False,  # Dissuade hubs
                  linLogMode=False,  # NOT IMPLEMENTED
                  adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
                  edgeWeightInfluence=1.0,

                  # Performance
                  jitterTolerance=1.0,  # Tolerance
                  barnesHutOptimize=True,
                  barnesHutTheta=2,
                  multiThreaded=False,  # NOT IMPLEMENTED

                  # Tuning
                  scalingRatio=1.0,
                  strongGravityMode=False,
                  gravity=0.05,
                  # Log
                  verbose=False)

    positions = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=num_force_iter)
    positions = np.array([positions[i] for i in sorted(positions.keys())])
    positions = positions / 5.0
    positions = positions - np.min(positions, axis = 0) - np.ptp(positions, axis = 0) / 2.0
    positions[:,0] = positions[:,0]  + 750
    positions[:,1] = positions[:,1]  + 250

    info_dict = {}
    info_dict['Date'] = '%s' %datetime.now()
    info_dict['Nodes'] = Epca.shape[0]
    info_dict['Filtered_Genes'] = n_filtered_genes
    info_dict['Gene_Var_Pctl'] = min_vscore_pctl
    info_dict['Min_Cells'] = min_cells
    info_dict['Min_Counts'] = min_counts
    info_dict['Num_Neighbors'] = k_neigh
    info_dict['Num_PCs'] = num_pc
    info_dict['Num_Force_Iter'] = num_force_iter

    return positions, links, info_dict
