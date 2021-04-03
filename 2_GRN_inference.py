import re
import pickle
import time
import pandas as pd
import numpy as np
import rcca
from scipy.linalg import norm, eigh
from scipy.stats import fisher_exact,ttest_ind
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import MinMaxScaler
from multiprocessing import Pool
import operator
import matplotlib.pyplot as plt
import pandas_profiling
import networkx as nx
import itertools
from collections import Counter
import mygene
from sklearn.cluster import AgglomerativeClustering
from community import community_louvain
import networkx.algorithms.community as nxcom
from collections import defaultdict
from scipy.spatial import distance
import pandas_profiling
import pickle
from pathlib import Path
import argparse



def weighted_kcca(data, reg=0.5, numCC=None, kernelcca=True, ktype='linear', gausigma=1.0, degree=2):
    """Set up and solve the kernel CCA eigenproblem"""
    if kernelcca:
        kernel = [rcca._make_kernel(d, ktype=ktype, gausigma=gausigma, degree=degree) for d in data]
    else:
        kernel = [d.T for d in data]
        
    nDs = len(kernel)
    nFs = [k.shape[0] for k in kernel]
    numCC = min([k.shape[1] for k in kernel]) if numCC is None else numCC

    # Get the auto- and cross-covariance matrices
    crosscovs = [np.dot(ki, kj.T) for ki in kernel for kj in kernel]

    # Allocate left-hand side (LH) and right-hand side (RH):
    LH = np.zeros((sum(nFs), sum(nFs)))
    RH = np.zeros((sum(nFs), sum(nFs)))

    # Fill the left and right sides of the eigenvalue problem
    for i in range(nDs):
        RH[sum(nFs[:i]) : sum(nFs[:i+1]),
           sum(nFs[:i]) : sum(nFs[:i+1])] = (crosscovs[i * (nDs + 1)]
                                             + reg * np.eye(nFs[i]))

        for j in range(nDs):
            if i != j:
                LH[sum(nFs[:j]) : sum(nFs[:j+1]),
                   sum(nFs[:i]) : sum(nFs[:i+1])] = crosscovs[nDs * j + i]

    LH = (LH + LH.T) / 2.
    RH = (RH + RH.T) / 2.

    maxCC = LH.shape[0]
    r, Vs = eigh(LH, RH, eigvals=(maxCC - numCC, maxCC - 1))
    r[np.isnan(r)] = 0
    rindex = np.argsort(r)[::-1]
    comp = []
    Vs = Vs[:, rindex]
    for i in range(nDs):
        comp.append(Vs[sum(nFs[:i]):sum(nFs[:i + 1]), :numCC])
    return np.sort(r)[::-1], comp
    
class _CCABase(rcca._CCABase):
    def __init__(self, numCV=None, reg=None, regs=None, numCC=None,
                 numCCs=None, kernelcca=True, ktype=None, verbose=False,
                 select=0.2, cutoff=1e-15, gausigma=1.0, degree=2):
        self.numCV = numCV
        self.reg = reg
        self.regs = regs
        self.numCC = numCC
        self.numCCs = numCCs
        self.kernelcca = kernelcca
        self.ktype = ktype
        self.cutoff = cutoff
        self.select = select
        self.gausigma = gausigma
        self.degree = degree
        if self.kernelcca and self.ktype == None:
            self.ktype = 'linear'
        self.verbose = verbose

    def train(self, data):
        nT = data[0].shape[0]
        if self.verbose:
            print('Training CCA, kernel = %s, regularization = %0.4f, '
                  '%d components' % (self.ktype, self.reg, self.numCC))

        eigenV, comps = weighted_kcca(data, self.reg, self.numCC, kernelcca=self.kernelcca, ktype=self.ktype, gausigma=self.gausigma, degree=self.degree)
        ###weighted sample weight, weighted by eigenvalue
        comps = [np.array([eigenV]*nT)*comps[i] for i in range(len(data))]
        
        self.cancorrs, self.ws, self.comps = rcca.recon(data, comps, kernelcca=self.kernelcca)
        if len(data) == 2:
            self.cancorrs = self.cancorrs[np.nonzero(self.cancorrs)]
        return self
    
class weighted_CCA(_CCABase):
    """Attributes:
        reg (float): regularization parameter. Default is 0.1.
        numCC (int): number of canonical dimensions to keep. Default is 10.
        kernelcca (bool): kernel or non-kernel CCA. Default is True.
        ktype (string): type of kernel used if kernelcca is True.
                        Value can be 'linear' (default) or 'gaussian'.
        verbose (bool): default is True.
    Returns:
        ws (list): canonical weights
        comps (list): canonical components
        cancorrs (list): correlations of the canonical components
                         on the training dataset
        corrs (list): correlations on the validation dataset
        preds (list): predictions on the validation dataset
        ev (list): explained variance for each canonical dimension
    """
    def __init__(self, reg=0.5, numCC=10, kernelcca=True, ktype=None,
                 verbose=False, cutoff=1e-15):
        super(weighted_CCA, self).__init__(reg=reg, numCC=numCC, kernelcca=kernelcca,
                                  ktype=ktype, verbose=verbose, cutoff=cutoff)

    def train(self, data):
        return super(weighted_CCA, self).train(data)

def kcca_embedding(TF_exp, TG_exp, normalize, n_comp=1, kernel='gaussian', reg = 0.5):
    kcca = weighted_CCA(kernelcca=True, ktype=kernel, reg = reg, numCC = n_comp) 
    kcca.train([TF_exp.to_numpy(), TG_exp.to_numpy()])
    
    #print('x_weight_dim', kcca.ws[0].shape)
    #print('y_weight_dim', kcca.ws[1].shape)
    
    TF_embed = pd.DataFrame(kcca.ws[0], index=TF_exp.columns)
    TG_embed = pd.DataFrame(kcca.ws[1], index=TG_exp.columns)
    
    if normalize == True:
        TF_embed_norm = TF_embed.apply(lambda x:x/norm(x), axis=1).fillna(0) #L2 norm
        TG_embed_norm = TG_embed.apply(lambda x:x/norm(x), axis=1).fillna(0) #L2 norm
        return (kcca, TF_embed_norm, TG_embed_norm)
    else:
        return (kcca, TF_embed, TG_embed)

def dotProduct(instance1, instance2):
    res = 0
    for x in range(len(instance1)):
        res += instance1[x] * instance2[x]
    return res
    

def getNeighbors(instance):
    '''Generate k-nearest neighbors of testInstance among trainingSet
    =================================================================
    '''
    test, idx1, trainingSet = instance
    testInstance = test.loc[idx1,:]
    li_distances = []
    length = len(testInstance)-1
    for i in range(len(trainingSet)):
        #dist = euclideanDistance(testInstance, trainingSet.iloc[i])
        #dist = correlation(testInstance, trainingSet.iloc[i])
        dist = dotProduct(testInstance, trainingSet.iloc[i])
        li_distances.append((trainingSet.index[i], dist))
    #sort list by the second item in each element from list
    #e.g.) li_distances = [(tg1,dist), (tg2,dist)]
    li_distances.sort(key = operator.itemgetter(1),reverse=True)
    dict_res = dict()
    dict_res[idx1] = li_distances
    #e.g.) dict_res = {tf1:[(tg1,dist),(tg2,dist),(tg3,dist)]}
    return dict_res

def inp_pair(df1,df2):
    li1 = set(df1.index)
    for tf in li1:
        yield (df1,tf,df2)

def TFTG_nwk(df_TF, df_TG):
    outs = {}
    for tf in set(df_TF.index):
        inst = (df_TF,tf,df_TG) 
        outs.update(getNeighbors(inst))
    return outs
    
def inp_pair_modularized_TFTG_nwk(nx_graph, dict_tf_comm, dict_tg_comm, edge_cutoff, df_exp, n_comp, reg, normalize):
    for comID_tf,comID_tg in itertools.product(dict_tf_comm.keys(),dict_tg_comm.keys()):
        yield (nx_graph, dict_tf_comm, dict_tg_comm, comID_tf, comID_tg, edge_cutoff, df_exp, n_comp, reg, normalize)
        
def mergeDicts(iter_dicts):
    dict_all=dict()
    for each_dict in iter_dicts:
        dict_all.update(each_dict)
    return dict_all



def modularized_TFTG_nwk(instance):
    def all_tftg_candidates(nx_obj, tfs, tgs):
        tfs_common, tgs_common = set(tfs).intersection(nx_obj.nodes()), set(tgs).intersection(nx_obj.nodes())
        dict_cent_allnodes=nx.betweenness_centrality_subset(nx_obj,tfs_common,tgs_common)
        set_nodes = set(pd.DataFrame.from_dict(dict_cent_allnodes,orient='index',columns=['paths']).loc[lambda x:x.paths !=0].index.to_list()) ##########
        nx_obj = nx_obj.subgraph(set_nodes.union(tfs).union(tgs)) ###########
        dict_tftgs = nx.to_dict_of_lists(nx_obj)
        dict_tftgs_filtered ={}
        for key in dict_tftgs:
            if len(dict_tftgs[key])==0:
                continue
            else:
                dict_tftgs_filtered[key] = dict_tftgs[key]
        return nx_obj, dict_tftgs_filtered

    def addEdges(dict_pair):
        li_res = []
        for tf in dict_pair.keys():
            for tg,dist in dict_pair[tf]:
                li_res.append((tf,tg,dist))
        return li_res

    def filterEdges(dict_pair,cutoff):
        dict_res = {}
        for tf in dict_pair.keys():
            dict_res[tf] = [(dist) for tg,dist in dict_pair[tf] if dist > cutoff]
        return dict_res

    def list2dict(li_pair):
        dict_res = {}
        for tf,tg in li_pair:
            if tf in dict_res:
                dict_res[tf].append(tg)
            else:
                dict_res[tf]=[tg]
        return dict_res

    def rearrange_dict(dict1):
        dict2 = {}
        for i1 in dict1:
            for i2,i3 in dict1[i1]:
                dict2[(i1,i2)]=i3
        return dict2
    nx_graph, dict_tf_comm, dict_tg_comm, comID_tf, comID_tg, edge_cutoff, df_exp, n_comp, reg, normalize = instance[0], instance[1], instance[2], instance[3], instance[4], instance[5], instance[6], instance[7], instance[8], instance[9]
    res = {}
    # print('TF_cluster: {}, \t TG_cluster:{}'.format(comID_tf,comID_tg))
    nx_graph_tmp, dict_tftg_candidates_tmp = all_tftg_candidates(nx_graph,dict_tf_comm[comID_tf],dict_tg_comm[comID_tg])
    if len(nx_graph_tmp.edges)==0:
        res[(comID_tf,comID_tg)] = []
        return res   
    else:    
        df_graph_tmp = pd.DataFrame(nx_graph_tmp.edges,columns=['TF','TG'])
    ### detecting valid TFTG pair with kCCA
    set_allNodes = set(nx_graph_tmp.nodes())
    set_outDegreeNodes = set(df_graph_tmp.loc[:,'TF'].to_list())
    set_visited = set()
    TFs = set(dict_tf_comm[comID_tf]).intersection(set_outDegreeNodes)
    set_visited.update(TFs)
    TGs = []
    for TF in TFs:
        if TF in dict_tftg_candidates_tmp:
            TGs.extend(dict_tftg_candidates_tmp[TF])
    TGs = set(TGs).intersection(set_allNodes)        
    if len(TGs) == 0:
        res[(comID_tf,comID_tg)] = []
        return res
    kcca, TF_embed, TG_embed = kcca_embedding(df_exp.loc[:,df_exp.columns.isin(TFs)],df_exp.loc[:,df_exp.columns.isin(TGs)],normalize=normalize,n_comp=n_comp,reg=reg)
    df_TFTG_pair_final = pd.DataFrame()
    dict_TFTG_pair = rearrange_dict(TFTG_nwk(TF_embed,TG_embed))
    if len(dict_TFTG_pair)==0:
        res[(comID_tf,comID_tg)] = []
        return res
    else:
        df_TFTG_pair = pd.DataFrame.from_dict(dict_TFTG_pair,orient='index',columns=['weight_kCCA'])
    li_TFTG_pair = list(set(nx_graph_tmp.edges()).intersection(set(df_TFTG_pair.index.to_list()))) ######
    df_TFTG_pair = df_TFTG_pair.loc[li_TFTG_pair]
    df_TFTG_pair_final = df_TFTG_pair
    df_TFTG_pair_2nd = df_TFTG_pair
    while True:
        TFs_2nd = set([tg for tf,tg in df_TFTG_pair_2nd.index.to_list()]).intersection(set_outDegreeNodes)
        set_visited.update(TFs_2nd)
        TGs_2nd = set(sum([dict_tftg_candidates_tmp[TF] for TF in TFs_2nd if TF in dict_tftg_candidates_tmp.keys()],[]))-set_visited
        if len(TGs_2nd) == 0:
            break
        else:
            kcca_2nd, TF_embed_2nd, TG_embed_2nd = kcca_embedding(df_exp.loc[:,df_exp.columns.isin(TFs_2nd)],df_exp.loc[:,df_exp.columns.isin(TGs_2nd)],normalize=normalize,n_comp=n_comp,reg=reg)
            dict_TFTG_pair_2nd = rearrange_dict(TFTG_nwk(TF_embed_2nd,TG_embed_2nd))
            if len(dict_TFTG_pair_2nd)==0:
                break 
            else:
                df_TFTG_pair_2nd = pd.DataFrame.from_dict(dict_TFTG_pair_2nd,orient='index',columns=['weight_kCCA'])
            li_TFTG_pair_2nd = list(set(nx_graph_tmp.edges()).intersection(set(df_TFTG_pair_2nd.index.to_list()))) 
            if len(li_TFTG_pair_2nd) == 0:
                break
            df_TFTG_pair_2nd = df_TFTG_pair_2nd.loc[li_TFTG_pair_2nd]
            df_TFTG_pair_final = pd.concat([df_TFTG_pair_final,df_TFTG_pair_2nd],axis=0)
            li_TFTG_pair += li_TFTG_pair_2nd
    if df_TFTG_pair_final.shape[0]!=0:
        scaler=MinMaxScaler()
        df_TFTG_pair_final['weight_kCCA_rescaled'] = scaler.fit_transform(df_TFTG_pair_final['weight_kCCA'].values.reshape(-1,1))
        #######edge_cutoff with rescaled edge weight 
        li_TFTG_pair_final = df_TFTG_pair_final.loc[lambda x:x.weight_kCCA_rescaled > edge_cutoff].index.to_list() 
        res[(comID_tf,comID_tg)] = li_TFTG_pair_final
        return res
    else:
        res[(comID_tf,comID_tg)] = []
        return res 

def all_tftg_candidates(instance):
    nx_obj, dict_tf_comm, dict_tg_comm, tf_id, tg_id = instance[0], instance[1], instance[2], instance[3], instance[4]
    tfs, tgs = dict_tf_comm[tf_id], dict_tg_comm[tg_id]
    tfs_common, tgs_common = set(tfs).intersection(nx_obj.nodes()), set(tgs).intersection(nx_obj.nodes())
    dict_cent_allnodes=nx.betweenness_centrality_subset(nx_obj,tfs_common,tgs_common)
    set_nodes = set(pd.DataFrame.from_dict(dict_cent_allnodes,orient='index',columns=['paths']).loc[lambda x:x.paths !=0].index.to_list())
    nx_obj = nx_obj.subgraph(set_nodes.union(tfs).union(tgs)) ###########
    dict_res = {}
    li_li_edges = pd.DataFrame(nx_obj.edges).values.tolist()
    set_tup_edges = [(i1,i2) for i1,i2 in li_li_edges]
    dict_res[(tf_id,tg_id)] = set_tup_edges
    return dict_res

def inp_pair_all_tftg_candidates(nx_GRN,dict_cluster_1,dict_cluster_2):
    for tf_id,tg_id in itertools.product(dict_cluster_1.keys(),dict_cluster_2.keys()):
        yield (nx_GRN,dict_cluster_1,dict_cluster_2,tf_id,tg_id)
        
def edgelist2nodeset(dict_edges):
    dict_nodeSet = {}
    for key in dict_edges:
        dict_nodeSet[key] = set(np.array(dict_edges[key]).flatten())
    return dict_nodeSet

def fisher_exact_test(dict_nodeSet, set_ref, nodeSet):
    dict_jaccard = {}
    set_ref = set_ref.intersection(nodeSet)
    for key in dict_nodeSet:
        if len(set(dict_nodeSet[key])) ==0:
            continue
        numDEGs = len(set(dict_nodeSet[key]).intersection(set_ref))
        coeff, pval= fisher_exact([[numDEGs,len(set_ref)-numDEGs],[len(dict_nodeSet[key])-numDEGs, len(nodeSet)-len(set_ref)+numDEGs]])
        dict_jaccard[key] = pval
    return dict_jaccard    


def GRN_inference(nx_GRN, df_exp, li_deg, tfModule, tgModule, nComp, thr, reg, normalize, nThreads=20):
    df_exp_filt=df_exp.loc[:,set(df_exp.columns.to_list()).intersection(set(li_deg))]
    
    set_nodes_tmp = set()
    for node in df_exp_filt.columns.to_list():
        if node in nx_GRN.nodes:
            set_nodes_tmp = set_nodes_tmp.union(node)
            set_nodes_tmp = set_nodes_tmp.union(set(nx_GRN.neighbors(node)))
    numDEGs = len(set_nodes_tmp.intersection(set(li_deg)))
    numTotal = len(set_nodes_tmp)
    nx_GRN_tmp = nx_GRN.subgraph(set_nodes_tmp).copy()
    tf_filt_keys = [key for key in dict_tf_community if len(set(nx_GRN_tmp.nodes()).intersection(set(dict_tf_community[key])))>1]
    tg_filt_keys = [key for key in dict_tg_community if len(set(nx_GRN_tmp.nodes()).intersection(set(dict_tg_community[key])))>1]
    dict_tf_community_filt = {}
    dict_tg_community_filt = {}
    for key in tf_filt_keys:
        dict_tf_community_filt[key] = dict_tf_community[key]
    for key in tg_filt_keys:
        dict_tg_community_filt[key] = dict_tg_community[key]
    inps = inp_pair_modularized_TFTG_nwk(nx_GRN_tmp, dict_tf_community_filt, dict_tg_community_filt, thr, df_exp, nComp, reg, normalize)
    with Pool(nThreads) as p:
        outs = p.map(modularized_TFTG_nwk,inps)
    return mergeDicts(outs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='python 2_GRN_inference.py nwk_GRN exp deg out -TFli [TF cataolog] -TFmodule [] -TGmodule [] -nComp [# components] -thr [thr] -nThreads [] -normalize []')
    parser.add_argument('nwk_GRN',help="GRN to be utilized as a guide TF-TG relation")
    parser.add_argument('exp',help='exp profile of genes')
    parser.add_argument('deg',help='deg list of genes')
    parser.add_argument('out',help='out file name')
    parser.add_argument('-TFli',required=True, help="TF catalog")
    parser.add_argument('-TFmodule',required=True)
    parser.add_argument('-TGmodule',required=True)
    parser.add_argument('-nComp',required=True, type=int, help='# components')
    parser.add_argument('-thr', required=True, type=float, help='')
    parser.add_argument('-reg',type=float)
    parser.add_argument('-nThreads',required=False, type=int)
    parser.add_argument('-normalize',required=False, default=True, type=bool)
    args = parser.parse_args()
    
    df_tfs = pd.read_csv(args.TFli,sep='\t')
    li_TFs = df_tfs.iloc[:,0].to_list()
    with open(args.TFmodule,'rb') as f:
        dict_tf_community =pickle.load(f)

    with open(args.TGmodule,'rb') as f:
        dict_tg_community =pickle.load(f)
    
    GRN = pd.read_csv(args.nwk_GRN,sep='\t')
    GRN = nx.from_pandas_edgelist(GRN,'TF','TG',create_using=nx.DiGraph())
   
    exp = pd.read_csv(args.exp,sep='\t',index_col=0 ).T 
    
    deg = pd.read_csv(args.deg,sep='\t',names=['gene']).iloc[:,0].to_list()
    dict_res = GRN_inference(GRN,exp,deg,dict_tf_community,dict_tg_community,args.nComp,args.thr,args.reg,args.normalize,args.nThreads)
    
    with open(args.out,'wb') as f:
        pickle.dump(dict_res, f)        
