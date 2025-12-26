import random

import numpy as np
import torch
from SpaGCN.models import simple_GC_DEC
from scipy.sparse import issparse
from sklearn.decomposition import PCA


class SpaGCN(object):
    def __init__(self):
        super(SpaGCN, self).__init__()
        self.l = None

    def set_l(self, l):
        self.l = l

    def train(
        self,
        adata,
        adj,
        num_pcs=50,
        lr=0.005,
        max_epochs=2000,
        weight_decay=0,
        opt="admin",
        init_spa=True,
        init="louvain",  # louvain or kmeans
        n_neighbors=10,  # for louvain
        n_clusters=None,  # for kmeans
        res=0.4,  # for louvain
        tol=1e-3,
    ):
        self.num_pcs = num_pcs
        self.res = res
        self.lr = lr
        self.max_epochs = max_epochs
        self.weight_decay = weight_decay
        self.opt = opt
        self.init_spa = init_spa
        self.init = init
        self.n_neighbors = n_neighbors
        self.n_clusters = n_clusters
        self.res = res
        self.tol = tol
        assert adata.shape[0] == adj.shape[0] == adj.shape[1]
        pca = PCA(n_components=self.num_pcs)
        if issparse(adata.X):
            pca.fit(adata.X.toarray())
            embed = pca.transform(adata.X.toarray())
        else:
            pca.fit(adata.X)
            embed = pca.transform(adata.X)
        ###------------------------------------------###
        if self.l is None:
            raise ValueError("l should not be set before fitting the model!")
        adj_exp = np.exp(-1 * (adj**2) / (2 * (self.l**2)))
        # ----------Train model----------
        self.model = simple_GC_DEC(embed.shape[1], embed.shape[1])
        self.model.fit(
            embed,
            adj_exp,
            lr=self.lr,
            max_epochs=self.max_epochs,
            weight_decay=self.weight_decay,
            opt=self.opt,
            init_spa=self.init_spa,
            init=self.init,
            n_neighbors=self.n_neighbors,
            n_clusters=self.n_clusters,
            res=self.res,
            tol=self.tol,
        )
        self.embed = embed
        self.adj_exp = adj_exp

    def predict(self):
        z, q = self.model.predict(self.embed, self.adj_exp)
        y_pred = torch.argmax(q, dim=1).data.cpu().numpy()
        # Max probability plot
        prob = q.detach().numpy()
        return y_pred, prob


def searchResolution(
    adata,
    adj,
    l,
    target_num,
    start=0.4,
    step=0.1,
    tol=5e-3,
    lr=0.05,
    max_epochs=10,
    r_seed=100,
    t_seed=100,
    n_seed=100,
    max_run=10,
):
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    res = start
    print("Start at res = ", res, "step = ", step)
    clf = SpaGCN()
    clf.set_l(l)
    clf.train(
        adata,
        adj,
        init_spa=True,
        init="louvain",
        res=res,
        tol=tol,
        lr=lr,
        max_epochs=max_epochs,
    )
    y_pred, _ = clf.predict()
    old_num = len(set(y_pred))
    print("Res = ", res, "Num of clusters = ", old_num)
    run = 0
    while old_num != target_num:
        random.seed(r_seed)
        torch.manual_seed(t_seed)
        np.random.seed(n_seed)
        old_sign = 1 if (old_num < target_num) else -1
        clf = SpaGCN()
        clf.set_l(l)
        clf.train(
            adata,
            adj,
            init_spa=True,
            init="louvain",
            res=res + step * old_sign,
            tol=tol,
            lr=lr,
            max_epochs=max_epochs,
        )
        y_pred, _ = clf.predict()
        new_num = len(set(y_pred))
        print("Res = ", res + step * old_sign, "Num of clusters = ", new_num)
        if new_num == target_num:
            res = res + step * old_sign
            print("recommended res = ", str(res))
            return res
        new_sign = 1 if (new_num < target_num) else -1
        if new_sign == old_sign:
            res = res + step * old_sign
            print("Res changed to", res)
            old_num = new_num
        else:
            step = step / 2
            print("Step changed to", step)
        if run > max_run:
            print("Exact resolution not found")
            print("Recommended res = ", str(res))
            return res
        run += 1
    print("recommended res = ", str(res))
    return res
