def get_predictions_kNN(
        self,
        embeddings: np.ndarray,
        k: int = 50,
        ef: int = 100,
        weighting: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    """Get predictions from kNN search results.

    Parameters
    ----------
    embeddings: numpy.ndarray
        Embeddings as a numpy array.
    k: int, default: 50
        The number of nearest neighbors.
    ef: int, default: 100
        The size of the dynamic list for the nearest neighbors.
        See https://github.com/nmslib/hnswlib/blob/master/ALGO_PARAMS.md
    weighting: bool, default: False
        Use distance weighting when getting the consensus prediction.

    Returns
    -------
    predictions: pandas.Series
        A pandas series containing celltype label predictions.
    nn_idxs: numpy.ndarray
        A 2D numpy array of nearest neighbor indices [num_cells x k].
    nn_dists: numpy.ndarray
        A 2D numpy array of nearest neighbor distances [num_cells x k].
    stats: pandas.DataFrame
        Prediction statistics columns:
        "hits" is a json string with the count for every class in k cells.
        "min_dist" is the minimum distance.
        "max_dist" is the maximum distance
        "vs2nd" is sum of best / (best + 2nd best).
        "vsAll" is sum of best / (all hits).
        "hits_weighted" is a json string with the weighted count for every class in k cells.
        "vs2nd_weighted" is weighted sum of best / (best + 2nd best).
        "vsAll_weighted" is weighted sum of best / (all hits).

    Examples
    --------
    >>> ca = CellAnnotation(model_path="/opt/data/model")
    >>> embedding = ca.get_embeddings(align_dataset(data, ca.gene_order).X)
    >>> predictions, nn_idxs, nn_dists, nn_stats = ca.get_predictions_kNN(embeddings)
    """

    start_time = time.time()
    nn_idxs, nn_dists = self.get_nearest_neighbors(
        embeddings=embeddings, k=k, ef=ef
    )
    end_time = time.time()
    print(
        f"Get nearest neighbors finished in: {float(end_time - start_time) / 60} min"
    )
    stats = {
        "hits": [],
        "hits_weighted": [],
        "min_dist": [],
        "max_dist": [],
        "vs2nd": [],
        "vsAll": [],
        "vs2nd_weighted": [],
        "vsAll_weighted": [],
        "dist_percelltype": [],
        "all_nn_indices": [],  # Store indices of NNs for each cell
        "all_nn_names": [],
        "min_dist_index": []

    }
    if k == 1:
        predictions = pd.Series(nn_idxs.flatten()).map(self.idx2label)
    else:
        predictions = []
        # (AT) here nn_idxs, nn_dists has dimension of (# of cells x k).
        # For Yang's data # of cells are ~104k and k is 50 (chosen) which means we have 50 closest neighbor information (indexes and distances)
        # That gives us nns and d_nns as array of 50.
        for nns, d_nns in tqdm(zip(nn_idxs, nn_dists), total=nn_idxs.shape[0]):
            # count celltype in nearest neighbors (optionally with distance weights)
            current_nn_indices = []
            current_nn_names = []

            celltype = defaultdict(float)
            celltype_weighted = defaultdict(float)
            celltype_dist = defaultdict(float)
            for neighbor, dist in zip(nns, d_nns):
                current_nn_indices.append(int(neighbor))
                current_nn_names.append(self.idx2label[neighbor])

                celltype[self.idx2label[
                    neighbor]] += 1  # (AT) length of dictionary idx2label is 7913892 which is the # of single cell profile used to train model.
                celltype_weighted[self.idx2label[
                    neighbor]] += 1 / dist  # (AT) IG this is the scimilarity score which is inverse of the distance (I might be wrong)
                celltype_dist[self.idx2label[neighbor]] += dist
            # predict based on consensus max occurrence
            if weighting:
                predictions.append(
                    max(celltype_weighted.items(), key=operator.itemgetter(1))[0]
                )
            else:
                predictions.append(
                    max(celltype.items(), key=operator.itemgetter(1))[0]
                )
            # compute prediction stats
            stats["hits"].append(json.dumps(celltype))
            stats["hits_weighted"].append(json.dumps(celltype_weighted))
            stats["min_dist"].append(np.min(d_nns))
            stats["max_dist"].append(np.max(d_nns))
            stats["dist_percelltype"].append(json.dumps(celltype_dist))
            stats["all_nn_indices"].append(json.dumps(current_nn_indices))
            stats["all_nn_names"].append(json.dumps(current_nn_names))
            stats["min_dist_index"].append(int(nns[np.argmin(d_nns)]))

            hits = sorted(celltype.values(), reverse=True)
            hits_weighted = sorted(celltype_weighted.values(), reverse=True)
            if len(hits) > 1:
                stats["vs2nd"].append(hits[0] / (hits[0] + hits[1]))
                stats["vsAll"].append(hits[0] / sum(hits))
                stats["vs2nd_weighted"].append(
                    hits_weighted[0] / (hits_weighted[0] + hits_weighted[1])
                )
                stats["vsAll_weighted"].append(
                    hits_weighted[0] / sum(hits_weighted)
                )
            else:
                stats["vs2nd"].append(1.0)
                stats["vsAll"].append(1.0)
                stats["vs2nd_weighted"].append(1.0)
                stats["vsAll_weighted"].append(1.0)
    return (
        pd.Series(predictions),
        nn_idxs,
        nn_dists,
        pd.DataFrame(stats),
    )