import numpy as np

''' Phylogeny
Phylogeny class defined by evolution rate model. Note that sequences
are added to the phylogeny individually and NOT at initialization.
'''
class Phylogeny:

    # Constant value.
    BASES = ['A','C','G','T']


    # _______________________________________________________________
    #                                             Object Construction

    ''' __init__
    Initialize a phylogeny from a given evolution rate model.
    '''
    def __init__(self, rates, nodes, tree, root, data, seqlen, branch_lengths):
        # Rates are a dictionary of the form,
        # {<rate class>: {
        #     'rate': <evolution rate>, 
        #     'prob': <rate equilibrium probabiliy>}}
        self.rates = rates

        # List of nodes in tree.
        self.nodes = nodes

        # Tree topology of the form,
        # {<node>: [<children of node>]}
        self.tree = tree

        # Will be some node in tree.
        self.root = root

        # Genentic data of the form,
        # {<sequence ID>: <sequence>}
        self.data = data

        self.seqlen = seqlen

        # Distance along branches between nodes of the form,
        # {<node1>: {
        #     <node2>: <distance from node1 to node2>}}
        self.branch_lengths = branch_lengths


    # _______________________________________________________________
    #                                                Graph Algorithms

    ''' _bfs
    Return the order in which nodes are traversed in a breadth-first 
    from the root node.
    '''
    def _bfs(self):
        # Begin with no nodes visited.
        visited = []
        queue = [self.root]
     
        # While there are still nodes to be visited.
        while queue:
            node = queue.pop(0)

            if node not in visited:
                visited.append(node)
                neighbours = self.tree[node]
                
                for neighbour in neighbours:
                    queue.append(neighbour)

        return visited


    # _______________________________________________________________
    #                                                Model Algorithms
    
    ''' _log_sum
    Calculate log(a + b).
    '''
    def _log_sum(self, a, b):
        # Avoid returning NaN.
        if a == -np.Inf and b == -np.Inf:
            return -np.Inf

        elif a > b: 
            return a + np.log1p(np.exp(b - a))

        else: 
            return b + np.log1p(np.exp(a - b))


    ''' _p_rate_trans
    Get the probability the evolution rate trasitions from i to j.
    '''
    def _p_rate_trans(self, i, j, auto_coef=np.log(0.7)):
        if i == j:
            return self._log_sum(auto_coef, 
                np.log(1 - np.exp(auto_coef)) + self.rates[j]['prob'])
        else:
            return (np.log(1 - np.exp(auto_coef)) + self.rates[j]['prob'])


    ''' _jcm
    Evaluate M_{ij}(time, rate) under the Jukes-Cantor model. In this 
    model, we are assuming that branch length is (time * rate).
    Arguments:
        i, j: DNA bases
        length: branch length
        rate: rate of evolution
    '''
    def _jcm(self, i, j, length, rate):
        if i == j:
            return np.log(0.25) + np.log(1 + 3 * np.exp((-4/3) * length * rate))
        else:
            return np.log(0.25) + np.log(1 - np.exp((-4/3) * length * rate))


    ''' _update_p_nodes
    Compute the node probabilities (ell values in paper) for each site.
    '''
    def _get_p_nodes(self, site, rate):
        p_nodes = {node:{b:-np.Inf for b in self.BASES} for node in self.nodes}

        order = self._bfs()

        for node in reversed(order):
            for basis in self.BASES:
                # If we are at a leaf.
                if not self.tree[node]:
                    # Kronecker delta function.
                    delta = 0 if (self.data[node][site] == basis) else -np.Inf
                    p_nodes[node][basis] = delta

                # If the node has children.
                else:
                    left, right = self.tree[node]
                    l_dist = self.branch_lengths[node][left]
                    r_dist = self.branch_lengths[node][right]

                    for x in self.BASES:
                        for y in self.BASES:
                            l_prob = self._jcm(basis, x, l_dist, self.rates[rate]['rate']) \
                                + p_nodes[left][x]

                            r_prob = self._jcm(basis, y, r_dist, self.rates[rate]['rate']) \
                                + p_nodes[right][y]

                            p_nodes[node][basis] = self._log_sum(p_nodes[node][basis], l_prob + r_prob)

        return p_nodes


    ''' likelihood
    Compute the likelihood of a tree.
    '''                            
    def likelihood(self):
        # Likelihood of the tree given data and that a site has some 
        # specfic rate category. Initialized outisde of loop so we 
        # can make use of recursive formula.
        ll_tree = {rate:-np.Inf for rate in self.rates}

        # Likelihood of the contribution of rates that maximizes 
        # likelihood at a given site.
        ll_rates = {rate:-np.Inf for rate in self.rates}

        site_rates = [{rate:None for rate in self.rates} for site in range(self.seqlen)]

        # Note that we start at the end of the sequence.
        for site in reversed(range(self.seqlen)):
            #print(ll_tree)
            ll_site = {rate:-np.Inf for rate in self.rates}

            # Compute the likelihood at each site.
            for rate in self.rates:
                p_nodes = self._get_p_nodes(site, rate)

                tail = self.root
                head = self.tree[tail][0]
                length = self.branch_lengths[tail][head]

                for x in self.BASES:
                    for y in self.BASES:
                        ll_site[rate] = self._log_sum(ll_site[rate],
                            np.log(0.25)
                            + p_nodes[tail][x] 
                            + p_nodes[head][y]
                            + self._jcm(x, y, length, self.rates[rate]['rate']))

            # Recursively calculate next.
            # Base case.
            if site == (self.seqlen - 1):
                for rate in self.rates:
                    ll_tree[rate] = ll_site[rate]
                    ll_rates[rate] = ll_site[rate]
                    site_rates[site][rate] = rate

            else:
                ll_tree_new = {rate:-np.Inf for rate in self.rates}
                ll_rates_new = {rate:-np.Inf for rate in self.rates}

                for i in self.rates:
                    rate_coef = -np.Inf

                    max_rate_contribution = -np.Inf
                    max_rate_category = i

                    for j in self.rates:
                        ll_contribution = self._p_rate_trans(i, j) + ll_tree[j]
                        rate_contribution = self._p_rate_trans(i, j) + ll_rates[j]

                        rate_coef = self._log_sum(rate_coef, ll_contribution)

                        if rate_contribution > max_rate_contribution:
                            max_rate_contribution = rate_contribution
                            max_rate_category = j

                    ll_tree_new[i] = ll_site[i] + rate_coef
                    ll_rates_new[i] = ll_site[i] + max_rate_contribution
                    site_rates[site][i] = max_rate_category

                ll_tree = ll_tree_new
                ll_rates = ll_rates_new
            # print(ll_rates)

        # Compute final likelihood of tree.
        ll = -np.Inf

        for rate in self.rates:
            ll = self._log_sum(ll, ll_tree[rate])

        # Compute maximal sequence of rates by backtracking through 
        # our maximal choices.
        final_rate_list = []
        for site in range(self.seqlen):
            if site == 0:
                final_rate_list.append(max(site_rates[site], key=site_rates[site].get))
            else:
                final_rate_list.append(site_rates[site][final_rate_list[site-1]])

        return final_rate_list, ll

