# A* Atom Mapper


The A* algorithm of [Heinonen et al.](https://pdfs.semanticscholar.org/f73a/6c2e7a7e248a8f3af16bb74653b1f46c2c82.pdf) is an atom mapping algorithm based on the A* search technique. A* calculates the “cost” in terms of graph editing operations that transform the input graph into the the target graph. The algorithm uses the principle of minimal chemical distance (PMCD) to identify optimal mappings, i.e. those which require fewer graph editing operations. The PMCD states that most chemical reactions follow the shortest path for converting reactants into products, i.e. the path involving the smallest possible number of bond transformation.

This implementation is the result of my master's thesis at Goethe University Frankfurt. The thesis can be found [here](https://raw.githubusercontent.com/zotko/a-star-atom-mapper/master/master_thesis.pdf).
