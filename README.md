# Data Science from Scratch
This section shows a snapshot of code in the book [Data Science from Scratch](https://github.com/joelgrus/data-science-from-scratch).

## A/B test
* Idea:
    * draw graph
    * build `normal_cdf(x, mu=0,sigma=1)`, `inverse_normal_cdf(p, mu=0, sigma=1, tolerance=0.00001)` (by *binary search*)
    * based on type of alternative hypothesis, discuss `normal_cdf`, `inverse_normal_cdf`
    * plug in number based on graphs
* Input: sample data, true value
* Output: p-value, power
* Representation: `R pnorm()`, `R qnorm()`
* Implementation:
```
def inverse_normal_cdf(p, mu=0, sigma=1, tolerance=0.00001):
    """find approximate inverse using binary search"""

    # if not standard, compute standard and rescale
    if mu != 0 or sigma != 1:
        return mu + sigma * inverse_normal_cdf(p, tolerance=tolerance)

    low_z, low_p = -10.0, 0            # normal_cdf(-10) is (very close to) 0
    hi_z,  hi_p  =  10.0, 1            # normal_cdf(10)  is (very close to) 1
    while hi_z - low_z > tolerance:
        mid_z = (low_z + hi_z) / 2     # consider the midpoint
        mid_p = normal_cdf(mid_z)      # and the cdf's value there
        if mid_p < p:
            # midpoint is still too low, search above it
            low_z, low_p = mid_z, mid_p
        elif mid_p > p:
            # midpoint is still too high, search below it
            hi_z, hi_p = mid_z, mid_p
        else:
            break

    return mid_z
```

## Gradient Descent
* Idea:  x_{m+1} = x_m - alpha * g(f(x_m)), where g() is gradient function
    1. define `minimize_batch(target_fn, gradient_fn, x_0, tolerance=0.000001)`
    2. need to pick `step_sizes` - alpha
    3. `gradient_fn(x)` returns `value`
    4. choose the `value` that minimizes the error function and see if it converges
* Input: `target_fn`, `x_0`
* Output `x` that minimizes `target_fn`
* Implementation:
```
def minimize_batch(target_fn, gradient_fn, theta_0, tolerance=0.000001):
    """use gradient descent to find theta that minimizes target function"""

    step_sizes = [100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001]

    theta = theta_0                           # set theta to initial value
    target_fn = safe(target_fn)               # safe version of target_fn
    value = target_fn(theta)                  # value we're minimizing

    while True:
        gradient = gradient_fn(theta)  
        next_thetas = [step(theta, gradient, -step_size)
                       for step_size in step_sizes]

        # choose the one that minimizes the error function        
        next_theta = min(next_thetas, key=target_fn)
        next_value = target_fn(next_theta)

        # stop if we're "converging"
        if abs(value - next_value) < tolerance:
            return theta
        else:
            theta, value = next_theta, next_value
```

## Principal component analysis
* Idea: extract one or more dimensions that capture as much of the variation in the data as possible
    1. def `de_mean_matrix(X)` to de-meaned matrix
    2. Specifically, given a direction `w` (a vector of magnitude 1), each row `x` in the matrix extends `dot(x, d)` in the `d` direction
    3. compute the variance of our data set in the direction determined by `w` by `directional_variance_i(x_i, w)` and `directional_variance(X, w)`
    4. find the direction that maximizes this variance by gradient descent with function `directional_variance_gradient_i(x_i, w)` and `directional_variance_gradient_i(x_i, w)`
* Input: data matrix `X`
* Output: direction vector `w`
* Implementation:
```
def first_principal_component(X):
    guess = [1 for _ in X[0]]
    unscaled_maximizer = maximize_batch(
        partial(directional_variance, X), # is now a function of w
        partial(directional_variance_gradient, X), # is now a function of w
        guess)
    return direction(unscaled_maximizer)

def project(v, w):
    """return the projection of v onto the direction w"""
    projection_length = dot(v, w)
    return scalar_multiply(projection_length, w)

def remove_projection_from_vector(v, w):
    """projects v onto w and subtracts the result from v"""
    return vector_subtract(v, project(v, w))
```

## KNN
* Idea: an observation is classified by a majority vote of its `k` nearest neighbors
* Input: `labeled_points`: `(point, label)`
* Output: predicted label
* Implementation:
    * def `majority_vote(labels)` that can handle with ties
    * Put in together:
    ```
    def knn_classify(k, labeled_points, new_point):
        # order the labeled points from nearest to farthest
        by_distance = sorted(labeled_points,
                             key=lambda (point, _): distance(point, new_point))

        # find the labels for the k closest
        k_nearest_labels = [label for _, label in by_distance[:k]]

        # and let them vote
        return majority_vote(k_nearest_labels)
    ```    

## Multiple Regression with Regularization
* Idea: estimate parameters by minimizing l2-loss function with penalty term using gradient descent
* Input: `(X, y)`
* Output: `(X_new, y_new)`
* Implementation:
    * def `ridge_penalty(beta, alpha)`
    * def `squared_error_ridge(x_i, y_i, beta, alpha)`
    * def `ridge_penalty_gradient(beta, alpha)`
    * def `squared_error_ridge_gradient(x_i, y_i, beta, alpha)`
    * def `estimate_beta_ridge(x, y, alpha)`

## Logistic Regression
* Idea: estimate parameters by minimizing log-loss function using gradient descent
* Input: `(X, y)`
* Output: `(X_new, y_new)`
* Implementation:
    * def `logistic_log_partial_ij(x_i, y_i, beta, j)`
    * def `def logistic_log_gradient_i(x_i, y_i, beta)`
    * def `logistic_log_gradient(x, y, beta)`

##  Naive Bayes
* Idea: p(spam|word_i, ..., word_n) = p(word_1|spam) * ... * p(word_n|spam) * p(spam)/p(word_i, ..., word_n)
    1. tokenize messages into distinct words by `tokenize(message)`
    2. def `count_words(training_set) """training set consists of pairs (message, is_spam)"""`
    3. def
    ```
    def word_probabilities(counts, total_spams, total_non_spams, k=0.5):
        """turn the word_counts into a list of triplets w, p(w | spam) and p(w | ~spam)"""
    ```
    4. def `spam_probability(word_probs, message)`
* Input: `(message, is_spam)`
* Output: `(new_message, is_spam)`
* Implementation:
```
def count_words(training_set):
    """training set consists of pairs (message, is_spam)"""
    counts = defaultdict(lambda: [0, 0])
    for message, is_spam in training_set:
        for word in tokenize(message):
            counts[word][0 if is_spam else 1] += 1
    return counts

def word_probabilities(counts, total_spams, total_non_spams, k=0.5):
    """turn the word_counts into a list of triplets
    w, p(w | spam) and p(w | ~spam)"""
    return [(w,
             (spam + k) / (total_spams + 2 * k),
             (non_spam + k) / (total_non_spams + 2 * k))
             for w, (spam, non_spam) in counts.iteritems()]

def spam_probability(word_probs, message):
    message_words = tokenize(message)
    log_prob_if_spam = log_prob_if_not_spam = 0.0

    for word, prob_if_spam, prob_if_not_spam in word_probs:

        # for each word in the message,
        # add the log probability of seeing it
        if word in message_words:
            log_prob_if_spam += math.log(prob_if_spam)
            log_prob_if_not_spam += math.log(prob_if_not_spam)

        # for each word that's not in the message
        # add the log probability of _not_ seeing it
        else:
            log_prob_if_spam += math.log(1.0 - prob_if_spam)
            log_prob_if_not_spam += math.log(1.0 - prob_if_not_spam)

    prob_if_spam = math.exp(log_prob_if_spam)
    prob_if_not_spam = math.exp(log_prob_if_not_spam)
    return prob_if_spam / (prob_if_spam + prob_if_not_spam)


class NaiveBayesClassifier:

    def __init__(self, k=0.5):
        self.k = k
        self.word_probs = []

    def train(self, training_set):

        # count spam and non-spam messages
        num_spams = len([is_spam
                         for message, is_spam in training_set
                         if is_spam])
        num_non_spams = len(training_set) - num_spams

        # run training data through our "pipeline"
        word_counts = count_words(training_set)
        self.word_probs = word_probabilities(word_counts,
                                             num_spams,
                                             num_non_spams,
                                             self.k)

    def classify(self, message):
        return spam_probability(self.word_probs, message)
```

## Decision Trees
* Idea:
    1. If the data all have the same label, then create a leaf node that predicts that label and then stop.
    2. If the list of attributes is empty (i.e., there are no more possible questions to ask), then create a leaf node that predicts the most common label and then stop.
    3. Otherwise, try partitioning the data by each of the attributes
    4. Choose the partition with the lowest partition entropy
    5. Add a decision node based on the chosen attribute
    6. Recur on each partitioned subset using the remaining attributes
* Input:
```
inputs = [
    ({'level':'Senior', 'lang':'Java', 'tweets':'no', 'phd':'no'}, False), ({'level':'Senior','lang':'Java','tweets':'no','phd':'yes'}, False), ({'level':'Mid', 'lang':'Python', 'tweets':'no', 'phd':'no'}, True), ({'level':'Junior','lang':'Python','tweets':'no','phd':'no'}, True)
]
```
* Output: `assignments`
* Implementation:
```
def group_by(items, key_fn):
    """returns a defaultdict(list), where each input item
    is in the list whose key is key_fn(item)"""
    groups = defaultdict(list)
    for item in items:
        key = key_fn(item)
        groups[key].append(item)
    return groups

def partition_by(inputs, attribute):
    """returns a dict of inputs partitioned by the attribute
    each input is a pair (attribute_dict, label)"""
    return group_by(inputs, lambda x: x[0][attribute])    

def partition_entropy_by(inputs,attribute):
    """computes the entropy corresponding to the given partition"""        
    partitions = partition_by(inputs, attribute)
    return partition_entropy(partitions.values())        

def classify(tree, input):
    """classify the input using the given decision tree"""

    # if this is a leaf node, return its value
    if tree in [True, False]:
        return tree

    # otherwise find the correct subtree
    attribute, subtree_dict = tree

    subtree_key = input.get(attribute)  # None if input is missing attribute

    if subtree_key not in subtree_dict: # if no subtree for key,
        subtree_key = None              # we'll use the None subtree

    subtree = subtree_dict[subtree_key] # choose the appropriate subtree
    return classify(subtree, input)     # and use it to classify the input

def build_tree_id3(inputs, split_candidates=None):

    # if this is our first pass,
    # all keys of the first input are split candidates
    if split_candidates is None:
        split_candidates = inputs[0][0].keys()

    # count Trues and Falses in the inputs
    num_inputs = len(inputs)
    num_trues = len([label for item, label in inputs if label])
    num_falses = num_inputs - num_trues

    if num_trues == 0:                  # if only Falses are left
        return False                    # return a "False" leaf

    if num_falses == 0:                 # if only Trues are left
        return True                     # return a "True" leaf

    if not split_candidates:            # if no split candidates left
        return num_trues >= num_falses  # return the majority leaf

    # otherwise, split on the best attribute
    best_attribute = min(split_candidates,
        key=partial(partition_entropy_by, inputs))

    partitions = partition_by(inputs, best_attribute)
    new_candidates = [a for a in split_candidates
                      if a != best_attribute]

    # recursively build the subtrees
    subtrees = { attribute : build_tree_id3(subset, new_candidates)
                 for attribute, subset in partitions.iteritems() }

    subtrees[None] = num_trues > num_falses # default case

    return (best_attribute, subtrees)

def forest_classify(trees, input):
    votes = [classify(tree, input) for tree in trees]
    vote_counts = Counter(votes)
    return vote_counts.most_common(1)[0][0]
```    

## Neural Networks
* Idea: a neural network for XOR (see Figure 18-3)
* Input: `neural_network`, `input_vector = (x, y)`
* Output: probability
* Representation: a neural network as a list (layers) of lists (neurons) of lists (weights)
* Implementation:
```
def neuron_output(weights, inputs):
    return sigmoid(dot(weights, inputs))

def feed_forward(neural_network, input_vector):
    """takes in a neural network (represented as a list of lists of lists of weights)
    and returns the output from forward-propagating the input"""

    outputs = []

    for layer in neural_network:

        input_with_bias = input_vector + [1]             # add a bias input
        output = [neuron_output(neuron, input_with_bias) # compute the output
                  for neuron in layer]                   # for this layer
        outputs.append(output)                           # and remember it

        # the input to the next layer is the output of this one
        input_vector = output

    return

# If neural_network is unknown, estimate it by minimizing errors with backpropagation
def backpropagate(network, input_vector, target):

    hidden_outputs, outputs = feed_forward(network, input_vector)

    # the output * (1 - output) is from the derivative of sigmoid
    output_deltas = [output * (1 - output) * (output - target[i])
                     for i, output in enumerate(outputs)]

    # adjust weights for output layer (network[-1])
    for i, output_neuron in enumerate(network[-1]):
        for j, hidden_output in enumerate(hidden_outputs + [1]):
            output_neuron[j] -= output_deltas[i] * hidden_output

    # back-propagate errors to hidden layer
    hidden_deltas = [hidden_output * (1 - hidden_output) *
                      dot(output_deltas, [n[i] for n in network[-1]])
                     for i, hidden_output in enumerate(hidden_outputs)]

    # adjust weights for hidden layer (network[0])
    for i, hidden_neuron in enumerate(network[0]):
        for j, input in enumerate(input_vector + [1]):
            hidden_neuron[j] -= hidden_deltas[i] * input
```

## KMeans
* Idea:
    1. Start with a set of k-means
    2. Assign each point to the mean to which it is closest
    3. If no point’s assignment has changed, stop and keep the clusters. Otherwise, go to previous step    
* Input: `(inputs, k)`
* Output: `assignments`
* Representation: define tree to be one of `True`, `False`, a tuple (`attribute`, `subtree_dict`)
* Implementation:
```
class KMeans:
    """performs k-means clustering"""

    def __init__(self, k):
        self.k = k          # number of clusters
        self.means = None   # means of clusters

    def classify(self, input):
        """return the index of the cluster closest to the input"""
        return min(range(self.k),
                   key=lambda i: squared_distance(input, self.means[i]))

    def train(self, inputs):

        self.means = random.sample(inputs, self.k)
        assignments = None

        while True:
            # Find new assignments
            new_assignments = map(self.classify, inputs)

            # If no assignments have changed, we're done.
            if assignments == new_assignments:                
                return

            # Otherwise keep the new assignments,
            assignments = new_assignments    

            for i in range(self.k):
                i_points = [p for p, a in zip(inputs, assignments) if a == i]
                # avoid divide-by-zero if i_points is empty
                if i_points:                                
                    self.means[i] = vector_mean(i_points)   
```

## Bottom-up Hierarchical Clustering
* Idea:  
    1. Make each input its own cluster of one.
    2. As long as there are multiple clusters remaining, find the two closest clusters and merge them.
* Representation:
    * `leaf1 = ([10, 20],) # to make a 1-tuple you need the trailing comma`
    * `merged = (1, [leaf1, leaf2]) # we will represent as 2-tuples (merge order, children)`
* Implementation:
```
def bottom_up_cluster(inputs, distance_agg=min):
    # start with every input a leaf cluster / 1-tuple
    clusters = [(input,) for input in inputs]

    # as long as we have more than one cluster left...
    while len(clusters) > 1:
        # find the two closest clusters
        c1, c2 = min([(cluster1, cluster2)
                     for i, cluster1 in enumerate(clusters)
                     for cluster2 in clusters[:i]],
                     key=lambda (x, y): cluster_distance(x, y, distance_agg))

        # remove them from the list of clusters
        clusters = [c for c in clusters if c != c1 and c != c2]

        # merge them, using merge_order = # of clusters left
        merged_cluster = (len(clusters), [c1, c2])

        # and add their merge
        clusters.append(merged_cluster)

    # when there's only one cluster left, return it
    return clusters[0]
```

## Recommender Systems (UserCF)
* Idea:  `Customers like you also bought ...`
    1. Create function `user_interest_matrix` by using `make_user_interest_vector(user_interests)`
    2. Create `user_similarities` matrix by some `cosine_similarity(w, v)` function
    3. Find most K similar users by `most_similar_users_to(user_id)`
    4. Make recommendation by `user_based_suggestions(user_id, include_current_interests=False)`
* Input:
```
users_interests = [
    ["Hadoop", "Big Data", "HBase", "Java", "Spark", "Storm", "Cassandra"],
    ["NoSQL", "MongoDB", "Cassandra", "HBase", "Postgres"]
]    
```
* Output:
```
[('MapReduce', 0.5669467095138409),
('MongoDB', 0.50709255283711),
('Postgres', 0.50709255283711),
('NoSQL', 0.3380617018914066),
('neural networks', 0.1889822365046136),
('deep learning', 0.1889822365046136), ('artificial intelligence', 0.1889822365046136),
#...
]
```
* Implementation:
```
#
# user-based filtering
#

def cosine_similarity(v, w):
    return dot(v, w) / math.sqrt(dot(v, v) * dot(w, w))

unique_interests = sorted(list({ interest
                                 for user_interests in users_interests
                                 for interest in user_interests }))

def make_user_interest_vector(user_interests):
    """given a list of interests, produce a vector whose i-th element is 1
    if unique_interests[i] is in the list, 0 otherwise"""
    return [1 if interest in user_interests else 0
            for interest in unique_interests]

user_interest_matrix = map(make_user_interest_vector, users_interests)

user_similarities = [[cosine_similarity(interest_vector_i, interest_vector_j)
                      for interest_vector_j in user_interest_matrix]
                     for interest_vector_i in user_interest_matrix]

def most_similar_users_to(user_id):
    pairs = [(other_user_id, similarity)                      # find other
             for other_user_id, similarity in                 # users with
                enumerate(user_similarities[user_id])         # nonzero
             if user_id != other_user_id and similarity > 0]  # similarity

    return sorted(pairs,                                      # sort them
                  key=lambda (_, similarity): similarity,     # most similar
                  reverse=True)                               # first


def user_based_suggestions(user_id, include_current_interests=False):
    # sum up the similarities
    suggestions = defaultdict(float)
    for other_user_id, similarity in most_similar_users_to(user_id):
        for interest in users_interests[other_user_id]:
            suggestions[interest] += similarity

    # convert them to a sorted list
    suggestions = sorted(suggestions.items(),
                         key=lambda (_, weight): weight,
                         reverse=True)

    # and (maybe) exclude already-interests
    if include_current_interests:
        return suggestions
    else:
        return [(suggestion, weight)
                for suggestion, weight in suggestions
                if suggestion not in users_interests[user_id]]

#
# Item-Based Collaborative Filtering
#

interest_user_matrix = [[user_interest_vector[j]
                         for user_interest_vector in user_interest_matrix]
                        for j, _ in enumerate(unique_interests)]

interest_similarities = [[cosine_similarity(user_vector_i, user_vector_j)
                          for user_vector_j in interest_user_matrix]
                         for user_vector_i in interest_user_matrix]

def most_similar_interests_to(interest_id):
    similarities = interest_similarities[interest_id]
    pairs = [(unique_interests[other_interest_id], similarity)
             for other_interest_id, similarity in enumerate(similarities)
             if interest_id != other_interest_id and similarity > 0]
    return sorted(pairs,
                  key=lambda (_, similarity): similarity,
                  reverse=True)

def item_based_suggestions(user_id, include_current_interests=False):
    suggestions = defaultdict(float)
    user_interest_vector = user_interest_matrix[user_id]
    for interest_id, is_interested in enumerate(user_interest_vector):
        if is_interested == 1:
            similar_interests = most_similar_interests_to(interest_id)
            for interest, similarity in similar_interests:
                suggestions[interest] += similarity

    suggestions = sorted(suggestions.items(),
                         key=lambda (_, similarity): similarity,
                         reverse=True)

    if include_current_interests:
        return suggestions
    else:
        return [(suggestion, weight)
                for suggestion, weight in suggestions
                if suggestion not in users_interests[user_id]]
```

## MapReduce
```
def wc_mapper(document):
    """for each word in the document, emit (word,1)

    >>> wc_mapper(document)
    { "data" : [1, 1],
      "science" : [1, 1],
      "big" : [1],
      "fiction" : [1] }
    """        

    for word in tokenize(document):
        yield (word, 1)

def wc_reducer(word, counts):
    """sum up the counts for a word

    >>> wc_reducer(word, counts)
    [("data", 2), ("science", 2), ("big", 1), ("fiction", 1)]
    """
    yield (word, sum(counts))

def word_count(documents):
    """count the words in the input documents using MapReduce"""

    # place to store grouped values
    collector = defaultdict(list)

    for document in documents:
        for word, count in wc_mapper(document):
            collector[word].append(count)

    return [output
            for word, counts in collector.iteritems()
            for output in wc_reducer(word, counts)]

def map_reduce(inputs, mapper, reducer):
    """runs MapReduce on the inputs using mapper and reducer

    >>> word_counts = map_reduce(documents, wc_mapper, wc_reducer)
    """
    collector = defaultdict(list)

    for input in inputs:
        for key, value in mapper(input):
            collector[key].append(value)

    return [output
            for key, values in collector.iteritems()
            for output in reducer(key,values)]
```

# cs 61a
##  Distributed data

* **Map phase**: Apply a `mapper` function to all inputs, emitting intermediate key-value pairs
    * The mapper yields zero or more key-value pairs for each input (word count example: `[('e', 1), ('o', 2)]` for each doc)
* **Reduce phase**: For each intermediate key, apply a `reducer` function to accumulate all
values associated with that key
    * All key-value pairs with the same key are processed together
```
def vowels(line):
    """Yield (vowel, count) pairs."""
    for v in 'aeiou':
        if v in line:
            yield (v, line.count(v))

def count_with_spark():
    """Count vowels in a text file.

    Run with spark-submit.
    """
    from pyspark import SparkContext
    from operator import add
    sc = SparkContext(appName="VowelCount")

    lines = sc.textFile('shakespeare.txt')
    vowels = lines.flatMap(vowels).reduceByKey(add).sortByKey().collect()
    print(vowels)
```

## Self Join
```
-- Siblings
select a.child as first, b.child as second
  from parents as a, parents as b
  where a.parent = b.parent and a.child < b.child;

-- Grandparents
create table grandparents as
  select a.parent as grandog, b.child as granpup
    from parents as a, parents as b
    where b.parent = a.child;

-- Numbers
create table odds as
  with
    odds(n) as (
      select 1 union
      select n+2 from odds where n < 15
    )
  select n from odds;

with
  i(n) as (
    select 1 union select n+1 from i where n < 20
  ) # i: table name; n column name
select a.n, b.n, c.n from i as a, i as b, i as c
       where a.n < b.n and a.n*a.n + b.n*b.n = c.n*c.n;
```
