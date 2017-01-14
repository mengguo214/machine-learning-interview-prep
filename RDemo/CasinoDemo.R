# Description: 
#   Demo of the occasionally dishonest casino
#
# Details:
#   1. State Space (suppose we have two: fair dice /loaded dice)
#   2. Switch between two states is according to a given matrix (Markov transition matrix)
#   3. Output is probabilistic but depends on the state (fair/loaded)
#   4. Want to guess the hidden state (fair/loaded) from the output observed (here tosses of the die).

# License: GNU v3
#
# Date: October 19, 2016
#
# Authors:
#    Rcode was written by Jiali Lin.
#    Depts. of Statistics, Virginia Tech,
#    Hutcheson Hall, 403K, Blacksburg, VA 24061 
#
# References:
#    http://web.stanford.edu/class/stats366/hmmR2.html
# 
# See Also:
#     https://cran.r-project.org/web/packages/HMM/HMM.pdf

rm(list = ls())

initHMM <-  function (States, Symbols, startProbs = NULL, transProbs = NULL, 
                      emissionProbs = NULL) 
{
  nStates = length(States)
  nSymbols = length(Symbols)
  S = rep(1/nStates, nStates)
  T = 0.5 * diag(nStates) + array(0.5/(nStates), c(nStates, 
                                                   nStates))
  E = array(1/(nSymbols), c(nStates, nSymbols))
  names(S) = States
  dimnames(T) = list(from = States, to = States)
  dimnames(E) = list(states = States, symbols = Symbols)
  if (!is.null(startProbs)) {
    S[] = startProbs[]
  }
  if (!is.null(transProbs)) {
    T[, ] = transProbs[, ]
  }
  if (!is.null(emissionProbs)) {
    E[, ] = emissionProbs[, ]
  }
  return(list(States = States, Symbols = Symbols, startProbs = S, 
              transProbs = T, emissionProbs = E))
}

simHMM <- function (hmm, length) 
{
  hmm$transProbs[is.na(hmm$transProbs)] = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  states = c()
  emission = c()
  states = c(states, sample(hmm$States, 1, prob = hmm$startProbs))
  for (i in 2:length) {
    state = sample(hmm$States, 1, prob = hmm$transProbs[states[i - 
                                                                 1], ])
    states = c(states, state)
  }
  for (i in 1:length) {
    emi = sample(hmm$Symbols, 1, prob = hmm$emissionProbs[states[i], 
                                                          ])
    emission = c(emission, emi)
  }
  return(list(states = states, observation = emission))
}

viterbi <-  function (hmm, observation) 
{
  hmm$transProbs[is.na(hmm$transProbs)] = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations = length(observation)
  nStates = length(hmm$States)
  v = array(NA, c(nStates, nObservations))
  dimnames(v) = list(states = hmm$States, index = 1:nObservations)
  for (state in hmm$States) {
    v[state, 1] = log(hmm$startProbs[state] * hmm$emissionProbs[state, 
                                                                observation[1]])
  }
  for (k in 2:nObservations) {
    for (state in hmm$States) {
      maxi = NULL
      for (previousState in hmm$States) {
        temp = v[previousState, k - 1] + log(hmm$transProbs[previousState, 
                                                            state])
        maxi = max(maxi, temp)
      }
      v[state, k] = log(hmm$emissionProbs[state, observation[k]]) + 
        maxi
    }
  }
  viterbiPath = rep(NA, nObservations)
  for (state in hmm$States) {
    if (max(v[, nObservations]) == v[state, nObservations]) {
      viterbiPath[nObservations] = state
      break
    }
  }
  for (k in (nObservations - 1):1) {
    for (state in hmm$States) {
      if (max(v[, k] + log(hmm$transProbs[, viterbiPath[k + 
                                                        1]])) == v[state, k] + log(hmm$transProbs[state, 
                                                                                                  viterbiPath[k + 1]])) {
        viterbiPath[k] = state
        break
      }
    }
  }
  return(viterbiPath)
}

forward <- function (hmm, observation) 
{
  hmm$transProbs[is.na(hmm$transProbs)] = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations = length(observation)
  nStates = length(hmm$States)
  f = array(NA, c(nStates, nObservations))
  dimnames(f) = list(states = hmm$States, index = 1:nObservations)
  for (state in hmm$States) {
    f[state, 1] = log(hmm$startProbs[state] * hmm$emissionProbs[state, 
                                                                observation[1]])
  }
  for (k in 2:nObservations) {
    for (state in hmm$States) {
      logsum = -Inf
      for (previousState in hmm$States) {
        temp = f[previousState, k - 1] + log(hmm$transProbs[previousState, 
                                                            state])
        if (temp > -Inf) {
          logsum = temp + log(1 + exp(logsum - temp))
        }
      }
      f[state, k] = log(hmm$emissionProbs[state, observation[k]]) + 
        logsum
    }
  }
  return(f)
}

backward <- function (hmm, observation) 
{
  hmm$transProbs[is.na(hmm$transProbs)] = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations = length(observation)
  nStates = length(hmm$States)
  b = array(NA, c(nStates, nObservations))
  dimnames(b) = list(states = hmm$States, index = 1:nObservations)
  for (state in hmm$States) {
    b[state, nObservations] = log(1)
  }
  for (k in (nObservations - 1):1) {
    for (state in hmm$States) {
      logsum = -Inf
      for (nextState in hmm$States) {
        temp = b[nextState, k + 1] + log(hmm$transProbs[state, 
                                                        nextState] * hmm$emissionProbs[nextState, observation[k + 
                                                                                                                1]])
        if (temp > -Inf) {
          logsum = temp + log(1 + exp(logsum - temp))
        }
      }
      b[state, k] = logsum
    }
  }
  return(b)
}

# The occasionally dishonest casino
#-------------
#require(HMM)
nSim = 2000
States = c("Fair", "Unfair")
Symbols = 1:6
transProbs = matrix(c(0.99, 0.01, 0.02, 0.98), c(length(States),
                                                 length(States)), byrow = TRUE)
emissionProbs = matrix(c(rep(1/6, 6), c(rep(0.1, 5), 0.5)),
                       c(length(States), length(Symbols)), byrow = TRUE)
hmm = initHMM(States, Symbols, transProbs = transProbs, emissionProbs = emissionProbs)
sim = simHMM(hmm, nSim)
vit = viterbi(hmm, sim$observation)
f = forward(hmm, sim$observation)
b = backward(hmm, sim$observation)
i <- f[1, nSim]
j <- f[2, nSim]
probObservations = (i + log(1 + exp(j - i)))
posterior = exp((f + b) - probObservations)
x = list(hmm = hmm, sim = sim, vit = vit, posterior = posterior)

# Plotting simulated throws at top
mn = "Fair and unfair die"
xlb = "Throw nr."
ylb = ""

plot(x$sim$observation, ylim = c(-7.5, 6), pch = 3, main = mn,
     xlab = xlb, ylab = ylb, bty = "n", yaxt = "n")
axis(2, at = 1:6)

# Simulated, which die was used (truth)
text(0, -1.2, adj = 0, cex = 0.8, col = "black", "True: green = fair die")
for (i in 1:nSim) {
  if (x$sim$states[i] == "Fair")
    rect(i, -1, i + 1, 0, col = "green", border = NA)
  else rect(i, -1, i + 1, 0, col = "red", border = NA)
}

# Most probable path (viterbi)
text(0, -3.2, adj = 0, cex = 0.8, col = "black", "Most probable path")
for (i in 1:nSim) {
  if (x$vit[i] == "Fair")
    rect(i, -3, i + 1, -2, col = "green", border = NA)
  else rect(i, -3, i + 1, -2, col = "red", border = NA)
}
##################Differences:
text(0, -5.2, adj = 0, cex = 0.8, col = "black", "Difference")
differing = !(x$sim$states == x$vit)
for (i in 1:nSim) {
  if (differing[i])
    rect(i, -5, i + 1, -4, col = rgb(0.3, 0.3, 0.3),
         border = NA)
  else rect(i, -5, i + 1, -4, col = rgb(0.9, 0.9, 0.9),
            border = NA)
}

#Posterior-probability
points(x$posterior[2, ] - 3, type = "l")

# Difference with classification by posterior-probability:
text(0, -7.2, adj = 0, cex = 0.8, col = "black", "Difference by posterior-probability")
differing = !(x$sim$states == x$vit)
for (i in 1:nSim) {
  if (posterior[1, i] > 0.5) {
    if (x$sim$states[i] == "Fair")
      rect(i, -7, i + 1, -6, col = rgb(0.9, 0.9, 0.9),
           border = NA)
    else rect(i, -7, i + 1, -6, col = rgb(0.3, 0.3, 0.3),
              border = NA)
  }
  else {
    if (x$sim$states[i] == "Unfair")
      rect(i, -7, i + 1, -6, col = rgb(0.9, 0.9, 0.9),
           border = NA)
    else rect(i, -7, i + 1, -6, col = rgb(0.3, 0.3, 0.3),
              border = NA)
  }
}