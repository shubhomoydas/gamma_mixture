# gamma_mixture

Mixture of Gamma Distributions
-------------
Let us assume we have anomaly scores from a homogeneous ensemble of anomaly detectors (such as an ensemble of Gaussian Mixture Models). Higher scores mean more anomalous. It is reasonable to expect that one data instance would have similar scores assigned to it across all ensemble members, but noisy; and true anomalies would likely be assigned higher scores by most ensemble members.

In this situation, we might model the distribution of scores from each anomaly detector as a mixture of two Gamma distributions. We then jointly infer the probability of an instance being anomalous along with the parameters of the Gamma distributions. We hope that the joint inference would be able to gather more robust evidence of an instance being anomalous (or not) based on stronger statistical support.

The derivations for E-M inference are at: https://github.com/shubhomoydas/rank_aggregation_mixtures/blob/master/documentation/ModelBasedRankAggregation-gamma.pdf

Mixture of Gaussian Distributions
-------------
The same idea can also be applied with the simpler assumption that the scores have mixtures of Gaussian distributions. The derivations are at https://github.com/shubhomoydas/rank_aggregation_mixtures/blob/master/documentation/ModelBasedRankAggregation-normal-diag-covs.pdf, if it helps.

Discussion
-------------
In practice this method does not work as well. It is more likely that the anomaly scores would fit Generalized Extreme Value distributions and not Gamma.
