\name{SEARS}
\alias{SEARS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A randomized distributed phase I-II seamless dose escalation/expansion schema for dose optimization and selection in early oncology clinical development
}
\description{
Implements seamless randomized phase I-II SEARS design for finding the optimal design. Practitioners can use a rich set of parameters to explore various real scenarios for their studies. The function can generate operating characteristics via simulation for practitioners to examine the design’s properties
}
\usage{
SEARS(p.p, p.d, p.tox, k1, k2, pi_t, pi_e, pT, eff_a = 0.5, eff_b = 0.5,
      plac_a = 0.5, plac_b = 0.5, tox_a = 1, tox_b = 1, csize, csize2,
      p.star = 0.2, q.star = 0.6, f.star = 0.06, p.star2 = 0.2,
      q.star2 = 0.98, d.cs, p.cs, phase1_size, n_earlystop,
      extrasafe_BOIN = FALSE, offset_BOIN = 0.05, Nsim, n_catchup,
      control_arm = "", power_c = 0.5, lower_bound = 0.05, weight1,
      weight2, seed = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.p}{
the true placebo response rate
  }
  \item{p.d}{
the true dose response rate vector
  }
  \item{p.tox}{
the true dose toxicity rate vector
  }
  \item{k1}{
the safety rule cutoff value in phase I
  }
  \item{k2}{
the safety rule cutoff value in phase II
  }
  \item{pi_t}{
the physician-specified upper toxicity rate threshold
  }
  \item{pi_e}{
the physician-specified lower response rate threshold
  }
  \item{pT}{
the target toxicity rate
  }
  \item{eff_a}{
the hyperparameter "a" for priors of the response rate for experimental doses. The default value is \code{eff_a = 0.5}
  }
  \item{eff_b}{
the hyperparameter "b" for priors of the response rate for experimental doses. The default value is \code{eff_b = 0.5}
  }
  \item{plac_a}{
the hyperparameter "a" for prior of the response rate for the control arm. The default value is \code{plac_a = 0.5}
  }
  \item{plac_b}{
the hyperparameter "b" for prior of the response rate for the control arm. The default value is \code{plac_b = 0.5}
  }
  \item{tox_a}{
the hyperparameter "a" for priors of the toxicity rates. The default value is \code{tox_a = 1}
  }
  \item{tox_b}{
the hyperparameter "b" for priors of the toxicity rates. The default value is \code{tox_b = 1}
  }
  \item{csize}{
the cohort size in phase I
  }
  \item{csize2}{
the cohort size in phase II
  }
  \item{p.star}{
the fixed cutoff probability for toxicity in phase I. The default value is \code{p.star = 0.2}
  }
  \item{q.star}{
the fixed cutoff probability for efficacy in phase I. The default value is \code{q.star = 0.6}
  }
  \item{f.star}{
the small probability cutoff. It will be used for futility dose exclusion. The default value is \code{f.star = 0.06}
  }
  \item{p.star2}{
the fixed cutoff probability for toxicity in phase II. The default value is \code{p.star2 = 0.2}
  }
  \item{q.star2}{
the fixed cutoff probability for efficacy in phase II. The default value is \code{q.star2 = 0.98}
  }
  \item{d.cs}{
the prespecified maximum allowable number of patients enrolled for each dose
  }
  \item{p.cs}{
the prespecified maximum allowable number of patients enrolled for the placebo
  }
  \item{phase1_size}{
the prespecified maximum sample size of phase I trial
  }
  \item{n_earlystop}{
the cutoff number in phase I. When the number of patients enrolled at a certain dose reaches this value in phase I, this dose will be graduated to phase II
  }
  \item{extrasafe_BOIN}{
the logical value which indicates whether a more stringent stopping rule will be applied to
phase I BOIN design. The default value is \code{extrasafe_BOIN = FALSE}
  }
  \item{offset_BOIN}{
the small positive number (between 0 and 0.5) to control how strict the stopping rule is when \code{extrasafe_BOIN = TRUE}. The default value is \code{offset_BOIN = 0.05}
  }
  \item{Nsim}{
the number of simulated trials
  }
  \item{n_catchup}{
the catch-up cutoff when employing the adaptive randomization in phase II
  }
  \item{control_arm}{
the argument is for phase II design only. If this argument is "fixed", then allocation probability of control arm (the first component of the allocation probability vector)
will be fixed to \eqn{\frac{1}{K}}. The default of this argument will return unfixed results; K indicates total number of arms (including control arm)
  }
  \item{power_c}{
the power correction parameter of the allocation probability. The default value is \code{power_c = 0.5}
  }
  \item{lower_bound}{
the lower bound of the allocation probability in phase II design. It must be a value between 0 and \eqn{\frac{1}{K}}. The default value is \code{lower_bound = 0.05}; K indicates total number of arms (including the control arm)
  }
  \item{weight1}{
the penalized weight in the utility function for the toxicity
  }
  \item{weight2}{
the additional penalized weight in utility function for dose(s) which has (have) the toxicity probability greater than the pre-specified DLT rate
  }
  \item{seed}{
the seed. The default value is \code{seed = 100}
  }
}
\value{
\code{SEARS()} returns a list with following elements (1) type I error (2) average sample size for the trial (3) average sample size for each dose (4) average sample size for placebo (5) selection percentage for each dose (6) average toxicity events for each dose
}
\author{
Chia-Wei Hsu, Haitao Pan
}
\examples{
SEARS(p.p = 0.2, p.d = c(0.2, 0.2, 0.2, 0.2, 0.2), p.tox = c(0.03, 0.06, 0.17, 0.3, 0.5),
      k1 = 0.95, k2 = 0.8, pi_t = 0.17, pi_e = 0.2, pT = 0.17, eff_a = 0.5, eff_b = 0.5,
      plac_a = 0.5, plac_b = 0.5, tox_a = 1, tox_b = 1, csize = 3, csize2 = 3,
      p.star = 0.2, q.star = 0.6, f.star = 0.06, p.star2 = 0.2, q.star2 = 0.98,
      d.cs = 36, p.cs = 36, phase1_size = 30, n_earlystop = 100, extrasafe_BOIN = FALSE,
      offset_BOIN = 0.05, Nsim = 10, n_catchup = 3, control_arm = "fixed", power_c = 0.5,
      lower_bound = 0.05, weight1 = 0.5, weight2 = 0.5, seed = 100)
}
