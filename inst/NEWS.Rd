
\name{NEWS}
\title{NEWS file for the \pkg{free1way.docreg} package}

\section{Changes in Version 1.1-0 (2026-10-29)}{
  \subsection{Bugfixes}{
    \itemize{
      \item \code{confint(<free1way>, parm)} returned confidence
            intervals for all parameters when \code{parm} requested a
            subset.
      \item \code{confint(<free1way>, what = "OVL")} now includes
            one for shift effect confidence intervals including zero.
    }
  }
  \subsection{Documentation}{
    \itemize{
      \item Add new \code{Wilcoxon} package vignette.
    }
  }
  \subsection{Features}{
    \itemize{
      \item \code{free1way(, MPL_Jeffreys = TRUE)} implements 
            maximum penalised likelihood estimation with Jeffreys prior
            for all model parameters simultaneously, along with Wald-type
            inference. 
      \item Added experimental \code{scale.test} function for testing
            scale parameters. Parameter estimation, however, is not
            available because the corresponding optimisation problem 
            is non-convex.
    }
  }
  \subsection{Other}{
    \itemize{
      \item Log-likelihood contributions are computed directly on the
            log-scale in a more stable way.
      \item \code{summary(<free1way>)} renamed p-values 
            to \code{Pr()} and thus \code{printCoefmat} recognises them.
      \item \code{rfree1way} also accepts \code{offset} terms equal to the
            number of simulated observations, allowing to sample from
            more complex models.
    }
  }
}

\section{Changes in Version 1.0-3 (2026-09-07)}{
  \subsection{Documentation}{
    \itemize{
      \item Prepare for upcoming changes in penalisation implementation:
            \code{free1way} will perform penalisation only when requested.
    }
  }
}

\section{Changes in Version 1.0-2 (2026-06-30)}{
  \subsection{Bugfixes}{
    \itemize{
      \item Fix missing links in vignette.
    }
  }
}

\section{Changes in Version 1.0-1 (2026-06-16)}{
  \subsection{Bugfixes}{
    \itemize{
      \item \code{strat_ratio != 1} in \code{rfree1way} didn't work as advertised.
      \item Zeros in \code{prob} argument to \code{power.free1way.test} triggered an
            error and are now handled correctly.
    }
  }
}

\section{Changes in Version 1.0-0 (2026-05-18)}{
  \itemize{
    \item Function \code{free1way} implementing semiparametrically efficient inference
          procedures (likelihood- and permutation-based) in distribution-free
          stratified K-sample oneway layouts was added to the \pkg{stats}
          package. Approximate (via \code{power.free1way.test}) and
          simulation-based (utilising \code{rfree1way} for random number
          generation) power-assessment and sample-size planning as well as 
          model diagnostics (available from a \code{plot} method and 
          probability-probability plots in \code{ppplot}) complement the
          functionality. Package \pkg{free1way.docreg} provides additional
          documentation and regression tests.
  }
}
