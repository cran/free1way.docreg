
\name{NEWS}
\title{NEWS file for the \pkg{free1way.docreg} package}

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
