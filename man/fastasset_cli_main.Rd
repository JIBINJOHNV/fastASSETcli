\name{fastasset_cli_main}
\alias{fastasset_cli_main}
\title{Run the fastASSET Command-Line Interface}
\usage{
fastasset_cli_main(arguments = commandArgs(trailingOnly = TRUE))
}
\arguments{
  \item{arguments}{Character vector of command-line arguments.}
}
\value{
Invisibly returns the output-path list, or \code{NULL} after displaying help or
version information.
}
\description{
Programmatic entry point used by the installed \code{fastasset} shell command.
All scientific and computational settings are supplied through CLI options; no
configuration file is read. The command can load an existing GenomicSEM LDSC
object or generate and save one from a validated per-trait summary-statistic
manifest before running FastASSET. The manifest is the only summary-statistic
input and its row order defines the trait order; the wide FastASSET table is
constructed internally. Manifest files may be tabular or GWAS-VCF;
VCF conversion uses bcftools, converts LP to raw P, derives case-control sample
prevalence from NC/NCO when population prevalence is supplied, and uses NEF
directly as total sample size for quantitative traits.
}
