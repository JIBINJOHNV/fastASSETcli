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
configuration file is read.
}
