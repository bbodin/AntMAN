## Test environments

* local R installation (f33), R 4.0.5
* win-builder (devel)

## R CMD check results

There are no ERRORs or WARNINGs.

## Resubmission

### Previous Feedback - 29th June 2021

No need to comment, you can \donttest{} them.

We also see:
   The Date field is over a month old.

#### Update - 1st July 2021

 * Update the date in DESCRIPTION
 * Used \donttest{} and \dontrun{} where necessary


### Last feedback - 1st July 2021
Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags in up to 30 .Rd files, e.g.:
      AM_mcmc_configuration.Rd: \value
      AM_mcmc_output.Rd: \value
      ...

You have examples for unexported functions.
    IAM_compute_stirling_ricor_abs() in:
       IAM_compute_stirling_ricor_abs.Rd
    IAM_compute_stirling_ricor_log() in:
       IAM_compute_stirling_ricor_log.Rd
       IAM_VnkDelta.Rd
       IAM_VnkNegBin.Rd
       IAM_VnkPoisson.Rd
    IAM_VnkDelta() in:
       IAM_VnkDelta.Rd
    IAM_VnkNegBin() in:
       IAM_VnkNegBin.Rd
    IAM_VnkPoisson() in:
       IAM_VnkPoisson.Rd
  Please either omit these examples or export the functions.

\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user.
Does not seem necessary.

Please unwrap the examples if they are executable in < 5 sec, or replace
\dontrun{} with \donttest{}.

Please make sure that you do not change the user's options, par or
working directory. If you really have to do so within functions, please
ensure with an *immediate* call of on.exit() that the settings are reset
when the function is exited. e.g.:
...
oldpar <- par(no.readonly = TRUE)    # code line i
on.exit(par(oldpar))            # code line i + 1
...
par(mfrow=c(2,2))            # somewhere after
...
e.g.: AM_plot.R, AM_prior.R
If you're not familiar with the function, please check ?on.exit. This
function makes it possible to restore options before exiting a function
even if the function breaks. Therefore it needs to be called immediately
after the option change within a function.


Please always make sure to reset to user's options(), working directory
or par() after you changed it in examples and vignettes and demos.
e.g.: tests
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)

Please do not set a seed to a specific number within a function. e.g.:
AM_demo.R


Additionally:
Have the issues why your package was archived been fixed?
Please explain this in the submission comments.

#### Update - 20th July 2021

 - We added return values everywhere it was needed.
 - We removed examples from unexported functions.
 - We removed dontrun{} or replaced it by donttest{} when necessary.
 - We made sure that any used of par() in the package functions is protected by on.exit()
 - We removed seed() from the package functions
 - To answer the question regarding the package being archived, we remove the include of openmp header and dont specify openmp options anymore. 