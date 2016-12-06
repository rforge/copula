---
title: Getting Help with R
---

## Helping Yourself

Before asking others for help, it's generally a good idea for you to try to
help yourself. R includes extensive facilities for accessing documentation
and searching for help. There are also specialized search engines for
accessing information about R on the internet, and general internet search
engines can also prove useful ([see below](#r-help-on-the-internet)).

### R Help: `help()` and `?`

The `help()` function and `?` help operator in R provide access to the
documentation pages for R functions, data sets, and other objects,
both for packages in the standard R distribution and for contributed
packages. To access documentation for the standard `lm` (linear model)
function, for example, enter the command `help(lm)` or `help("lm")`,
or `?lm` or `?"lm"` (i.e., the quotes are optional).

To access help for a function in a package that's *not* currently loaded,
specify in addition the name of the package: For example, to obtain
documentation for the `rlm()` (robust linear model) function in the
**MASS** package, `help(rlm, package="MASS")`.

Standard names in R consist of upper- and lower-case letters, numerals
(`0-9`), underscores (`_`), and periods (`.`), and must begin with a letter
or a period. To obtain help for an object with a *non-standard* name (such
as the help operator `?`), the name must be quoted: for example,
`help('?')` or `?"?"`.

You may also use the `help()` function to access information about a
package in your library --- for example, `help(package="MASS")` --- which
displays an index of available help pages for the package along with some
other information.

Help pages for functions usually include a section with executable examples
illustrating how the functions work. You can execute these examples in the
current R session via the `example()` command: e.g., `example(lm)`.

### Vignettes and Code Demonstrations: `browseVignettes()`, `vignette()` and `demo()`

Many packages include *vignettes*, which are discursive documents
meant to illustrate and explain facilities in the package. You can
discover vignettes by accessing the help page for a package, or via
the `browseVignettes()` function: the command `browseVignettes()`
opens a list of vignettes from *all* of your installed packages in
your browser, while `browseVignettes(package=package-name)` (e.g.,
`browseVignettes(package="survival")`) shows the vignettes, if any,
for a particular package. `vignette()` is employed similarly, but
displays a list of vignettes in text form.

You can also use the `vignette("vignette-name")` command to view a vignette
(possibly specifying the name of the package in which the vignette resides,
if the vignette name is not unique): for example, `vignette("timedep")` or
`vignette("timedep", package="survival")` (which are, in this case,
equivalent).

Vignettes may also be accessed from the CRAN page for the package
(e.g. [survival](https://cran.r-project.org/web/packages/survival/index.html)),
if you wish to review the vignette for a package prior to installing and/or
using it.

Packages may also include extended code demonstrations ("demos"). The
command `demo()` lists all demos for all packages in your library, while
`demo(package="package-name")` (e.g., `demo(package="stats")`) lists demos
in a particular package. To run a demo, call the `demo()` function with the
quoted name of the demo (e.g., `demo("nlm")`), specifying the name of the
package if the name of the demo isn't unique (e.g., `demo("nlm",
package="stats")`, where, in this case, the package name need not be given
explicitly).

### Searching for Help Within R

The `help()` function and `?` operator are useful only if you already know
the name of the function that you wish to use. There are also facilities in
the standard R distribution for discovering functions and other
objects. The following functions cast a progressively wider net. Use the
help system to obtain complete documentation for these functions: for
example, `?apropos`.

#### `apropos()`

The `apropos()` function searches for objects, including functions,
directly accessible in the current R session that have names that include a
specified character string. This may be a literal string or a *regular
expression* to be used for pattern-matching (see `?"regular
expression"`). By default, string matching by `apropos()` is
case-insensitive. For example, `apropos("^glm")` returns the names of all
accessible objects that start with the (case-insensitive) characters
`"glm"`.

#### `help.search()` and `??`

The `help.search()` function scans the documentation for packages installed
in your library. The (first) argument to `help.search()` is a character
string or regular expression. For example, `help.search("^glm")` searches
for help pages, vignettes, and code demos that have help "aliases,"
"concepts," or titles that begin (case-insensitively) with the characters
`"glm"`. The `??` operator is a synonym for `help.search()`: for example,
`??"^glm"`.

#### `RSiteSearch()`

`RSiteSearch()` uses an
[internet search engine](http://search.r-project.org) (also see
[below](#r-help-on-the-internet)) to search for information in function
help pages and vignettes for all CRAN packages, and in CRAN task views
(described [below](#cran-task-views)). Unlike the `apropos()` and
`help.search()` functions, `RSiteSearch()` requires an active internet
connection and doesn't employ regular expressions. Braces may be used to
specify multi-word terms; otherwise matches for individual words are
included. For example, `RSiteSearch("{generalized linear model}")` returns
information about R functions, vignettes, and CRAN task views related to
the term `"generalized linear model"` without matching the individual words
`"generalized"`, `"linear"`, or `"model"`.

`findfn()` and `???` in the **sos** package, which is *not* part of the
standard R distribution but is available on CRAN, provide an alternative
interface to `RSiteSearch()`.

#### `help.start()`

`help.start()` starts and displays a hypertext based version of R's online
documentation in your default browser that provides links to locally
installed versions of the R manuals, a listing of your currently installed
packages and other documentation resources.

## R Help on the Internet

There are internet search sites that are specialized for R searches,
including [search.r-project.org](http://search.r-project.org/) (which is
the site used by `RSiteSearch`) and [Rseek.org](http://www.rseek.org/).

It is also possible to use a general search site like
[Google](http://google.com/), by qualifying the search with "R" or the name
of an R package (or both). It can be particularly helpful to paste an error
message into a search engine to find out whether others have solved a
problem that you encountered.

### CRAN Task Views

CRAN Task Views are documents that summarize R resources on CRAN in
particular areas of application, helping your to navigate the maze of
thousands of CRAN packages. A
[list of available Task Views](https://cran.r-project.org/web/views/) may
be found on CRAN.

### R FAQs (Frequently Asked Questions)

There are three primary FAQ listings which are periodically updated to
reflect very commonly asked questions by R users. There is a
[Main R FAQ](https://cran.r-project.org/doc/FAQ/R-FAQ.html), a
[Windows specific R FAQ](https://cran.r-project.org/bin/windows/base/rw-FAQ.html)
and a
[Mac OS (OS X) specific R FAQ](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html).


## Asking for Help

If you find that you can't answer a question or solve a problem yourself,
you can ask others for help, either locally (if you know someone who is
knowledgeable about R) or on the internet. In order to ask a question
effectively, it helps to phrase the question clearly, and, if you're trying
to solve a problem, to include a small, self-contained, reproducible
example of the problem that others can execute. For information on how to
ask questions, see, e.g., the R mailing list
[posting guide](https://www.r-project.org/posting-guide.html), and the
document about
[how to create reproducible examples for R](http://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example)
on Stack Overflow.

### Stack Overflow

[Stack Overflow](http://stackoverflow.com) is a well organized and
formatted site for help and discussions about programming. It has
excellent searchability. Topics are tagged, and ["r" is a very popular
tag on the site](http://stackoverflow.com/tags/r/info) with almost
150,000 questions (as of summer 2016).  To go directly to R-related
topics, visit 
[http://stackoverflow.com/questions/tagged/r](http://stackoverflow.com/questions/tagged/r).
For an example both of the value of the site's organization and
information that is very useful to R users, see 
["How to make a great R reproducible example?"](http://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example), 
which is also mentioned above.



### R Email Lists

The R Project maintains a number of subscription-based
[email lists](https://www.r-project.org/mail.html) for posing and answering
questions about R, including the general
[R-help](https://stat.ethz.ch/mailman/listinfo/r-help) email list, the
[R-devel](https://stat.ethz.ch/mailman/listinfo/r-devel) list for R code
development, and
[R-package-devel](https://stat.ethz.ch/mailman/listinfo/r-package-devel)
list for developers of CRAN packages; lists for announcements about
[R](https://stat.ethz.ch/mailman/listinfo/r-announce) and
[R packages](https://stat.ethz.ch/mailman/listinfo/r-packages); and a
variety of more specialized lists. Before posing a question on one of these
lists, please read the
[R mailing list instructions](https://www.r-project.org/mail.html) and the
[posting guide](https://www.r-project.org/posting-guide.html).

