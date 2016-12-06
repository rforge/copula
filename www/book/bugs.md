---
title: Bug Reporting in R
---

This document explains what to do if you find a bug in the R project, or have a patch that you'd like to submit. It covers:

 1. Making sure your problem is a bug;
 2. Where to submit bug reports (see spam protection **Note**!)
 3. How to write useful and productive bug reports;
 4. How to submit patches;
 5. What to do if something goes wrong.

The intent is to make the most of your time and the R developers' time, by ensuring that bug reports are clear-cut and easy for the developers to respond to.

## Making sure it's a bug

There are a couple of situations where something is definitely a bug; one example of those is when the R session unexpectedly terminates, or *segfaults*. This looks something like:

        *** caught segfault ***
    address (nil), cause 'memory not mapped'

If you're seeing an error like this, unless you have written your own call to compiled code or an internal function (e.g. via `.C` or `.Internal`), it's definitely a bug[^1].

[^1]: If you are making a call to compiled code, check that you have not crashed R by using the wrong argument types (modes).

Other examples of obvious bugs are situations where code does not do what the documentation says it should: either the code is wrong, or the documentation is. One way or another something needs fixing.

Code doing something *unexpected* is not necessarily a bug - make sure to carefully review the documentation for the function you are calling to see if the behaviour it exhibits is what it was designed to do, even if it's not what you want. Similarly, issues with [seemingly-identical numbers not being equal](https://cran.r-project.org/doc/FAQ/R-FAQ.html#Why-doesn_0027t-R-think-these-numbers-are-equal_003f) are known, documented and intractable problems - not bugs.

In all cases, if you think it might be a bug, try launching R from the command line with the `--vanilla` option, to make sure it's a clean session, and see if the bug still appears then. Reduce your code to the minimum required to run the function call where the bug occurred, in particular, only attaching packages required by that call (if any).

If, rather than a bug, you have a suggestion for a new feature in R, you can submit this to the same place you would submit a bug, see the next section.

If your issue doesn't fall into any of these categories - if it's R running slower than expected, say, or something working but not being defined in the way you think would be best, you should consult someone. If you don't know anyone who can look over your code and see if it can be sped up, or if another function would suit your needs better, some useful places to ask for help are:

 1. [The r-help mailing list](https://stat.ethz.ch/mailman/listinfo/r-help); or [other R mailing lists](https://www.r-project.org/mail.html).
 2. [The Bioconductor support site](https://support.bioconductor.org/).
 3. [The R section of StackOverflow](http://stackoverflow.com/questions/tagged/r).
 4. [The R Internet Relay Chat (IRC) channel](http://webchat.freenode.net/?channels=#R).

## Where to submit bug reports and patches

If your problem is definitely a bug - either because it falls into one of the bug categories above, or because you've asked people for help and they've confirmed it's an issue - it's time to submit a report so that it can be fixed.

Depending on the problem, you might need to submit bug reports in different places. The first step is to see which package the function with a bug comes from. The R Core team only maintains the core language and the R packages
labelled with `Maintainer: R Core Team <R-core@r-project.org>`. You can see
this label by running e.g. `maintainer("graphics")` in R.

If your bug is in `somePkg` and that is not maintained by the R Core team, you should submit your report to the package maintainer. Running `bug.report(package = "somePkg")`
directs you to the right place,
either opening the relevant bug tracking web page or helping you to compose an email to the package maintainer.

The `bug.report` function is disabled in some IDEs (e.g. RStudio) to avoid misuse; to identify the right place to submit a report yourself, start by looking at the output of `packageDescription("somePkg")`,
the R help for the package, or the webpage for the package from the relevant repository, e.g. [CRAN](https://cran.r-project.org/web/packages/available_packages_by_name.html) or [Bioconductor](https://www.bioconductor.org/packages/release/BiocViews.html#___Software). Some packages have a bug submission page, such as an issue tracker on GitHub, listed under the `BugReports` field in the package description. If you follow this link you may find your bug has already been reported, otherwise you can submit your report there, following the guidelines on bug reporting discussed below. If there is no bug submission page, you should email your bug report to the package `Maintainer` via the address in the package description.

If your bug is in the language, though, or the Core-supported packages, you should submit your report to R's [Bugzilla](https://bugs.r-project.org/bugzilla3/).  
**NOTE: due to abuse by spammers, since 2016-07-09 only users who have previously submitted bugs can submit new ones on R's Bugzilla.  We're working on a better system...  In the mean time, post (e-mail) to `R-devel` or ask an R Core member to add you manually to R's Bugzilla members.**  
It is important to try to make sure that the report isn't extraneous. The easiest way to do this is to first look at the [upcoming changes in R](https://svn.r-project.org/R/trunk/doc/NEWS.Rd), to see if the bug has already been patched (just not released yet), and to [browse the latest bug reports](https://bugs.r-project.org/bugzilla/buglist.cgi?bug_file_loc_type=allwordssubstr&bug_status=NEW&bug_status=ASSIGNED&bug_status=REOPENED&bug_status=UNCONFIRMED&bugidtype=include&chfieldto=Now&cmdtype=doit&emailassigned_to1=1&emailassigned_to2=1&emailcc2=1&emailreporter2=1&emailtype1=substring&emailtype2=substring&field0-0-0=noop&long_desc_type=substring&order=bugs.delta_ts%20desc&query_format=advanced&short_desc_type=allwordssubstr&type0-0-0=noop) or [search for the bug](https://bugs.r-project.org/bugzilla/query.cgi) in Bugzilla to see if (even if it hasn't been patched yet) it has been reported. If your bug has not yet been reported or fixed, you can report the bug following the guidelines in the section [Writing a good bug report](#writing-a-good-bug-report). If you have a patch accompanying your bug, see the section [How to submit patches](#how-to-submit-patches).


If you wish to submit a feature request, rather than a bug report, your best bet is to ask about it first on the [r-devel](https://stat.ethz.ch/mailman/listinfo/r-devel) mailing list. If the feedback is positive, you can submit your suggestion using the bug reporting form on Bugzilla, where you should select `Wishlist` in the `Component` field and start your summary with `Wishlist:`.

Issues related to message translations should be sent to the last translator or to the relevant [translation team](https://developer.r-project.org/TranslationTeams.html). To find the last translator, you will need to look at the comments at the top of the relevant `.po` file in the R source code, for example German translations of messages in the base package are in `src/library/base/po/R-de.po`. You can download the R source code from CRAN, or otherwise [browse the R-devel sources](https://svn.r-project.org/R/trunk/) or [their mirror on GitHub](https://github.com/wch/r-source).

## Writing a good bug report

Bug reports should include a way of [reproducing](https://en.wikipedia.org/wiki/Reproducibility) the bug. This should be as simple as possible. If the person trying to fix the bug can't work out how to make it appear, or has to jump through a lot of unnecessary hoops to make it appear, you're going to waste a lot of their time.

Bugzilla is maintained by a small number of people, so it's best to make sure your bug report is clear and well-written. If it's not, it will suck in more energy from the maintainers and take longer for the bug to get fixed - or it may end up not getting handled at all. In particular, you should:

 1. Write a clear and unique summary for the bug. “Stopping a run of lm() causes a crash” is good; “software crashes” is not.
 2. Include, in the description, the steps to reproduce the bug mentioned above.
 3. Identify both what happened (“the software crashed”) and what you expected to happen (“lm() should stop running”).
 4. Identify the platform, architecture, and version of R where you found this bug. You can retrieve that from within R by typing `R.version`.
 5. Just focus on the facts of what happened, rather than on your theories of what the bug is and where it comes from.

At that point, you've written a good bug report! Sit back and wait for a
developer to respond to it. If you run into any problems with that
response, see the section
[What to do if there's an issue](#what-to-do-if-there-is-an-issue).

## How to submit patches

Sometimes you'll find a bug and also see, from looking at the code, how to fix it. When this happens you've got the opportunity to submit a patch, which can reduce the workload on the R developers: they get to test, tweak and include the code instead of having to write it all from scratch.

To prepare a patch, you're going to need the latest developer version of R. This is maintained in a [Subversion](http://subversion.apache.org/) (SVN) repository. Once you've got SVN installed on your system, open the command line and type:

    svn checkout https://svn.r-project.org/R/trunk/ R-devel

This should create a directory, `R-devel`, in your current working directory. This contains the source code for the newest version of R.

Go through and make the changes you need to make in order to patch the bug - try to keep to whatever coding style and conventions the functions you're changing use, just to make things easier. Once you're done, go to the `R-devel` directory in your terminal and type:

    svn update
    svn diff > patch.diff

This updates the code then creates a new file, `patch.diff`, that contains the changes between the latest version of R, and your alterations. And that's a patch! Just attach that to the bug report you're writing, note in the report that there's an associated patch, and you're done.

## What to do if there is an issue

In an ideal world you write an informative bug report (and maybe submit a patch), someone comes along promptly and fixes it, and everyone is happy. In the world we've got, the people maintaining R have a lot of responsibilities, and all of them are doing this work as volunteers. This means that, practically speaking, bugs may take a very long time to get fixed, accidentally get missed, or result in an unexpected or unpleasant outcome - not out of any maliciousness but simply because the people responsible for the software can get pretty stressed.

If you experience technical issues with R's Bugzilla that do not resolve themselves after a period of time, you should contact the current maintainer:  [simon.urbanek@R-project.org](mailto:simon.urbanek@R-project.org). If you feel like your bug has been missed (e.g. because a new release of R has come out, and it was not fixed), you can bring attention to it by simply adding a comment like "This is still present in the x.y.z release" on Bugzilla.
Even better would be to install a pre-release alpha or beta version to confirm it is still present, and report that.

If you feel it has been assessed wrongly, you can leave a comment to that effect on Bugzilla.  If you are personally acquainted with a member of R Core you could contact them directly.  In either case, present your case clearly, and respect the fact that the R Core members may judge the importance of the issue (or even whether it is a bug or not) differently than you do.
