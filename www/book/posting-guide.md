---
title: "Posting Guide: How to ask good questions that prompt useful answers"
---

This guide is intended to help you get the most out of the R mailing lists, and to avoid embarrassment. Like many responses posted on the list, it is written in a concise manner. This is not intended to be unfriendly - it is more a consequence of allocating the limited available time and space to technical issues rather than to social niceties.

**The list:** Remember that R is free software, constructed and maintained by volunteers. They have various reasons for contributing software and participating on the mailing lists, but often have limited time.

**Good manners:** Remember that customs differ. Some people are very direct. Others surround everything they say with hedges and apologies. Be tolerant. Rudeness is never warranted, but sometimes \`read the manual' *is* the appropriate response. Don't waste time discussing such matters on the list. Ad hominem comments are absolutely out of place.

**Questions about statistics:** The R mailing lists are primarily intended for questions and discussion about the R software. However, questions about statistical methodology are sometimes posted. If the question is well-asked and of interest to someone on the list, it *may* elicit an informative up-to-date answer. See also the Usenet groups sci.stat.consult (applied statistics and consulting) and sci.stat.math (mathematical stat and probability).

**Basic statistics and classroom homework:** R-help is not intended for these.

**Which list: R-help, R-package-devel, R-devel, or Bioconductor?**
There have been two, now three widely used mailing lists for questions and discussion about R and a list dedicated more specifically to issues in the use of bioconductor packages and bioinformatics.
**R-help** is intended to be comprehensible to people who want to use R to solve problems but who are not necessarily interested in or knowledgeable about programming.
**R-package-devel** has been created in May 2015 specifically to help R package authors in writing and checking their R packages, notably those meant for CRAN.
**R-devel** is intended for questions and discussion about R development and programming.
Questions likely to prompt discussion unintelligible to non-programmers should rather
go to R-devel than R-help.
Questions about _package development_, however, e.g., passing `R CMD check` should go
to `R-package-devel` rather than R-devel.
For example, questions involving C, C++, etc. code should go to R-devel. More general
questions involving pure R code and questions likely to be of interest to
the large and diverse set of subscribers to R-help should go to R-help.\
**Bioconductor** is for announcements about the development of [Bioconductor](http://www.bioconductor.org/),
availability of new code, questions and answers about problems and solutions using
Bioconductor, etc. See [Bioconductor mailing lists](http://bioconductor.org/help/mailing-list/)
for details. See [below](#contrib_pkgs) for questions on *contributed packages*.

**Platform-specific questions**: There are lists **R-sig-Mac**, **R-sig-Debian** and **R-sig-Fedora** for R on Mac OS X, Debian/Ubuntu and Fedora/Redhat respectively. Questions specific to those platforms (especially *re* installation and the **R.app** GUI on Mac OS X) are more likely to get informed responses on the appropriate list, and that is certainly the place to discuss possible bugs.

**Do your homework before posting:** If it is clear that you have done basic background research, you are far more likely to get an informative response. See also [Further Resources](#further) further down this page.

-   Do `help.search("keyword")` and `apropos("keyword")` with different keywords (type this at the R prompt).
-   Do `RSiteSearch("keyword")` with different keywords (at the R prompt) to search R functions, contributed packages and R-Help postings. See `?RSiteSearch` for further options and to restrict searches.
-   Read the online help for relevant functions (type `?functionname`, e.g., `?prod`, at the R prompt)
-   If something seems to have changed in R, look in the latest [NEWS](http://cran.r-project.org/src/base/NEWS) file on CRAN for information about it.
-   Search the R-faq and the R-windows-faq if it might be relevant (<http://cran.r-project.org/faqs.html>)
-   Read at least the relevant section in [An Introduction to R](http://cran.r-project.org/doc/manuals/R-intro.pdf)
-   If the function is from a package accompanying a book, e.g., the MASS package, consult the book before posting
-   The R Wiki has a [section on finding functions and documentation](http://rwiki.sciviews.org/doku.php?id=getting-started:reference-cards:getting-help)

**Technical details of posting**: See [General Instructions](http://www.r-project.org/mail.html#instructions) for more details of the following:

-   No HTML posting (harder to detect spam) (note that this is the default in some mail clients - you may have to turn it off). Note that chances have become relatively high for 'HTMLified' e-mails to be completely intercepted (without notice to the sender).
-   No binary attachments except for PS, PDF, and some image and archive formats (others are automatically stripped off because they can contain malicious software). Files in other formats and larger ones should rather be put on the web and have only their URLs posted. This way a reader has the option to download them or not.
-   Use an informative subject line (not something like \`question')
-   For new subjects, compose a new message and include the 'r-help@R-project.org' (or 'r-devel@R-project.org') address specifically. (Replying to an existing post and then changing the subject messes up the threading in the archives and in many people's mail readers.)
-   If you can't send from an email address that simply accepts replies, then say so in your posting so that people are not inconvenienced when they try to respond to your message
-   Some consider it good manners to include a concise signature specifying affiliation

**Surprising behavior and bugs:** What you think is a bug may be many other things, such as a default behavior that you do not like, a feature, an undocumented feature, or a bug in the documentation. You do not need to commit yourself to one of these in order to ask a question. If it is a real bug, someone will notice it. Before you post a real bug report, make sure you read [R Bugs](http://CRAN.R-project.org/doc/FAQ/R-FAQ.html#R-Bugs) in the R-faq. If you're not completely and utterly sure something is a bug, post a question to r-help, not a bug report to r-bugs - every bug report requires manual action by one of the R-core members. Also, see Simon Tatham's essay on [How to Report Bugs Effectively](http://www.chiark.greenend.org.uk/~sgtatham/bugs.html).

For questions about unexpected behavior or a possible bug, you should, at a minimum, copy and paste the output from `sessionInfo()` into your message. When mentioning version numbers, always use the full version number, e.g., \`2.6.1', not just \`2.6', and also mention the platform (Windows, Linux, MacOS X, with their versions). Other potentially relevant details include the locale (type `Sys.getlocale()` at the R prompt), and whether you installed a pre-compiled binary version of R or compiled it yourself. If the function is in a package other than \`base', include the header output from `library(help=thatPackage)`. If you are using an old version of R and think it does not work properly, upgrade to the latest version and try that, before posting. If possible, try the current R-patched or R-devel version of R (see the FAQ for details), to see if the problem has already been addressed.

For questions about functions in standard packages distributed with R (see the FAQ [Add-on packages in R](http://cran.r-project.org/doc/FAQ/R-FAQ.html#Add-on-packages-in-R)), ask questions on R-help.\
 If the question relates to a *contributed package* , e.g., one downloaded from CRAN, try contacting the package maintainer first. You can also use `find("functionname")` and `packageDescription("packagename")` to find this information. **Only** send such questions to R-help or R-devel if you get no reply or need further assistance. This applies to both requests for help and to bug reports.

Don't say *\`R crashed'*, which you may take to mean that R gave an error and terminated your piece of code, but most people will take to mean abnormal termination of the R program. Say exactly what happened, including any error messages you received.

**Examples:** Sometimes it helps to provide a *small* example that someone can actually run. For example:

      If I have a matrix x as follows:
      > x <- matrix(1:8, nrow=4, ncol=2,
                    dimnames=list(c("A","B","C","D"), c("x","y"))
      > x
        x y
      A 1 5
      B 2 6
      C 3 7
      D 4 8
      >

      how can I turn it into a dataframe with 8 rows, and three
      columns named `row', `col', and `value', which have the
      dimension names as the values of `row' and `col', like this:
      > x.df
         row col value
      1    A   x      1
       ...
      (To which the answer might be:
      > x.df <- reshape(data.frame(row=rownames(x), x), direction="long",
                        varying=list(colnames(x)), times=colnames(x),
                        v.names="value", timevar="col", idvar="row")
      )

When providing examples, it is best to give an R command that constructs the data, as in the `matrix()` expression above. For more complicated data structures, `dump("x", file=stdout())` will print an expression that will recreate the object `x`.

**Further resources:** not always consulting these before posting is OK, except perhaps if you have a very general question and the topic of the document indicates immediate relevancy. A response might just point you towards one of these.

-   [Writing R Extensions](http://cran.r-project.org/doc/manuals/R-exts.pdf) covers packages, to and from C & C++, R help files.
-   [R Data Import/Export](http://cran.r-project.org/doc/manuals/R-data.pdf) covers reading from and exporting to other file formats
-   [R Language Definition](http://cran.r-project.org/doc/manuals/R-lang.pdf) (aka \`R Language Manual') describes the R language, data objects, etc. This document specifies many details of R that are not covered in the help pages or in "An Introduction to R".
-   The article \`R Help Desk' by Uwe Ligges in the [R Newsletter Vol 3 No 1](http://cran.r-project.org/doc/Rnews/Rnews_2003-1.pdf) has useful tips about where to find information about R.
-   Two books are often cited in answers to questions: [Venables & Ripley, 2002, Modern Applied Statistics with S (4th ed)](http://www.stats.ox.ac.uk/pub/MASS4) ("MASS") and
 [Pinheiro & Bates, 2000](http://link.springer.com/book/10.1007%2Fb98882)
<!-- this is DOI:10.1007/b98882 -->
 (about lme in particular). If you can get these books, look at them. If you cannot, you might indicate this in your question. Further books relevant to R are listed at [Publications](http://www.r-project.org/doc/bib/R-publications.html).
-   The [Contributed Documentation](http://cran.r-project.org/other-docs.html) page lists a number of useful introductory and other documents available online.
-   The R-help archives (see under \`Archives and Search Facilities' at <http://www.r-project.org/mail.html>)
-   A Google search of the web in general or the r-project site specifically (put \`site:r-project.org' at the beginning of the search terms)

**Responding to other posts:**

-   In general, respond both to the list and to the poster. This will result in your response being archived and the poster will see your answer.
-   Rudeness and ad hominem comments are not acceptable. Brevity is OK.
-   Take care when you quote other people's comments to respect their rights, e.g., as summarized [here](http://www.jiscmail.ac.uk/policyandsecurity/copyrightissues.html). In particular:
    -   Private messages should never be quoted without permission.
    -   The original authorship and meaning should always be clear.
-   It is common practice to quote prior messages in entirety.
-   If you edit the quoted text from prior messages, include sufficient context so that the content is clear even if your archived message is referenced at a much later date and in isolation from other messages.
-   When responding to a very simple question, use the following algorithm:
    1.  compose your response
    2.  type `4*runif(1)` at the R prompt, and wait this many hours
    3.  check for new posts to R-help; if no similar suggestion, post your response

    (This is partly in jest, but if you know immediately why it is suggested, you probably should use it! Also, it's a nice idea to replace 4 by the number of years you have been using R or S-plus.)

**Common posting mistakes:**

-   Missing or non-specific message subject
-   Not doing basic homework before posting a question
-   Not being more specific than \`function xxx doesn't work'
-   Being overly specific and not stating your real goal.
-   Not including version and OS information in question regarding unexpected behavior.
-   Finding a bug in an old version of R that has been fixed in the most recent version.
-   Claiming that something is a bug when in fact the software is working as intended and documented, just not in the way you first expected.
-   Claiming that some commonly used function is not behaving in a sensible manner (if you find the behavior odd, it's more productive and polite to ask why it behaves the way it does - after reading the relevant documentation.)
-   Threatening not to use the software if you cannot get your question answered. Even when intended as a statement of fact, this tends to create negative attitudes.

**Final words:**

It is a skill to ask good questions. If at first you don't get the answers that are useful to you, don't get discouraged. A response that is concise and technically accurate may be just that, and not an intended putdown. If you feel insulted by some response to a post of yours, don't make any hasty response in return - you're more likely than not to regret it. Read Eric Raymond's essay [How To Ask Questions The Smart Way](http://www.catb.org/~esr/faqs/smart-questions.html) for more suggestions, and for insight into people's behavior on technical mailing lists (but don't try asking people at catb.org questions about R).

Posters should be aware that the R lists are *public* discussion lists and anything you post will be **archived and accessible** via several websites for many years.

*Compiled by Tony Plate (tplate at acm dot org), most recently updated July, 2013*

