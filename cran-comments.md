## Test environments

* macOS 11.1 Big Sur, R 4.0.3 (2020-10-10) (local)
* Ubuntu Linux 16.04 LTS, R-release, GCC (r-hub)
* Fedora Linux, R-devel, clang, gfortran (r-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (r-hub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (r-hub)
* Windows Server 2008, R-oldrel (win-builder)
* Windows Server 2008, R-release (win-builder)

## R CMD check results

`0 errors | 0 warnings | 2 notes`

## Comments

This is a resubmission. 
We received the following notes when we tested the package:

```
N  Maintainer: 'Alexandra Sarafoglou <alexandra.sarafoglou@gmail.com>'
  
  New submission
  
  Possibly mis-spelled words in DESCRIPTION:
    Sarafoglou (19:28)
    al (19:42)
    et (19:39)
  
N  checking for GNU extensions in Makefiles
   GNU make is a SystemRequirements.
```

## Changes to initial submission (January 21st, 2021)

- Bug fixes: Bayes factors BFr0 and BF0r are displayed correctly
- Enhanced the display output in function bayes_factor corrected for cases in which users test the restricted against the null hypothesis 

- replaced print()/cat() functions with message() in restriction_list() function
- in summary() methods we made the output visible (it previously returned an invisible object)
- added S3 method print.summary.bmult_bridge()
- we reset to user's options(), working directory or par() after changing them (e.g., in the plot() function)
