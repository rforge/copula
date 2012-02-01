## Copyright (C) 2010 Marius Hofert and Martin Maechler
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

####  Generating Exponentially Tilted Stable Random Vars (r ET stable)
####  ================================================================
####  Experiments with retstable*() versions
###
### More (computer intensive) experiments are in
### ../demo/retstable.R
### ~~~~~~~~~~~~~~~~~~~

require(nacopula)
source(system.file("Rsource", "utils.R", package="nacopula"))
##--> tryCatch.W.E(), canGet()

###---- TODO --- *some* good "regression tests"
### using both retstableR() and retstable()
