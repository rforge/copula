## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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

debye1 <- function(x) {
    d <- debye_1(abs(x))
    ## ifelse(x >= 0, d, d - x / 2) ## k = 1, Frees & Valdez 1998, p.9
    d - (x<0) * x / 2
}


debye2 <- function(x) {
    d <- debye_2(abs(x))
    ## ifelse(x >= 0, d, d - x * 2/3) ## k = 2, Frees & Valdez 1998, p.9
    d - (x<0) * x * 2/3
}
