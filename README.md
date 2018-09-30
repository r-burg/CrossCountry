# CrossCountry
Compute the distance between two gliding turning points, using both the Haversine method and Vincenty algorithm

Copyright (C) Roger Burghall 2018
Released under GPL v3.0

There are (at least) two methods of calculating the distance between two points of given latitude and longitude: the rough and ready "Haversine method" and the more difficult but accurate "Vincenty Algorithm".

In order to compare them, and to determine if the simpler method is accurate enough, this software has been written and tested. (In general it seems as though it underestimates distances a little, at least in the UK, and thus can be used for planning. If used to evaluate claims it would occasionally reject valid ones.)

The file "TurnPoints.dat" contains 2016 UK turning points.
