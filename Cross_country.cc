/***************************************************************

                    Cross_country.cc

    Copyright (c) R. Burghall 2017 all rights reserved

	Issued under GPL v3.0. If you copy and/or distribute
	this program, copy and/or distribute this licence with it.

    Version 0.1 RB 2017-12-21

    Written to investigate the merits of Haversine and 
    Vincenty methods of determining the distance between
    two points on the surface of the Earth.

    The Haversine formula is said to be a crude and inaccurate 
    method while the Vincenty algorithm is complex, difficult
    and highly accurate - except that the Earth isn't really 
    an ellipsoid, WGS84 or otherwise. There are hills!

***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <assert.h>

#include "Vincenty.h"

using namespace std;

#define VINCENTY true
#define HALVERSINE false

#define VERSION 0.1

double toRad(double degrees)
{
	return(degrees * M_PI / 180.0);
}



double square(double x)
{
	return(x * x);
}



double arcsin(double x)
{
	return(2.0 * atan(x / (1.0 + sqrt(1.0 - x*x))));
}



double hav(double theta)
{
	return(sin(theta/2.0) * sin(theta/2.0));
}



double haversine(Position p1, Position p2)
{
	double R = 6371.0;
	double d = 2.0 * R * arcsin(sqrt(hav(p1.GetLat( )-p2.GetLat( ))+cos(p1.GetLat( ))*cos(p2.GetLat( ))*hav(p2.GetLong( )-p1.GetLong( ))));
	return(d);
}



double Vincenty(Position x1, Position x2, double &bearing)
{
	double a, b;                                  /// Major and minor semiaxes of the ellipsoid.
	a = 6378137.0;
	b = 6356752.314245;
	double f = (a - b) / a;                       /// Flattening.
	assert((f > 1.0 / 298.257224) && (f < 1.0 / 298.257223));
	double long1, long2;
	double Phi1 = x1.GetLat( ),
	Phi2 = x2.GetLat( );                         /// Geodetic latitude.
	double L = x2.GetLong( )-x1.GetLong( );       /// Difference in longitude.
	double U1 = atan((1.0 - f) * tan(Phi1));      /// "Reduced latitude" = atan((1 - f) * tan(Phi1)).
	double U2 = atan((1.0 - f) * tan(Phi2));      /// "Reduced latitude" = atan((1 - f) * tan(Phi2)).
	double Lambda = L;                                /// = L as first approximation.
	double Lambda1, diff;
	double sigma, sinsigma, cossigma, deltasigma;
	double sinalpha, cossqalpha, cos2sigmam, C;
	double alpha1, alpha2;
	double usq, A, B, s;
	double sinlambda, coslambda;
	double sinU1 = sin(U1), sinU2 = sin(U2), cosU1 = cos(U1), cosU2 = cos(U2);
	double i = 0, n = 100;

	do {
		sinlambda = sin(Lambda), coslambda = cos(Lambda);
		sinsigma = sqrt(square(cosU2 * sinlambda) + square(cosU1 * sinU2 - sinU1 * cosU2 * coslambda));
		if(sinsigma == 0.0) return(0.0);
		cossigma = sinU1 * sinU2 + cosU1 * cosU2 * coslambda;
		sigma = atan2(sinsigma, cossigma);
		sinalpha = cosU1 * cosU2 * sinlambda / sinsigma;
		cossqalpha = 1.0 - sinalpha * sinalpha;
		cos2sigmam = cossigma - 2.0 * sinU1 * sinU2 / cossqalpha;
		C = f / 16.0 * cossqalpha * (4.0 + f *(4.0 - 3.0 * cossqalpha));
		Lambda1 = L + (1.0 - C) * f * sinalpha * (sigma + C * sinsigma * (cos2sigmam + C * cossigma * (-1.0 + 2.0 * cos2sigmam * cos2sigmam)));
		diff = Lambda - Lambda1;
		Lambda = Lambda1;
//    cout << "Lambda = " << Lambda << "\n";
		if(++i > n) break;
	} while(fabs(diff) > 1e-13);
	  
	usq = cossqalpha * (a * a - b * b) / (b * b);
	A = 1.0 + usq / 16384.0 * (4096.0 + usq * (-768.0 + usq * (320.0 -175.0 * usq)));
	B = usq / 1024.0 * (256.0 + usq * (-128.0 + usq * (74.0 - 47.0 * usq)));
	deltasigma = B * sinsigma * (cos2sigmam + B / 4.0 * (cossigma * (-1.0 + 2.0 * cos2sigmam * cos2sigmam) - B / 6.0 * cos2sigmam * (-3.0 + 4.0 * sinsigma * sinsigma) * (-3.0 + 4.0 * cos2sigmam * cos2sigmam)));
	s = b * A * (sigma - deltasigma);
	alpha1 = atan2(cosU2 * sin(Lambda), cosU1 * sinU2 - sinU1 * cosU2 * cos(Lambda));
	alpha2 = atan2(cosU1 * sin(Lambda), -sinU1 * cosU2 + cosU1 * sinU2 * cos(Lambda));
    bearing = (alpha1+alpha2)/2.0;
	return(s);
}



/*!
Convert latitude or longitude string to a double.
*/
double str_to_val(char *ptr)
{
	int n;
	int dd = 0.0;
	double mm = 0.0;
	while(*ptr == ' ') ptr++;
	strcat(ptr, "\n");
	n = sscanf(ptr, "%d%lf", &dd, &mm);
	return(dd + mm / 60.0);
}



void clean(char* line)
{
    static char *ptr;
    for(ptr = line; (*ptr != '\r') && (*ptr != '\n'); ptr++) {
        if(*ptr == '\r') *ptr = '\0';
        if(*ptr == '\n') *ptr = '\0';
    }
    *ptr = '\0';
}



bool cfound(char *ptr, char c)
{
    int i;
    int n = strlen(ptr);
    for(i=0; i < n; i++) if(*ptr++ == c) return(true);
    return(false);
}



double leg(char *a1, char *a2, double &h, double &v)
{
    char c1, c2;
    int d1, d2;
    double m1, m2;
    char name[258], bfr[258];
    char tp1[258], tp2[258];
    bool found = false;
    int l;
    FILE *tpfile;
	Position x1, x2;

    // Open TP file.
    tpfile = fopen("./TurnPoints.dat", "r");
    if (tpfile == NULL) { perror ("Error opening file"); return(-1); }
    else {
        while(!found) {
            // Read a line; it should be the name of the turning point.
            if(fgets(name, 256, tpfile) == NULL) ;
            clean(name);
            // Get the next line, which should be the trigraph.
            if(fgets(bfr, 256, tpfile) == NULL) ;
            clean(bfr);
            l = 1;
            // Have we found it?
            if(strncmp(bfr, a1, 256) == 0) {
//                puts(name);
                found = true;
            }
            // Look for location.
            if(found) for(;;) {
                bfr[0] = '\0';
                if(fgets(bfr, 256, tpfile) == NULL) ;
                clean(bfr);
                l++;
                if(found && (l == 9)) {
                    strncpy(tp1, bfr, 256) ;
                    break;
                }
            }
            // Look for a blank line between TPs
            do {
                if(fgets(bfr, 256, tpfile) == NULL) ;
                clean(bfr);
                if(strlen(bfr) == 0) break;
            } while((bfr[0] != '\r') && (bfr[0] != '\n'));
            if(found) break;           // We found it - don't keep looking
            if(feof(tpfile)) exit(-3); // End of file? What should we do?
        }
        if(found) {  printf("%s: %s\n", name, tp1); }
    }
    fclose(tpfile);

    tpfile = NULL;
    tpfile = fopen("./TurnPoints.dat", "r");
    found = false;
    if (tpfile == NULL) { perror ("Error opening file"); return(-1); }
    else {
        while(!found) {
            // Read a line; it should be the name of the turning point.
            if(fgets(name, 256, tpfile) == NULL) ;
            clean(name);
            // Get the next line, which should be the trigraph.
            if(fgets(bfr, 256, tpfile) == NULL) ;
            clean(bfr);
            l = 1;
            // Have we found it?
            if(strncmp(bfr, a2, 256) == 0) {
//                puts(name);
                found = true;
            }
            // Look for the location
            if(found) for(;;) {
                bfr[0] = '\0';
                if(fgets(bfr, 256, tpfile) == NULL) ;
                clean(bfr);
                l++;
                if(found && (l == 9)) {
                    // This should be the lat & long
                    strncpy(tp2, bfr, 256) ;
                    break;
                }
            }
            // Look for a blank line between TPs
            do {
                if(fgets(bfr, 256, tpfile) == NULL) ;
                clean(bfr);
                if(strlen(bfr) == 0) break;
            } while((bfr[0] != '\r') && (bfr[0] != '\n'));
            if(found) break;        // We found it - don't keep looking
            if(feof(tpfile)) exit(-2); // End of file? What should we do?
        }
        if(found) {  printf("%s: %s\n", name, tp2); }
    }
    fclose(tpfile);

    strncat(tp1, "\n", 256);
    sscanf(tp1, "%d %lf%c %d %lf%c", &d1, &m1, &c1, &d2, &m2, &c2);
    if(cfound(tp1, 'S')) d1 = -d1, m1 = -m1;
    if(cfound(tp1, 'E')) d2 = -d2, m2 = -m2;
    x1.SetLat(d1*M_PI/180.0 + m1*M_PI/(180.0*60.0));
    x1.SetLong(d2*M_PI/180.0 + m2*M_PI/(180.0*60.0));
//    printf("<%d %lf>, ", d1, m1);
//    printf("<%d %lf>\n", d2, m2);

    strncat(tp2, "\n", 256);
    sscanf(tp2, "%d %lf%c %d %lf%c", &d1, &m1, &c1, &d2, &m2, &c2);
    if(cfound(tp2, 'S')) d1 = -d1, m1 = -m1;
    if(cfound(tp2, 'E')) d2 = -d2, m2 = -m2;
    x2.SetLat(d1*M_PI/180.0 + m1*M_PI/(180.0*60.0));
    x2.SetLong(d2*M_PI/180.0 + m2*M_PI/(180.0*60.0));
//    printf("<%d %lf>, ", d1, m1);
//    printf("<%d %lf>\n", d2, m2);

    double hdg;
    double s = Vincenty(x1, x2, hdg);
    printf("Vincenty --> %0.3f km\n", s/1000.0);
    hdg = -hdg*180.0/M_PI;
    if(hdg < 0.0) hdg = 360.0 + hdg;
    printf("Track: %0.1lf\n", hdg);
    v = s / 1000.0;

    s = haversine(x1, x2);
#if HALVERSINE
    printf("Haversine --> %0.3f km\n\n", s);
#endif
    h = s;
}



void upper_case(char *ptr)
{
    int i;
    while(*ptr) {
        *ptr = toupper(*ptr);
        ++ptr;
    }
}



/*!
Try 50 3.979, 5 42.885
58 38.64133, 3 4.205667
*/
int main(int argc, char **argv)
{
    char c1, c2;
    char a1[256], a2[256];
    double h, v, htotal = 0.0, vtotal = 0.0;

    puts("Enter 'END' as t.p. to end cross-country.\n\n");
    if(argc == 1) {
        // No parameters after the program name, so enter TPs manually.
        puts("Start TP: ");
        scanf("%s", a1);
        upper_case(a1);
        while(1) {
            puts("Next TP: ");
            a2[0] = '\0';
            scanf("%s", a2);
            upper_case(a2);
            if(strcmp(a2, "END") == 0) break;
            leg(a1, a2, h, v);
            htotal += h;
            vtotal += v;
//            leg(a1, a2);
            strcpy(a1, a2);
        }
#if HALVERSINE
        printf("Total distance (Haversine) = %0.3f; ", htotal);
#endif
        printf("Total distance (Vincenty) = %0.3f", vtotal);
    } else if(argc == 3) {
        strncpy(a1, argv[1], 5);
        strncpy(a2, argv[2], 5);
        leg(a1, a2, htotal, vtotal);
    } else {
        printf("*** Too many arguments!");
        exit(-1);
    }

    printf("\nBye!\n");

/*    getchar( );         // I seem to need to remove a character already waiting
    // Now wait to be dismissed
    do {
        c1 = getchar( );
    } while((c1 != ' ') && (c1 != '\r') && (c1 != '\n'));
*/
}




