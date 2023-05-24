#include <iostream>
#include <iomanip>
#include "utm.h"
#include <cmath>
#include "geo_ned.h"
#include <math.h>

using namespace utmconv;

void utmconv::tm_set_params(const double &a, const double &f, const double &origin_latitude, const double &central_meridian, const double &false_easting, const int &false_northing, const double &scale_factor)
{
    a_tm = a;
    f_tm = f;
    central_meridian_tm = central_meridian;
    origin_lat_tm = origin_latitude;
    false_e_tm = false_easting;
    false_n_tm = false_northing;
    scale_tm = scale_factor;

    origin_lon_tm = central_meridian;
    es_tm = 2*f_tm - pow(f_tm,2); //eccentricity squared
    ebs_tm = es_tm / (1.0 - es_tm); // second eccentricity squared
    b_tm = a_tm * (1.0 - f_tm);

    double tn = (a_tm - b_tm) / (a_tm + b_tm);
    double tn2 = pow(tn,2);
    double tn3 = pow(tn,3);
    double tn4 = pow(tn,4);
    double tn5 = pow(tn,5);

    ap_tm = a_tm * (1.0 - tn + 5.0 * (tn2 - tn3)/4.0 + 81.0 * (tn4 - tn5)/64.0);
    bp_tm = 3.0 * a_tm * (tn -tn2 + 7.0 * (tn3 - tn4)/8.0 + 55.0 * tn5/64.0)/2.0;
    cp_tm = 15.0 * a_tm * (tn2 - tn3 + 3.0 * (tn4 - tn5)/4.0)/16.0;
    dp_tm = 35.0 * a_tm * (tn3 - tn4 + 11.0 * tn5 / 16.0)/ 48.0;
    ep_tm = 315.0 * a_tm * (tn4 - tn5) / 512.0;

}

double utmconv::sphsn (const double &lat)
{
    return a_tm/sqrt(1.0 - es_tm*pow(sin(lat),2));
}


double utmconv::sphtmd (const double &lat)
{
    return ap_tm*lat - bp_tm*sin(2.0*lat) + cp_tm*sin(4.0*lat) - dp_tm*sin(6.0*lat)+ ep_tm*sin(8.0*lat);
}


double utmconv::sphsr (const double &lat)
{
    return a_tm*(1.0 - es_tm)/pow(denom(lat), 3);
}


double utmconv::denom (const double &lat)
{
    return sqrt(1.0 - es_tm*pow(sin(lat),2));
}


void utmconv::tm_geodetic_to_tranmerc(const double &lat, const double &lon, utm_coords &coords_utm)
{
    double dlam = lon - origin_lon_tm;

    if (dlam > M_PI) dlam -= 2.0*M_PI;
    if (dlam < -M_PI) dlam += 2.0*M_PI;
    if (fabs(dlam) < 2.e-10) dlam = 0.0;

    double s = sin(lat);
    double c = cos(lat);
    double c2 = pow(c,2);
    double c3 = pow(c,3);
    double c5 = pow(c,5);
    double c7 = pow(c,7);
    double t = tan(lat);
    double tan2 = pow(t, 2);
//    double tan3 = pow(t,3);
    double tan4 = pow(t,4);
//    double tan5 = pow(t,5);
    double tan6 = pow(t,6);
    double eta = ebs_tm * c2;
    double eta2 = pow(eta,2);
    double eta3 = pow(eta,3);
    double eta4 = pow(eta,4);

    double sn = sphsn(lat);

    double tmd = sphtmd(lat);
    double tmdo = sphtmd(origin_lat_tm);

    // northing
    double t1 = (tmd - tmdo) * scale_tm;
    double t2 = sn * s * c * scale_tm/2.0;
    double t3 = sn * s * c3 * scale_tm * (5.0 - tan2 + 9.0 * eta + 4.0 * eta2)/24.0;
    double t4 = sn * s * c5 * scale_tm * (61.0 - 58.0 * tan2 + tan4 + 270.0 * eta - 330.0 * tan2 * eta + 445.0 * eta2 + 324.0 * eta3 -680.0 * tan2 * eta2 + 88.0 * eta4 -600.0 * tan2 * eta3 - 192.0 * tan2 * eta4) / 720.0;
    double t5 = sn * s * c7 * scale_tm * (1385.0 - 3111.0 * tan2 + 543.0 * tan4 - tan6) / 40320.0;
    double northing = false_n_tm + t1 + pow(dlam,2.0) * t2 + pow(dlam,4.0) * t3 + pow(dlam,6.0) * t4 + pow(dlam,8.0) * t5;

    // easting
    double t6 = sn * c * scale_tm;
    double t7 = sn * c3 * scale_tm * (1.0 - tan2 + eta ) /6.0;
    double t8 = sn * c5 * scale_tm * (5.0 - 18.0 * tan2 + tan4 + 14.0 * eta - 58.0 * tan2 * eta + 13.0 * eta2 + 4.0 * eta3 - 64.0 * tan2 * eta2 - 24.0 * tan2 * eta3 )/ 120.0;
    double t9 = sn * c7 * scale_tm * ( 61.0 - 479.0 * tan2 + 179.0 * tan4 - tan6 ) /5040.0;
    double easting = false_e_tm + dlam * t6 + pow(dlam,3.0) * t7 + pow(dlam,5.0) * t8 + pow(dlam,7.0) * t9;;

    coords_utm.easting = easting;
    coords_utm.northing = northing;
}


void utmconv::tranmerc_to_geodetic(const double &easting, const double &northing, wgs84_coords &coords_wgs)
{
    double t10 = 0.0;


    double tmdo = sphtmd(origin_lat_tm);
    double tmd = tmdo + (northing - false_n_tm)/scale_tm;
    double sr = sphsr(0.0); // first estimate
    double ftphi = tmd/sr;

    for (int i = 0; i <= 4; i++)
    {
        t10 = sphtmd(ftphi);
        sr = sphsr(ftphi);
        ftphi += (tmd - t10)/sr;
    }

    sr = sphsr(ftphi); // radius of Curvature in the meridian
    double sn = sphsn(ftphi); // radius of Curvature in the meridian

//    double s = sin(ftphi); // sine cosine terms
    double c = cos(ftphi);
    double t = tan(ftphi);

    double tan2 = pow(t,2);
    double tan4 = pow(t,4);
    double eta = ebs_tm * pow(c,2);
    double eta2 = pow(eta,2);
    double eta3 = pow(eta,3);
    double eta4 = pow(eta,4);
    double de = easting - false_e_tm;
    if (fabs(de)<0.0001) de = 0.0;

    // calculate the latitude
    t10 = t / (2.0 * sr * sn * pow(scale_tm,2));
    double t11 = t * (5.0 + 3.0 * tan2 + eta - 4.0 * pow(eta,2) -9.0 * tan2 * eta) / (24.0 * sr * pow(sn,3) * pow(scale_tm, 4));
    double t12 = t * (61.0 + 90.0 * tan2 + 46.0 * eta + 45.0 * tan4 - 252.0 * tan2 * eta - 3.0 * eta2 + 100.0 * eta3 -66.0 * tan2 * eta2 - 90.0 * tan4 * eta + 88.0 * eta4 + 225.0 * tan4 * eta2 + 84.0 * tan2 * eta3 - 192.0 * tan2 * eta4) / (720.0 * sr * pow(sn,5) * pow(scale_tm, 6));
    double t13 = t * (1385.0 + 3633.0 * tan2 + 4095.0 * tan4 + 1575.0 * pow(t,6)) / (40320.0 * sr * pow(sn,7) * pow(scale_tm, 8));
    double lat = ftphi -pow(de, 2) * t10 + pow(de,4) * t11 - pow(de,6) * t12 + pow(de,8)*t13;

    // calculate the longitude
    double t14 = 1.0 / (sn * c * scale_tm);
    double t15 = (1.0 + 2.0 * tan2 + eta) / (6.0 * pow(sn,3) * c * pow(scale_tm, 3));
    double t16 = (5.0 + 6.0 * eta + 28.0 * tan2 - 3.0 * eta2 + 8.0 * tan2 * eta + 24.0 * tan4 - 4.0 * eta3 + 4.0 * tan2 * eta2 + 24.0 * tan2 * eta3) / (120.0 * pow(sn,5) * c * pow(scale_tm,5));
    double t17 = (61.0 +  662.0 * tan2 + 1320.0 * tan4 + 720.0 * pow(t,6)) / (5040.0 * pow(sn,7) * c * pow(scale,7));
    double dlam = de * t14 - pow(de,3) * t15 + pow(de,5) * t16 - pow(de,7) * t17; // difference in longitude
    double lon = origin_lon_tm + dlam;

    while (lat> M_PI/2.0)
    {
        lat = M_PI - lat;
        lon += M_PI;
        if (lon > M_PI) lon -= 2.0 * M_PI;
    }

    while (lat < -M_PI/2.0)
    {
        lat = -(lat+M_PI);
        lon += M_PI;
        if (lon > M_PI) lon -= 2*M_PI;

    }

    if (lon > 2.0*M_PI) lon -= 2.0*M_PI;
    if (lon < -M_PI) lon += 2.0 * M_PI;

    coords_wgs.latitude = lat;
    coords_wgs.longitude = lon;

}

void utmconv::geodetic_to_utm(const double &latitude, const double &longitude, utm_coords &coords){

    int zone = 0;
    int false_northing = 0;
    char hemisphere;
    char zlet;

    double lat = latitude*deg_to_rad;
    double lon = longitude*deg_to_rad;
//    double lat = deg2rad(latitude);
//    double lon = deg2rad(longitude);

    int lat_deg_int = static_cast<int>(latitude);
    int lon_deg_int = static_cast<int>(longitude);

    if (zone_override>0)
        zone = zone_override;

    else
    {
        zone = static_cast<int>((longitude + 180.0)/ 6.0 + 1.0);

        // Handle areas with special conventions (Denmark & South West Norway)
        if (lat_deg_int > 55 && lat_deg_int < 64 && lon_deg_int > -1 && lon_deg_int < 3)
            zone = 31;

        if (lat_deg_int > 55 && lat_deg_int < 64 && lon_deg_int > 2 && lon_deg_int < 12)
            zone = 32;

        // Handle areas with special conventions (Svalbard)
        if (lat_deg_int > 71)
        {
            if (lon_deg_int > -1 && lon_deg_int < 9)
                zone = 31;

            if (lon_deg_int > 8 && lon_deg_int < 21)
                zone = 33;

            if (lon_deg_int > 20 && lon_deg_int < 33)
                zone = 35;

            if (lon_deg_int > 32 && lon_deg_int < 42)
                zone = 37;

        }
    }

    // Calculate central_meridian for this zone
    double central_meridian = ((zone - 1)*6.0 - 180.0 + 3.0) * deg_to_rad;

    // Set falso northing based on hemisphere
    if (latitude >= 0.0)
    {
        false_northing = 0;
        hemisphere = 'N';

        // Determine to UTM zone letter
        if(latitude>=72.0) zlet = 'X';

        else if(latitude>=64.0) zlet = 'W';
        else if(latitude>=56.0) zlet = 'V';
        else if(latitude>=48.0) zlet = 'U';
        else if(latitude>=40.0) zlet = 'T';
        else if(latitude>=32.0) zlet = 'S';
        else if(latitude>=24.0) zlet = 'R';
        else if(latitude>=16.0) zlet = 'Q';
        else if(latitude>=8.0) zlet = 'P';

        else zlet = 'N';
    }

    else
    {
        false_northing = 10000000;
        hemisphere = 'S';

        // Determine to UTM zone letter
        if(latitude>= -8.0) zlet = 'M';

        else if(latitude>= -16.0) zlet = 'L';
        else if(latitude>= -24.0) zlet = 'K';
        else if(latitude>= -32.0) zlet = 'J';
        else if(latitude>= -40.0) zlet = 'H';
        else if(latitude>= -48.0) zlet = 'G';
        else if(latitude>= -56.0) zlet = 'F';
        else if(latitude>= -64.0) zlet = 'E';
        else if(latitude>= -72.0) zlet = 'D';

        else zlet = 'C';

    }

    // set parameters for WGS-84, UTM, the false northing and the zone central meridian
    tm_set_params(wgs84_a, wgs84_f, utm_origin_latitude, central_meridian, utm_false_easting, false_northing, utm_scale_factor);

    // perform conversion
    tm_geodetic_to_tranmerc(lat, lon, coords);

    coords.hemisphere = hemisphere;
    coords.zone = zone;
    coords.zone_letter = zlet;

}

void utmconv::utm_to_geodetic(const char &hemisphere, const int &zone, const double &easting, const double &northing, wgs84_coords &coords)
{
    int false_northing = 0;

    //calculate the central meridian for the zone
    double central_meridian = ((zone - 1)*6.0 - 180.0 + 3.0)*deg_to_rad;

    // determine the false northing based on the hemisphere
    if (hemisphere == 'N') false_northing = 0;
    else false_northing = 10000000;

    // set parameters for WGS-84, UTM, the false northing and the zone central meridian
    tm_set_params(wgs84_a, wgs84_f, utm_origin_latitude, central_meridian, utm_false_easting, false_northing, utm_scale_factor);

    // perform conversion
    tranmerc_to_geodetic(easting, northing, coords);

    coords.latitude *= rad_to_deg;
    coords.longitude *= rad_to_deg;

}