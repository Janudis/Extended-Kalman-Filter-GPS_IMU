#ifndef UTM_HPP_
#define UTM_HPP_

namespace utmconv {
// Define global variables
    static double deg_to_rad = 0.01745329251994329575f;
    static double rad_to_deg = 57.29577951308232087721f;

    static double wgs84_a = 6378137.0; // WGS84 semi-major axis of ellipsoid [m]
    static double wgs84_f = 1.0 / 298.257223563; // WGS84 flattening of ellipsoid

    static double utm_false_easting = 500000.0;
    static double utm_scale_factor = 0.9996;
    static double utm_origin_latitude = 0.0;

    static double false_e = 500000.00;
    static double false_n = 0.00;
    static double scale = 0.9996;
    static int zone_override = 0;

// Define global variables for Transverse Mercator (tm)
    static double a_tm = 0.0;
    static double f_tm = 0.0;
    static double central_meridian_tm = 0.0;
    static double origin_lat_tm = 0.0;
    static double false_e_tm = 0.0;
    static double false_n_tm = 0.0;
    static double scale_tm = 0.0;

    static double origin_lon_tm = 0.0;
    static double es_tm = 0.0;
    static double ebs_tm = 0.0;
    static double b_tm = 0.0;

    static double ap_tm = 0.0;
    static double bp_tm = 0.0;
    static double cp_tm = 0.0;
    static double dp_tm = 0.0;
    static double ep_tm = 0.0;


    struct utm_coords {
        char hemisphere;
        int zone;
        char zone_letter;
        double easting; // X
        double northing; // Y
    };

    struct wgs84_coords {
        double latitude;
        double longitude;
    };


// Declaration of functions
    void geodetic_to_utm(const double &lat, const double &lon, utm_coords &coords);

    void utm_to_geodetic(const char &hemisphere, const int &zone, const double &easting, const double &northing,
                         wgs84_coords &coords);

    void tm_set_params(const double &a, const double &f, const double &origin_latitude, const double &central_meridian,
                       const double &false_easting, const int &false_northing, const double &scale_factor);

    void tm_geodetic_to_tranmerc(const double &lat, const double &lon, utm_coords &coords_utm);

    void tranmerc_to_geodetic(const double &easting, const double &northing, wgs84_coords &coords_wgs);

    double sphsn(const double &lat);

    double sphtmd(const double &lat);

    double sphsr(const double &lat);

    double denom(const double &lat);

}
#endif
