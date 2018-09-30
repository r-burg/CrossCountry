
#ifndef VINCENTYH

#define NDEBUG true

class Position {
  double longitude, latitude;
public:
  bool SetLong(double long1) { longitude = long1; };
  bool SetLat(double lat1) { latitude = lat1; };
  double GetLong( ) { return(longitude); };
  double GetLat( ) { return(latitude); };
};

#define VINCENTYH

#endif


