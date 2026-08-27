// For cartesian coordinates, the ring average is not needed, but we still need to define the functions to avoid compilation errors. The functions will be empty or return default values.

// Ring average related functions in curvilinear coordinates
void gather_Etackle_RingAverage_arrays(){}
void broadcast_Etackle_RingAverage_arrays(){}
int RingAverage_gather_data_3d_array(double ***field_array_scr, double ***field_array_dst){return 0;}
int RingAverage_broadcast_data_3d_array(double ***field_array_scr, double ***field_array_dst){return 0;}

// Special operation in pole regions for curvilinear coordinates
void tackle_Efield_pole(){}
void tackle_Magfield_pole(){}

// Set configuration for ring average in curvilinear coordinates
void init_RingAverage_variables(){}
void set_RingAverage_geo_data(){}
// Spherical: ring average in k direction. The chunk number depends on j direction
void set_RingAverage_config(){}

void gather_value_RingAverage_arrays(){}
void HydroRingAverage(){}
void MagneticRingAverage(){}
void broadcast_value_RingAverage_arrays(){}

