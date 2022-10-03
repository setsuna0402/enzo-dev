struct PhotonBuffer {
  float		Photons;
  int		Type;
  float		Energy;                   
  float		ColumnDensity;
  double	CrossSection;
  FLOAT		EmissionTimeInterval;
  FLOAT		EmissionTime;
  FLOAT		CurrentTime;
  FLOAT		Radius;
  long		ipix;
  int		level;
  FLOAT		SourcePosition[3];
  float		SourcePositionDiff;
  // KH 2022/10/2
  float     SourceCreationTime;     // When Radiation source is formed in code units
  int		SuperSourceID;
};

struct GroupPhotonList {
  char PausedPhoton;
  int ToLevel;
  int ToGrid;
  //  int FromLevel;
  //  int FromGrid;
  PhotonBuffer buffer;
};

