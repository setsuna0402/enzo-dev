#define MAX_SOURCES 3000
// Photons
int    NumberOfPhotonPackages;
//int    NumberOfRenderingPackages;

//  start pointer for linked list of packages
PhotonPackageEntry *PhotonPackages;

// linked list of packages with its work already finished
PhotonPackageEntry *FinishedPhotonPackages;  

// linked list of packages that are "paused", waiting to be merged
PhotonPackageEntry *PausedPhotonPackages;

// linked list of packages used for projections or volume renderings
//PhotonPackageEntry *RenderingPhotonPackages;

// Sources
//int    NumberOfRadiationSources;
grid  **SubgridMarker; // pointers to first array elements of subgrids for each cell
int HasRadiation;

// For adaptive timestep control while restricting the change in HII
// to 10%, we need to record the minimum kph for which a ray passes
// through with a cumulative optical depth >0.1.
float MaximumkphIfront;
// For adaptive timestep control while restricting the change in HeII
// to 10%, we need to record the minimum kph for which a ray passes
// through with a cumulative optical depth >0.1.
float MaximumkpHeIIfront;

int IndexOfMaximumkph;
int IndexOfMaximumkpHeII;
int OriginalProcessorNumber;
