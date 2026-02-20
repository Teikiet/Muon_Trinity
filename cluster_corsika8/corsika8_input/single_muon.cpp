using namespace corsika;

// 1. Define the particle type (Muon)
auto const type = Code::MuMinus;

// 2. Define the energy (1 TeV)
auto const energy = 1_TeV;

// 3. Define the starting position (5 km height)
// Note: Heights are typically relative to the Earth's center or a starting point.
// If using the standard coordinate system where Sea Level is ~6371 km:
Point const injectionPos(rootCS, 0_m, 0_m, 6371_km + 5_km);

// 4. Define the direction (downwards)
DirectionVector const direction(rootCS, {0, 0, -1});
