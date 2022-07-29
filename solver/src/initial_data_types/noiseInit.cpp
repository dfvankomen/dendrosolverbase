// NOTE: before, it was calling rand() * 2 - 1, which actually gave enormous
// numbers, not floats I assumed that these were supposed to be random float
// calls between 0 and 2 - 1 so that we get numbers uniformly between -1 and 1

double amplitude = dsolve::BSSN_NOISE_AMPLITUDE;                                      // 10^(-6)

var[VAR::U_ALPHA] = 1.0 + amplitude * dsolve::rand_uniform_mid_range(0, 1);

var[VAR::U_CHI] = 1.0 + amplitude * dsolve::rand_uniform_mid_range(0, 1);

var[VAR::U_K] = amplitude * dsolve::rand_uniform_mid_range(0, 1);

var[VAR::U_GT0] = amplitude * dsolve::rand_uniform_mid_range(0, 1);
var[VAR::U_GT1] = amplitude * dsolve::rand_uniform_mid_range(0, 1);
var[VAR::U_GT2] = amplitude * dsolve::rand_uniform_mid_range(0, 1);

var[VAR::U_BETA0] = amplitude * dsolve::rand_uniform_mid_range(0, 1);
var[VAR::U_BETA1] = amplitude * dsolve::rand_uniform_mid_range(0, 1);
var[VAR::U_BETA2] = amplitude * dsolve::rand_uniform_mid_range(0, 1);

var[VAR::U_B0] = amplitude * dsolve::rand_uniform_mid_range(0, 1);
var[VAR::U_B1] = amplitude * dsolve::rand_uniform_mid_range(0, 1);
var[VAR::U_B2] = amplitude * dsolve::rand_uniform_mid_range(0, 1);

var[VAR::U_GT00] = 1.0 + amplitude * dsolve::rand_uniform_mid_range(0, 1);  // XX
var[VAR::U_GT01] = amplitude * dsolve::rand_uniform_mid_range(0, 1);        // XY
var[VAR::U_GT02] = amplitude * dsolve::rand_uniform_mid_range(0, 1);        // XZ
var[VAR::U_GT11] = 1.0 + amplitude * dsolve::rand_uniform_mid_range(0, 1);  // YY
var[VAR::U_GT12] = amplitude * dsolve::rand_uniform_mid_range(0, 1);        // YZ
var[VAR::U_GT22] = 1.0 + amplitude * dsolve::rand_uniform_mid_range(0, 1);  // ZZ

var[VAR::U_AT00] = amplitude * dsolve::rand_uniform_mid_range(0, 1);  // XX
var[VAR::U_AT01] = amplitude * dsolve::rand_uniform_mid_range(0, 1);  // XY
var[VAR::U_AT02] = amplitude * dsolve::rand_uniform_mid_range(0, 1);  // XZ
var[VAR::U_AT11] = amplitude * dsolve::rand_uniform_mid_range(0, 1);  // YY
var[VAR::U_AT12] = amplitude * dsolve::rand_uniform_mid_range(0, 1);  // YZ
var[VAR::U_AT22] = amplitude * dsolve::rand_uniform_mid_range(0, 1);  // ZZ
